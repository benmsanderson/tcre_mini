"""Joint land process model (A3): beta_L, gamma_L tropics and gamma_L extratropics in one PyMC model.

Re-implements upstream build_beta_land_model and build_gamma_land_model (patched copy in
upstream_patched) with one change: the per-model biases bias_m of the three components are correlated
through an LKJ correlation matrix, instead of independent. The common biases stay independent, as
upstream (their correlation is not identified from one ensemble; a first version that also correlated
them did not converge: r-hat 1.58, 75 divergences).
Priors, process structure, the IAV emergent constraint on gamma_L tropics, the nitrogen
pseudo-observations on eta and the definition of the reported truth
(X_true = baseline * prod(eta) + sum(delta), excluding biases) are as upstream.
"""
from __future__ import annotations

from copy import deepcopy

import numpy as np
import pymc as pm
import pytensor.tensor as pt
import xarray as xr

COMPS = ("beta_land", "gamma_land_tropics", "gamma_land_extratropics")


def _baseline(spec, name):
    if spec["dist"] == "LogNormal":
        return pm.LogNormal(name, spec["mu"], spec["sigma"])
    return pm.Normal(name, spec["mu"], spec["sigma"])


def _delta(spec, name):
    if spec["dist"] == "Normal":
        return pm.Normal(name, spec["mu"], spec["sigma"])
    if spec["dist"] == "LogNormal":
        return pm.LogNormal(name, spec["mu"], spec["sigma"])
    if spec["dist"] == "NegLogNormal":
        return pm.Deterministic(name, -pm.LogNormal(f"negative_{name}", spec["mu"], spec["sigma"]))
    raise ValueError(spec["dist"])


def build_joint_land_model(beta_land: xr.DataArray, gamma_t: xr.DataArray, gamma_x: xr.DataArray,
                           lookup_table: xr.DataArray, iav_df, iav_obs_mean, iav_obs_sigma,
                           pri_beta: dict, pri_gamma: dict, lkj_eta: float = 2.0) -> pm.Model:
    """beta_L and gamma_L (tropics, extratropics) fitted together, with per-model biases correlated across
    the three through an LKJ prior; otherwise as the development code's separate builders."""
    models = list(beta_land.model.values)
    assert list(gamma_t.model.values) == models == list(gamma_x.model.values)
    lt = lookup_table.sel(model=models)
    M = len(models)
    Y = np.column_stack([beta_land.values, gamma_t.values, gamma_x.values])
    add_p = list(lt.process.values)
    mult_p = list(pri_beta["eta"].keys())
    L_add = np.asarray(lt.sel(process=add_p)).astype(float)
    L_mult = np.asarray(lt.sel(process=mult_p)).astype(float)
    pri_t = deepcopy(pri_gamma)
    pri_t["delta"]["permafrost"] = {"dist": "Normal", "mu": 0.0, "sigma": 1e-6}
    with pm.Model(coords={"model": models, "comp": list(COMPS)}) as model:
        # shared correlation of biases across the three components
        r = pm.LKJCorr("bias_r", n=3, eta=lkj_eta)  # packed upper triangle (r01, r02, r12)
        R = pm.Deterministic("bias_R", pt.stack([pt.stack([1.0, r[0], r[1]]), pt.stack([r[0], 1.0, r[2]]),
                                                  pt.stack([r[1], r[2], 1.0])]))
        LR = pt.slinalg.cholesky(R)
        sig_struct = pt.stack([
            pm.HalfNormal("beta_land_sigma_struct", pri_beta["sigma_struct"]["sigma"]),
            pm.HalfNormal("gamma_land_tropics_sigma_struct", pri_t["sigma_struct"]["sigma"]),
            pm.HalfNormal("gamma_land_extratropics_sigma_struct", pri_gamma["sigma_struct"]["sigma"])])
        rho = pt.stack([pm.Beta(f"{c}_rho", 5, 5) for c in COMPS])
        sig_m, sig_c = sig_struct * pt.sqrt(1 - rho), sig_struct * pt.sqrt(rho)
        z_m = pm.Normal("bias_m_raw", 0, 1, shape=(M, 3))
        z_c = pm.Normal("bias_common_raw", 0, 1, shape=3)
        bias_m = pm.Deterministic("bias_m", (z_m @ LR.T) * sig_m, dims=("model", "comp"))
        bias_c = pm.Deterministic("bias_common", z_c * sig_c, dims="comp")

        # beta_L: multiplicative processes (as upstream)
        b0 = _baseline(pri_beta["baseline"], "beta_land_unscaled")
        eta_mu = np.array([pri_beta["eta"][p]["mu"] for p in mult_p])
        eta_sd = np.array([pri_beta["eta"][p]["sigma"] for p in mult_p])
        log_eta = pm.Deterministic("beta_land_log_eta", eta_mu + pm.Normal("beta_land_log_eta_offset", 0, 1, shape=len(mult_p)) * eta_sd)
        eta = pm.Deterministic("beta_land_eta", pt.exp(log_eta))
        s_mult = np.array([pri_beta["sigma_mult_struct"].get(p, 0.05) for p in mult_p])
        eps_mult = pm.Normal("beta_land_eps_mult", 0, s_mult[:, None], shape=(len(mult_p), M))
        mu_b = b0 * pt.exp(pt.sum(L_mult * (log_eta[:, None] + eps_mult), axis=0))
        pm.Normal("beta_land_eta_nitrogen_1", mu=eta[mult_p.index("nitrogen")],
                  sigma=pri_beta["nitrogen_obs"]["sigma"], observed=pri_beta["nitrogen_obs"]["values"])
        pm.Deterministic("beta_land", b0 * pt.prod(eta))

        mus = [mu_b]
        for c, pri in (("gamma_land_tropics", pri_t), ("gamma_land_extratropics", pri_gamma)):
            g0 = _baseline(pri["baseline"], f"{c}_unscaled")
            d = pt.stack([_delta(pri["delta"][p], f"{c}_{p}") for p in add_p])
            s_add = np.array([pri["sigma_process_struct"].get(p, 5.0) for p in add_p])
            eps_add = pm.Normal(f"{c}_eps_add", 0, s_add[:, None], shape=(len(add_p), M))
            mus.append(g0 + pt.sum(L_add * (d[:, None] + eps_add), axis=0))
            pm.Deterministic(c, g0 + pt.sum(d))
        mu = pt.stack(mus, axis=1) + bias_m + bias_c[None, :]
        noise = np.array([pri_beta["likelihood_noise"], pri_t["likelihood_noise"], pri_gamma["likelihood_noise"]])
        pm.Normal("land_CMIP", mu=mu, sigma=noise[None, :], observed=Y)

        # IAV emergent constraint on gamma_L tropics (as upstream add_emergent_constraint)
        iav = iav_df.set_index("Model").reindex(models)
        mask = iav.gamma_IAV.notna().values
        m_ec, b_ec = pm.Normal("m_IAV", 0, 5), pm.Normal("b_IAV", 0, 5)
        s_ec = pm.HalfNormal("regression_sigma_IAV", 1.0)
        pm.Normal("IAV_CMIP_lik", mu=m_ec * Y[mask, 1] + b_ec,
                  sigma=pt.sqrt(s_ec ** 2 + iav.sigma_iav.values[mask] ** 2), observed=iav.gamma_IAV.values[mask])
        pm.Normal("IAV_obs_lik", mu=m_ec * model["gamma_land_tropics"] + b_ec, sigma=iav_obs_sigma, observed=iav_obs_mean)
        pm.Deterministic("gamma_land", model["gamma_land_tropics"] + model["gamma_land_extratropics"])
    return model
