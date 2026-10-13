"""Perfect-model (leave-one-out and leave-one-family-out) harness for the upstream framework (A4).

For a held-out model m:
1. remove m (or its whole family) from the land training data (sheet), the IAV regression, the ocean
   training data (Terhaar) and the ocean predictors; drop the paired nitrogen experiments that come from
   m's family (ACCESS-ESM C/CN for the UM-ACCESS family, UKESM1-1 ctrl/nonlim likewise);
2. replace the real-world observations with m's own simulated observables: gamma_IAV (with the real
   observational sigma), AMOC, SSS, CUC (with the real observational sigmas), restandardised on the
   training ensemble; where m has no simulated IAV the constraint is made uninformative (sigma 1e3);
3. refit the three upstream component models unchanged (patched builders);
4. combine sample-wise as upstream, with alpha set to m's own TCR / dCa (so the test isolates the
   carbon-cycle part; the raw-ensemble baseline gets the same alpha).
Two land outputs are scored: the reported truth (baseline * prod(eta) + sum(delta), "all processes on",
what upstream reports for the real world) and the flag-conditional predictive for a new model with m's
own process flags (tests the process regression out of sample).

Families approximate shared components: UM-ACCESS (UM/JULES or CABLE lineage), NCAR (CLM land),
NEMO-IPSL (NEMO ocean), each other model on its own.
"""
from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd
import xarray as xr

from . import upstream as up
from . import upstream_data as ud

FAMILIES = {
    "ACCESS-ESM1-5": "UM-ACCESS", "UKESM1-0-LL": "UM-ACCESS",
    "CESM2": "NCAR", "NorESM2-LM": "NCAR", "CMCC-ESM2": "NCAR",
    "IPSL-CM6A-LR": "NEMO-IPSL", "CNRM-ESM2-1": "NEMO-IPSL", "EC-Earth3-CC": "NEMO-IPSL",
}
NITROGEN_OBS_MODEL = ["ACCESS-ESM1-5", "UKESM1-0-LL"]  # sources of the paired experiments ACCESS-ESM C/CN, UKESM1-1 ctrl/nonlim


def family(m: str) -> str:
    return FAMILIES.get(m, m)


def holdout_set(m: str, mode: str, models) -> list[str]:
    if mode == "loo":
        return [m]
    if mode == "lfo":
        return [x for x in models if family(x) == family(m)]
    raise ValueError(mode)


def build_inputs(m: str, mode: str = "loo") -> dict:
    """Upstream inputs with m (or its family) removed and m's observables as pseudo-observations."""
    base = up.inputs()
    all_models = set(base["gamma_land_data"].model.values) | set(base["beta_ocean"].model.values)
    out = set(holdout_set(m, mode, all_models))
    keep = lambda da: da.sel(model=[x for x in da.model.values if x not in out])
    inp = dict(base)
    inp["gamma_land_data"] = keep(base["gamma_land_data"])
    inp["beta_land_data"] = keep(base["beta_land_data"])
    inp["beta_ocean"], inp["gamma_ocean"] = keep(base["beta_ocean"]), keep(base["gamma_ocean"])
    # IAV pseudo-observation
    iav = base["iav_df"]
    row = iav.loc[iav.Model == m]
    if len(row) and np.isfinite(row.gamma_IAV.values[0]):
        inp["iav_obs_mean"] = float(row.gamma_IAV.values[0])
    else:
        inp["iav_obs_mean"], inp["iav_obs_sigma"] = -5.0, 1e3
    inp["iav_df"] = iav.loc[~iav.Model.isin(out)].reset_index(drop=True)
    # ocean predictors: training ensemble, m's values as pseudo-obs with real obs sigmas
    pred = ud.terhaar_predictors()
    if m not in pred.model.values:
        raise KeyError(f"{m} has no ocean predictors")
    train = pred.sel(model=[x for x in pred.model.values if x not in out])
    assert list(train.model.values) == list(inp["beta_ocean"].model.values)
    pseudo = {v: (float(pred[v].sel(model=m)), ud.OCEAN_EC_OBS[v][1]) for v in ("AMOC", "SSS", "CUC")}
    inp["ec"] = ud.standardize_ocean_ec(train, pseudo)
    # paired nitrogen experiments from the held-out family
    pri_b = deepcopy(up.priors.beta_land)
    drop = (lambda src: src in out) if mode == "loo" else (lambda src: family(src) in {family(x) for x in out})
    vals = [v for v, src in zip(pri_b["nitrogen_obs"]["values"], NITROGEN_OBS_MODEL) if not drop(src)]
    if vals:
        pri_b["nitrogen_obs"]["values"] = vals
    else:
        pri_b["nitrogen_obs"] = {"values": [1.0], "sigma": 1e3}
    inp["priors_beta_land"] = pri_b
    inp["heldout"] = sorted(out)
    return inp


def _post(idata, v):
    x = idata.posterior[v].values
    return x.reshape(x.shape[0] * x.shape[1], *x.shape[2:])


def flag_predictive(m: str, gl, bl, pri_g: dict, pri_b: dict, seed: int = 0) -> dict:
    """Posterior predictive of a new model with m's process flags (upstream generative structure):
    X_m = baseline * exp(sum_p L_mp (log eta_p + e_p)) + sum_p L_mp (delta_p + e'_p) + bias_common + b_m + noise,
    with fresh model-level effects e, e' ~ N(0, sigma_p), b_m ~ N(0, sigma_m), noise ~ N(0, likelihood_noise)."""
    rng = np.random.default_rng(seed)
    flags = ud.lookup_table().sel(model=m)
    out = {}
    # beta_land
    procs = list(pri_b["eta"].keys())
    L = np.array([float(flags.sel(process=p)) for p in procs])
    log_eta = _post(bl.idata, "beta_land_log_eta")
    n = log_eta.shape[0]
    e = rng.normal(0, [pri_b["sigma_mult_struct"].get(p, 0.05) for p in procs], (n, len(procs)))
    out["beta_land"] = (_post(bl.idata, "beta_land_unscaled") * np.exp(((log_eta + e) * L).sum(1))
                        + _post(bl.idata, "beta_land_bias_common") + rng.normal(0, 1, n) * _post(bl.idata, "beta_land_sigma_m")
                        + rng.normal(0, pri_b["likelihood_noise"], n))
    # gamma_land tropics + extratropics
    g = 0.0
    procs_a = list(ud.lookup_table().process.values)
    La = np.array([float(flags.sel(process=p)) for p in procs_a])
    for c in ("gamma_land_tropics", "gamma_land_extratropics"):
        d = _post(gl.idata, f"{c}_delta")
        e = rng.normal(0, [pri_g["sigma_process_struct"].get(p, 5.0) for p in procs_a], d.shape)
        g = g + (_post(gl.idata, f"{c}_unscaled") + ((d + e) * La).sum(1) + _post(gl.idata, f"{c}_bias_common")
                 + rng.normal(0, 1, d.shape[0]) * _post(gl.idata, f"{c}_sigma_m") + rng.normal(0, pri_g["likelihood_noise"], d.shape[0]))
    out["gamma_land"] = g
    return out


def fit_heldout(m: str, mode: str, alpha: float, seed: int = 0, sampler=None) -> dict:
    inp = build_inputs(m, mode)
    gl = up.fit_gamma_land(inp, seed=seed + 1, sampler=sampler)
    bl = up.fit_beta_land(inp, pri=inp["priors_beta_land"], seed=seed + 2, sampler=sampler)
    oc = up.fit_ocean(inp, seed=seed + 3, sampler=sampler)
    b = up.flat(bl.idata, "beta_land") + up.flat(oc.idata, "beta_ocean")
    g = up.flat(gl.idata, "gamma_land") + up.flat(oc.idata, "gamma_ocean")
    fp = flag_predictive(m, gl, bl, up.priors.gamma_land, inp["priors_beta_land"], seed=seed + 9)
    bf = fp["beta_land"] + up.flat(oc.idata, "beta_ocean")
    gf = fp["gamma_land"] + up.flat(oc.idata, "gamma_ocean")
    return dict(
        beta_land_flag=fp["beta_land"], gamma_land_flag=fp["gamma_land"],
        af_flag=ud.K / (ud.K + bf + alpha * gf), tcre_flag=1e3 * alpha / (ud.K + bf + alpha * gf),
        model=m, mode=mode, heldout=inp["heldout"],
        beta_land=up.flat(bl.idata, "beta_land"), gamma_land=up.flat(gl.idata, "gamma_land"),
        beta_ocean=up.flat(oc.idata, "beta_ocean"), gamma_ocean=up.flat(oc.idata, "gamma_ocean"),
        af=ud.K / (ud.K + b + alpha * g), tcre=1e3 * alpha / (ud.K + b + alpha * g),
        diagnostics={k: f.diagnostics for k, f in dict(gamma_land=gl, beta_land=bl, ocean=oc).items()},
        n_train_land=int(inp["gamma_land_data"].sizes["model"]), n_train_ocean=int(inp["beta_ocean"].sizes["model"]),
    )
