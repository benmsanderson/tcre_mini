"""Run the upstream TCRE framework (notebooks/TCRE_total.ipynb) as a function.

Uses the minimally patched builders in tcre_evidence.upstream_patched (the pinned
commit does not construct; see that package's docstring) and the frozen data
snapshot. Each component is sampled separately and combined sample-wise, exactly
as in notebook cells 15-36.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field

import arviz as az
import numpy as np
import pymc as pm
import xarray as xr

from . import upstream_data as ud
from .upstream_patched import model_builders
from .upstream_priors import priors

SAMPLER_DEFAULTS = dict(draws=1000, tune=1000, chains=4, target_accept=0.99)


@dataclass
class FitResult:
    """One fitted component: its name, the ArviZ InferenceData and sampler diagnostics."""
    name: str
    idata: az.InferenceData
    diagnostics: dict = field(default_factory=dict)


def diagnostics(idata, var_names) -> dict:
    """Max r-hat, min bulk/tail ESS over the given variables, and divergence count."""
    post = idata.posterior[var_names]
    rhat = az.rhat(post)
    ess_b = az.ess(post, method="bulk")
    ess_t = az.ess(post, method="tail")
    as_max = lambda d: float(max(float(d[v].max()) for v in d.data_vars))
    as_min = lambda d: float(min(float(d[v].min()) for v in d.data_vars))
    div = int(idata.sample_stats["diverging"].sum()) if "diverging" in idata.sample_stats else -1
    return dict(rhat_max=as_max(rhat), ess_bulk_min=as_min(ess_b), ess_tail_min=as_min(ess_t),
                divergences=div, n_draws=int(post.sizes["chain"] * post.sizes["draw"]))


def _sample(model, seed, sampler, progressbar):
    kw = {**SAMPLER_DEFAULTS, **(sampler or {})}
    with model:
        return pm.sample(random_seed=seed, progressbar=progressbar, **kw)


def inputs(iav_obs=None, alpha_prior=None) -> dict:
    """Assemble the notebook's model inputs from the snapshot."""
    cmip = ud.cmip_sheet_dataset()
    iav_df, iav_mean, iav_sigma = ud.iav()
    if iav_obs is not None:
        iav_mean, iav_sigma = iav_obs
    ocean = ud.terhaar_ocean(69)
    ec = ud.standardize_ocean_ec(ud.terhaar_predictors())
    return dict(
        gamma_land_data=cmip[["gamma_land_tropics", "gamma_land_extratropics", "lookup_table"]].dropna(dim="model"),
        beta_land_data=cmip[["beta_land", "lookup_table"]].dropna(dim="model"),
        iav_df=iav_df, iav_obs_mean=iav_mean, iav_obs_sigma=iav_sigma,
        beta_ocean=ocean.beta_ocean, gamma_ocean=ocean.gamma_ocean, ec=ec,
        alpha=alpha_prior or deepcopy(priors.alpha),
    )


def fit_gamma_land(inp, pri=None, seed=1, sampler=None, progressbar=False) -> FitResult:
    """Land carbon-climate feedback: additive process model for tropics and extratropics, with the IAV constraint."""
    d = inp["gamma_land_data"]
    m = model_builders.build_gamma_land_model(
        d["gamma_land_tropics"], d["gamma_land_extratropics"], d["lookup_table"],
        inp["iav_df"], inp["iav_obs_mean"], inp["iav_obs_sigma"], pri or priors.gamma_land)
    idata = _sample(m, seed, sampler, progressbar)
    return FitResult("gamma_land", idata, diagnostics(idata, ["gamma_land", "gamma_land_tropics", "gamma_land_extratropics"]))


def fit_beta_land(inp, pri=None, seed=2, sampler=None, progressbar=False) -> FitResult:
    """Land carbon-concentration feedback: multiplicative process model with the paired nitrogen runs."""
    d = inp["beta_land_data"]
    m = model_builders.build_beta_land_model(d["beta_land"], d["lookup_table"], pri or priors.beta_land)
    idata = _sample(m, seed, sampler, progressbar)
    return FitResult("beta_land", idata, diagnostics(idata, ["beta_land"]))


def fit_ocean(inp, pri=None, seed=3, sampler=None, progressbar=False) -> FitResult:
    """Ocean beta and gamma: joint emergent constraint on AMOC, sea-surface salinity and chemical uptake capacity."""
    ec = inp["ec"]
    m = model_builders.build_joint_ocean_model(
        inp["beta_ocean"], inp["gamma_ocean"], ec["O_m"], ec["O_obs"], ec["O_obs_unc"], pri or priors.joint_ocean)
    idata = _sample(m, seed, sampler, progressbar)
    return FitResult("ocean", idata, diagnostics(idata, ["beta_ocean", "gamma_ocean"]))


def flat(idata, var) -> np.ndarray:
    """All posterior draws of a scalar variable as one array (chains concatenated)."""
    return idata.posterior[var].values.reshape(-1)


def combine(gl: FitResult, bl: FitResult, oc: FitResult, alpha: dict, seed=4) -> dict:
    """Notebook cells 31-35: sample-wise sum of land and ocean posteriors, independent TCR draws."""
    beta = flat(bl.idata, "beta_land") + flat(oc.idata, "beta_ocean")
    gamma = flat(gl.idata, "gamma_land") + flat(oc.idata, "gamma_ocean")
    rng = np.random.default_rng(seed)
    tcr = rng.normal(alpha["TCR_mu"], alpha["TCR_sigma"], beta.size)
    a = tcr / alpha["pco2_ref"]
    afk = 1.0 / (ud.K + beta + a * gamma)  # upstream "AF" (= AF/k, ppm PgC^-1)
    return dict(beta=beta, gamma=gamma, tcr=tcr, alpha=a, af_over_k=afk, tcre=1e3 * a * afk)


def summarise_tcre(tcre: np.ndarray) -> dict:
    """Mean, sd, 95 % HDI and 5-95 % range of TCRE draws (K/EgC)."""
    lo, hi = az.hdi(tcre, hdi_prob=0.95)
    return dict(mean=float(tcre.mean()), sd=float(tcre.std()), hdi95_lo=float(lo), hdi95_hi=float(hi),
                p05=float(np.percentile(tcre, 5)), p95=float(np.percentile(tcre, 95)))


def run(seed=0, sampler=None, iav_obs=None, progressbar=False) -> dict:
    """Full upstream pipeline. Component seeds derived from `seed`."""
    inp = inputs(iav_obs=iav_obs)
    gl = fit_gamma_land(inp, seed=seed + 1, sampler=sampler, progressbar=progressbar)
    bl = fit_beta_land(inp, seed=seed + 2, sampler=sampler, progressbar=progressbar)
    oc = fit_ocean(inp, seed=seed + 3, sampler=sampler, progressbar=progressbar)
    comb = combine(gl, bl, oc, inp["alpha"], seed=seed + 4)
    return dict(inputs=inp, fits=dict(gamma_land=gl, beta_land=bl, ocean=oc), combined=comb,
                tcre=summarise_tcre(comb["tcre"]))
