"""M6 / A5: prior versus data information in the upstream framework.

(1) Prior predictive of each component's reported truth (beta_land, gamma_land, beta_ocean, gamma_ocean)
    and of TCRE, against the posterior: sd ratio, mean shift, and 1D KL(posterior || prior) by KDE.
(2) Sensitivity: refit one component with one prior changed, recombine with the default fits of the
    other components, report TCRE mean and 5-95 % range.
Writes results/m6_prior_info.csv, results/m6_prior_sensitivity.csv, results/m6_a5.json,
results/fig_m6_prior_posterior.png
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json
import math
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc as pm
from joblib import Parallel, delayed
from scipy import stats

from tcre_evidence import plotting as pl
from tcre_evidence import upstream as up
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS
from tcre_evidence.upstream_patched import model_builders as mb

SEED = 0
SAMPLER = dict(cores=1)


def kl_kde(p_samples, q_samples, n=512):
    """KL(P || Q) from samples via KDE on a common grid."""
    lo = min(np.percentile(p_samples, 0.1), np.percentile(q_samples, 0.1))
    hi = max(np.percentile(p_samples, 99.9), np.percentile(q_samples, 99.9))
    x = np.linspace(lo, hi, n)
    p = stats.gaussian_kde(p_samples)(x) + 1e-12
    q = stats.gaussian_kde(q_samples)(x) + 1e-12
    p, q = p / np.trapezoid(p, x), q / np.trapezoid(q, x)
    return float(np.trapezoid(p * np.log(p / q), x))


def models(inp, pri_g=None, pri_b=None, pri_o=None):
    d, b = inp["gamma_land_data"], inp["beta_land_data"]
    ec = inp["ec"]
    return dict(
        gamma_land=lambda: mb.build_gamma_land_model(d["gamma_land_tropics"], d["gamma_land_extratropics"], d["lookup_table"],
                                                     inp["iav_df"], inp["iav_obs_mean"], inp["iav_obs_sigma"], pri_g or up.priors.gamma_land),
        beta_land=lambda: mb.build_beta_land_model(b["beta_land"], b["lookup_table"], pri_b or up.priors.beta_land),
        ocean=lambda: mb.build_joint_ocean_model(inp["beta_ocean"], inp["gamma_ocean"], ec["O_m"], ec["O_obs"],
                                                 ec["O_obs_unc"], pri_o or up.priors.joint_ocean))


def sens_variants():
    g, b, o = up.priors.gamma_land, up.priors.beta_land, up.priors.joint_ocean
    V = []

    def gv(name, f):
        p = deepcopy(g); f(p); V.append(("gamma_land", name, p))

    def bv(name, f):
        p = deepcopy(b); f(p); V.append(("beta_land", name, p))

    def ov(name, f):
        p = deepcopy(o); f(p); V.append(("ocean", name, p))

    gv("default", lambda p: None)
    gv("permafrost=0", lambda p: p["delta"].__setitem__("permafrost", {"dist": "Normal", "mu": 0, "sigma": 1e-6}))
    gv("permafrost -LogN(log 10, 1)", lambda p: p["delta"]["permafrost"].update(mu=math.log(10)))
    gv("permafrost -LogN(log 60, 1)", lambda p: p["delta"]["permafrost"].update(mu=math.log(60)))
    gv("nitrogen Normal(0, 10)", lambda p: p["delta"].__setitem__("nitrogen", {"dist": "Normal", "mu": 0, "sigma": 10}))
    gv("sigma_struct 10", lambda p: p["sigma_struct"].update(sigma=10.0))
    gv("sigma_struct 100", lambda p: p["sigma_struct"].update(sigma=100.0))
    bv("default", lambda p: None)
    bv("eta_N mu 0", lambda p: p["eta"]["nitrogen"].update(mu=0.0))
    bv("eta sigma 1", lambda p: [p["eta"][k].update(sigma=1.0) for k in p["eta"]])
    bv("no nitrogen pseudo-obs", lambda p: p.__setitem__("nitrogen_obs", {"values": [1.0], "sigma": 1e3}))
    bv("sigma_struct 0.2", lambda p: p["sigma_struct"].update(sigma=0.2))
    bv("sigma_struct 1", lambda p: p["sigma_struct"].update(sigma=1.0))
    ov("default", lambda p: None)
    ov("regression sigma 1", lambda p: p["regression"].update(sigma=1.0))
    ov("regression sigma 100", lambda p: p["regression"].update(sigma=100.0))
    ov("LKJ eta 1", lambda p: p["lkj"].update(eta=1.0))
    ov("LKJ eta 10", lambda p: p["lkj"].update(eta=10.0))
    ov("theta sigma 2", lambda p: p["theta"].update(sigma=2.0))
    ov("gamma_s_true sigma 3", lambda p: p["gamma_s_true"].update(sigma=3.0))
    return V


def fit_variant(comp, name, pri, seed):
    inp = up.inputs()
    kw = {"gamma_land": dict(pri_g=pri), "beta_land": dict(pri_b=pri), "ocean": dict(pri_o=pri)}[comp]
    with models(inp, **kw)[comp]():
        idata = pm.sample(random_seed=seed, progressbar=False, **up.SAMPLER_DEFAULTS, **SAMPLER)
    var = {"gamma_land": ["gamma_land"], "beta_land": ["beta_land"], "ocean": ["beta_ocean", "gamma_ocean"]}[comp]
    return comp, name, {v: up.flat(idata, v) for v in var}, up.diagnostics(idata, var)


def main():
    pl.style()
    inp = up.inputs()
    V = sens_variants()
    res = Parallel(n_jobs=4, verbose=5)(delayed(fit_variant)(c, n, p, SEED + i) for i, (c, n, p) in enumerate(V))
    post = {(c, n): (s, d) for c, n, s, d in res}
    # prior predictive of the reported truths
    prior = {}
    for comp, build in models(inp).items():
        with build():
            pp = pm.sample_prior_predictive(draws=4000, random_seed=SEED)
        for v in {"gamma_land": ["gamma_land"], "beta_land": ["beta_land"], "ocean": ["beta_ocean", "gamma_ocean"]}[comp]:
            prior[v] = pp.prior[v].values.reshape(-1)
    rng = np.random.default_rng(SEED)
    n = 4000
    tcr = rng.normal(up.priors.alpha["TCR_mu"], up.priors.alpha["TCR_sigma"], n)
    a = tcr / 280.0

    def tcre(bl, gl, bo, go):
        return 1e3 * a / (ud.K + bl + bo + a * (gl + go))

    D = {v: post[(c, "default")][0][v] for c, v in [("beta_land", "beta_land"), ("gamma_land", "gamma_land"),
                                                     ("ocean", "beta_ocean"), ("ocean", "gamma_ocean")]}
    with np.errstate(all="ignore"):
        prior["tcre"] = tcre(prior["beta_land"][:n], prior["gamma_land"][:n], prior["beta_ocean"][:n], prior["gamma_ocean"][:n])
    D["tcre"] = tcre(D["beta_land"], D["gamma_land"], D["beta_ocean"], D["gamma_ocean"])
    cmip = dict(beta_land=inp["beta_land_data"].beta_land.values,
                gamma_land=(inp["gamma_land_data"].gamma_land_tropics + inp["gamma_land_data"].gamma_land_extratropics).values,
                beta_ocean=inp["beta_ocean"].values, gamma_ocean=inp["gamma_ocean"].values)
    rows = []
    for v in ["beta_land", "gamma_land", "beta_ocean", "gamma_ocean", "tcre"]:
        p = prior[v][np.isfinite(prior[v])]
        p_trim = p[(p > np.percentile(p, 0.5)) & (p < np.percentile(p, 99.5))]
        q = D[v]
        rows.append(dict(quantity=v, prior_median=float(np.median(p)), prior_iqr=float(np.subtract(*np.percentile(p, [75, 25]))),
                         post_mean=float(q.mean()), post_sd=float(q.std()),
                         post_iqr=float(np.subtract(*np.percentile(q, [75, 25]))),
                         iqr_ratio_post_prior=float(np.subtract(*np.percentile(q, [75, 25])) / np.subtract(*np.percentile(p, [75, 25]))),
                         cmip_mean=float(np.mean(cmip[v])) if v in cmip else np.nan,
                         cmip_sd=float(np.std(cmip[v], ddof=1)) if v in cmip else np.nan,
                         kl_post_prior=kl_kde(q, p_trim)))
    info = pd.DataFrame(rows)
    info.to_csv(RESULTS / "m6_prior_info.csv", index=False, float_format="%.4g")
    print(info.round(3).to_string(index=False))

    srows = []
    for c, nme, _ in V:
        s, d = post[(c, nme)]
        X = dict(D)
        X.update(s)
        t = tcre(X["beta_land"], X["gamma_land"], X["beta_ocean"], X["gamma_ocean"])
        srows.append(dict(component=c, variant=nme, tcre_mean=float(t.mean()), tcre_p05=float(np.percentile(t, 5)),
                          tcre_p95=float(np.percentile(t, 95)),
                          **{f"{k}_mean": float(v.mean()) for k, v in s.items()},
                          **{f"{k}_sd": float(v.std()) for k, v in s.items()}, **d))
    S = pd.DataFrame(srows)
    S.to_csv(RESULTS / "m6_prior_sensitivity.csv", index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    print(S.drop(columns=["ess_tail_min", "n_draws"]).round(3).to_string(index=False))
    json.dump(dict(sampler={**up.SAMPLER_DEFAULTS, **SAMPLER}, pymc=pm.__version__), open(RESULTS / "m6_a5.json", "w"))

    fig, axes = plt.subplots(1, 4, figsize=(9, 2.4))
    for ax, v, lab in zip(axes, ["beta_land", "gamma_land", "beta_ocean", "gamma_ocean"],
                          ["β$_L$ (PgC ppm$^{-1}$)", "γ$_L$ (PgC K$^{-1}$)", "β$_O$ (PgC ppm$^{-1}$)", "γ$_O$ (PgC K$^{-1}$)"]):
        q = D[v]
        lo, hi = np.percentile(np.r_[q, cmip[v]], [0.5, 99.5])
        pad = 0.6 * (hi - lo)
        bins = np.linspace(lo - pad, hi + pad, 60)
        ax.hist(prior[v], bins, density=True, histtype="step", color=pl.SERIES[1], label="prior")
        ax.hist(q, bins, density=True, histtype="step", color=pl.SERIES[0], label="posterior")
        ax.plot(cmip[v], np.zeros_like(cmip[v]), "|", ms=10, color=pl.TEXT, label="CMIP6 inputs")
        ax.set_xlabel(lab)
        ax.set_yticks([])
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12))
    pl.save(fig, RESULTS / "fig_m6_prior_posterior.png")


if __name__ == "__main__":
    main()
