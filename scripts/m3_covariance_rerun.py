"""M3 / A3: recombination and covariance.

(1) Deck bootstrap (11 C4MIP models, alpha from the same runs, COU-BGC): paired / upstream-pattern /
    independent recombination.
(2) Upstream pipeline variants (seed 0):
    U    upstream as published (separate land fits, independent TCR prior);
    J    joint land model (correlated biases across beta_L, gamma_L tropics/extratropics);
    U+a, J+a   TCR prior samples rank-correlated with total beta and gamma (Iman-Conover), target rank
               correlations from the deck models (alpha aligned by name);
    U+h  as U+a but target correlations from upstream's CMIP reference (hard-coded alpha order).
Writes results/m3_covariance.json, results/m3_samples.npz, results/fig_m3_covariance.png
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pymc as pm
from scipy import stats

from tcre_evidence import covariance as cv
from tcre_evidence import feedbacks as fb
from tcre_evidence import joint_land as jl
from tcre_evidence import plotting as pl
from tcre_evidence import upstream as up
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS
from tcre_evidence.rankcorr import fisher_ci, iman_conover

SEED = 0


def target_corr(tcr, beta, gamma):
    return stats.spearmanr(np.c_[tcr, beta, gamma]).statistic


def tcre_stats(tcr, beta, gamma):
    a = tcr / 280.0
    t = 1e3 * a / (ud.K + beta + a * gamma)
    hdi = az.hdi(t, hdi_prob=0.95)
    return dict(mean=float(t.mean()), sd=float(t.std()), p05=float(np.percentile(t, 5)),
                p50=float(np.percentile(t, 50)), p95=float(np.percentile(t, 95)),
                hdi95_lo=float(hdi[0]), hdi95_hi=float(hdi[1])), t


def main():
    pl.style()
    out = {}
    # (1) deck bootstrap
    t = fb.table("2x", "picontrol")
    vals = dict(alpha=t.alpha.values, beta_land=t.beta_land_coubgc.values, gamma_land=t.gamma_land_coubgc.values,
                beta_ocean=t.beta_ocean_coubgc.values, gamma_ocean=t.gamma_ocean_coubgc.values)
    out["deck_correlations"] = cv.pair_correlations(vals)
    out["deck_bootstrap"] = {s: cv.summary(cv.resample(vals, s, 400_000, SEED)) for s in cv.SCHEMES}
    print("deck correlations:", {k: round(v, 2) for k, v in out["deck_correlations"].items()})
    for s, r in out["deck_bootstrap"].items():
        print(f"  deck {s:12s} sd AF {r['sd_af']:.3f} TCRE {r['tcre_mean']:.2f} ± {r['tcre_sd']:.3f} [{r['tcre_p05']:.2f}, {r['tcre_p95']:.2f}]")

    # target rank correlations (TCR, beta_tot, gamma_tot)
    C_deck = target_corr(t.tcr.values, (t.beta_land_coubgc + t.beta_ocean_coubgc).values,
                         (t.gamma_land_coubgc + t.gamma_ocean_coubgc).values)
    cp = ud.carbon_parameters("2xCO2")
    C_hard = target_corr(ud.HARDCODED_ALPHA * 280, (cp.beta_land + cp.beta_ocean).values,
                         (cp.gamma_land + cp.gamma_ocean).values)
    n = len(t)
    out["target_rank_corr"] = dict(
        deck=C_deck.tolist(), hardcoded=C_hard.tolist(),
        deck_ci90={"tcr_beta": fisher_ci(C_deck[0, 1], n), "tcr_gamma": fisher_ci(C_deck[0, 2], n),
                   "beta_gamma": fisher_ci(C_deck[1, 2], n)})
    print("target rank corr deck (TCR,b,g):\n", np.round(C_deck, 2), "\nhard-coded:\n", np.round(C_hard, 2))

    # (2) upstream variants
    inp = up.inputs()
    fits = dict(gamma_land=up.fit_gamma_land(inp, seed=SEED + 1), beta_land=up.fit_beta_land(inp, seed=SEED + 2),
                ocean=up.fit_ocean(inp, seed=SEED + 3))
    gd, bd = inp["gamma_land_data"], inp["beta_land_data"]
    jm = jl.build_joint_land_model(bd.beta_land, gd.gamma_land_tropics, gd.gamma_land_extratropics, gd.lookup_table,
                                   inp["iav_df"], inp["iav_obs_mean"], inp["iav_obs_sigma"],
                                   up.priors.beta_land, up.priors.gamma_land)
    with jm:
        jid = pm.sample(random_seed=SEED + 5, progressbar=False, **up.SAMPLER_DEFAULTS)
    fits["joint_land"] = up.FitResult("joint_land", jid, up.diagnostics(jid, ["beta_land", "gamma_land", "bias_R"]))
    out["diagnostics"] = {k: f.diagnostics for k, f in fits.items()}
    for k, d in out["diagnostics"].items():
        print(f"  {k:10s}", {a: round(b, 4) for a, b in d.items()})
    R = jid.posterior["bias_R"].values.reshape(-1, 3, 3)
    out["joint_land_bias_corr_posterior_mean"] = R.mean(0).tolist()
    bl_s, gl_s = up.flat(fits["beta_land"].idata, "beta_land"), up.flat(fits["gamma_land"].idata, "gamma_land")
    bj_s, gj_s = up.flat(jid, "beta_land"), up.flat(jid, "gamma_land")
    out["land_truth_posterior_corr"] = dict(separate=float(np.corrcoef(bl_s, gl_s)[0, 1]),
                                            joint=float(np.corrcoef(bj_s, gj_s)[0, 1]))
    bo, go = up.flat(fits["ocean"].idata, "beta_ocean"), up.flat(fits["ocean"].idata, "gamma_ocean")
    rng = np.random.default_rng(SEED + 4)
    tcr = rng.normal(up.priors.alpha["TCR_mu"], up.priors.alpha["TCR_sigma"], bo.size)
    variants = {"U": (bl_s + bo, gl_s + go), "J": (bj_s + bo, gj_s + go)}
    res, samples = {}, {}
    for name, (b, g) in variants.items():
        res[name], samples[name] = tcre_stats(tcr, b, g)
        for tag, C in (("a", C_deck), ("h", C_hard)):
            if name == "J" and tag == "h":
                continue
            X = iman_conover(np.c_[tcr, b, g], C, seed=SEED + 6)
            res[f"{name}+{tag}"], samples[f"{name}+{tag}"] = tcre_stats(*X.T)
    out["tcre_variants"] = res
    for k, r in res.items():
        print(f"  {k:5s} TCRE {r['mean']:.3f} ± {r['sd']:.3f}  5-95 [{r['p05']:.2f}, {r['p95']:.2f}]  HDI95 [{r['hdi95_lo']:.2f}, {r['hdi95_hi']:.2f}]")
    out["sampler"] = up.SAMPLER_DEFAULTS
    out["versions"] = dict(pymc=pm.__version__, arviz=az.__version__)
    json.dump(out, open(RESULTS / "m3_covariance.json", "w"), indent=2)
    np.savez_compressed(RESULTS / "m3_samples.npz", **{k.replace("+", "_"): v.astype("float32") for k, v in samples.items()})

    # figure: TCRE distributions for the variants
    fig, ax = plt.subplots(figsize=(6, 3))
    bins = np.linspace(0, 3.6, 73)
    for i, k in enumerate(["U", "J+a", "U+h"]):
        ax.hist(samples[k], bins, histtype="step", lw=2, color=pl.SERIES[i], density=True,
                label={"U": "upstream as published", "J+a": "joint land + TCR–carbon correlation (deck)",
                       "U+h": "TCR–carbon correlation (hard-coded α order)"}[k])
    ax.set_xlabel("TCRE (K EgC$^{-1}$)")
    ax.set_ylabel("density")
    ax.legend(fontsize=7, loc="upper right")
    pl.save(fig, RESULTS / "fig_m3_covariance.png")


if __name__ == "__main__":
    main()
