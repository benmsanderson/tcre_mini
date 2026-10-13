"""M6 / A6: process identifiability in the upstream additive/multiplicative process design.

Synthetic truth experiments on the real lookup table (13 models x 4 flags), refitted with the upstream
process model (patched CMIP_process.build_vectorized_process_model, upstream gamma_land priors):
  S1 'clean'      data generated from the model's own structure with known delta_p;
  S2 'family'     as S1 plus a land-model-family offset (CLM models CESM2, NorESM2-LM, CMCC-ESM2 share
                  +/-15 PgC/K) with delta_permafrost = 0 in truth;
  S3 'eps'        nitrogen scales beta by eta_N = 0.7 and gamma contains kappa * beta * dCa / dT
                  (the COU-BGC eps channel), with no direct nitrogen effect on gamma (delta_N = 0);
Plus (4) real data: the extratropical gamma model with and without the unused paired-experiment nu
evidence (nu_nitrogen from ACCESS-ESM CN-C and UKESM1-1 ctrl-nonlim, nu_fire, nu_vegetation from
UKESM1-1) as pseudo-observations on the summed tropical + extratropical delta_p.
Writes results/m6_identifiability.csv, results/m6_nu_evidence.json, results/fig_m6_identifiability.png
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc as pm
import xarray as xr
from joblib import Parallel, delayed

from tcre_evidence import plotting as pl
from tcre_evidence import upstream as up
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS, SNAPSHOT_DATA
from tcre_evidence.upstream_patched import CMIP_process, model_builders as mb

SEED = 0
NREP = 12
PROC = ["nitrogen", "permafrost", "fire", "veg"]
CLM = ["CESM2", "NorESM2-LM", "CMCC-ESM2"]
TRUE_DELTA = {"nitrogen": 10.0, "permafrost": -30.0, "fire": -5.0, "veg": 5.0}
SAMPLER = dict(draws=1000, tune=1000, chains=4, cores=1, target_accept=0.99, progressbar=False)


def generate(scenario, rng, lt):
    M = lt.sizes["model"]
    L = lt.sel(process=PROC).values.astype(float)  # (P, M)
    pri = up.priors.gamma_land
    base = -15.0
    s_proc = np.array([pri["sigma_process_struct"][p] for p in PROC]) * 0.25  # moderate model-level scatter
    d = dict(TRUE_DELTA)
    extra = np.zeros(M)
    if scenario == "family":
        d["permafrost"] = 0.0
        fam = rng.choice([-15.0, 15.0])
        extra = np.where(np.isin(lt.model.values, CLM), fam, 0.0)
    if scenario == "eps":
        d["nitrogen"] = 0.0
        beta0, eta, kappa, dCa_over_dT = 1.2, 0.7, -0.02, 280 / 1.8  # kappa: d(beta)/dT-type degradation
        beta_m = beta0 * eta ** L[0] * np.exp(rng.normal(0, 0.15, M))
        extra = kappa * beta_m * dCa_over_dT - kappa * beta0 * dCa_over_dT
    eps = rng.normal(0, s_proc[:, None], L.shape)
    delta = np.array([d[p] for p in PROC])
    y = base + (L * (delta[:, None] + eps)).sum(0) + rng.normal(0, 8.0) + rng.normal(0, 8.0, M) + extra + rng.normal(0, 5.0, M)
    return xr.DataArray(y, coords=dict(model=lt.model.values)), d


def fit_synth(scenario, rep):
    rng = np.random.default_rng(1000 * {"clean": 1, "family": 2, "eps": 3}[scenario] + rep)
    lt = ud.lookup_table()
    y, d = generate(scenario, rng, lt)
    pri = up.priors.gamma_land
    m = CMIP_process.build_vectorized_process_model(
        y, lt, additive_processes=PROC, noprocess_prior=mb._noprocess_prior(pri["baseline"], "g"),
        sigma_struct=pri["sigma_struct"], sigma_process_struct=pri["sigma_process_struct"],
        delta_priors=mb._delta_priors(pri["delta"]), likelihood_noise=pri["likelihood_noise"], var_name="g")
    with m:
        idata = pm.sample(random_seed=SEED + rep, **SAMPLER)
    out = []
    dg = up.diagnostics(idata, ["g_delta"])
    for i, p in enumerate(PROC):
        s = idata.posterior["g_delta"].values[..., i].reshape(-1)
        lo, hi = np.percentile(s, [5, 95])
        out.append(dict(scenario=scenario, rep=rep, process=p, truth=d[p], post_mean=float(s.mean()),
                        post_sd=float(s.std()), cov90=bool(lo <= d[p] <= hi), **dg))
    return out


def nu_evidence():
    """Paired-experiment gamma differences (as upstream data_utils.load_process_evidence, 2xCO2)."""
    df = pd.read_csv(SNAPSHOT_DATA / "TCREsource_betagamma.csv").set_index("model")
    g = df.gamma_L_2xCO2
    return {"nitrogen": [float(g["ACCESS-ESM_CN"] - g["ACCESS-ESM_C"]), float(g["UKESM1-1_ctrl"] - g["UKESM1-1_nonlim"])],
            "fire": [float(g["UKESM1-1_fire"] - g["UKESM1-1_ctrl"])],
            "veg": [float(g["UKESM1-1_ctrl"] - g["UKESM1-1_nodgvm"])]}


def fit_real(with_nu, sigma_nu=5.0):
    inp = up.inputs()
    d = inp["gamma_land_data"]
    m = mb.build_gamma_land_model(d["gamma_land_tropics"], d["gamma_land_extratropics"], d["lookup_table"],
                                  inp["iav_df"], inp["iav_obs_mean"], inp["iav_obs_sigma"], up.priors.gamma_land)
    nu = nu_evidence()
    if with_nu:
        procs = list(d["lookup_table"].process.values)
        with m:
            for p, vals in nu.items():
                i = procs.index(p)
                tot = m["gamma_land_tropics_delta"][i] + m["gamma_land_extratropics_delta"][i]
                pm.Normal(f"nu_{p}", mu=tot, sigma=sigma_nu, observed=np.array(vals))
    with m:
        idata = pm.sample(random_seed=SEED + 7, **SAMPLER)
    procs = list(d["lookup_table"].process.values)
    tot = idata.posterior["gamma_land_tropics_delta"] + idata.posterior["gamma_land_extratropics_delta"]
    res = dict(with_nu=with_nu, nu=nu, diagnostics=up.diagnostics(idata, ["gamma_land"]),
               gamma_land_mean=float(idata.posterior.gamma_land.mean()), gamma_land_sd=float(idata.posterior.gamma_land.std()))
    for i, p in enumerate(procs):
        s = tot.values[..., i].reshape(-1)
        res[f"delta_{p}_mean"], res[f"delta_{p}_sd"] = float(s.mean()), float(s.std())
    return res, up.flat(idata, "gamma_land")


def main():
    pl.style()
    tasks = [(s, r) for s in ("clean", "family", "eps") for r in range(NREP)]
    res = Parallel(n_jobs=4, verbose=5)(delayed(fit_synth)(s, r) for s, r in tasks)
    R = pd.DataFrame([x for r in res for x in r])
    R.to_csv(RESULTS / "m6_identifiability.csv", index=False, float_format="%.4g")
    pri = up.priors.gamma_land["delta"]
    prior_sd = {"nitrogen": float(np.std(np.random.default_rng(0).lognormal(pri["nitrogen"]["mu"], 1, 100000))),
                "permafrost": float(np.std(np.random.default_rng(0).lognormal(pri["permafrost"]["mu"], 1, 100000))),
                "fire": 10.0, "veg": 10.0}
    agg = R.groupby(["scenario", "process"]).agg(truth=("truth", "first"), mean_est=("post_mean", "mean"),
                                                 bias=("post_mean", lambda s: 0), rmse=("post_mean", "std"),
                                                 post_sd=("post_sd", "mean"), cov90=("cov90", "mean"),
                                                 rhat_max=("rhat_max", "max"), divergences=("divergences", "sum")).reset_index()
    agg["bias"] = agg.mean_est - agg.truth
    agg["post_sd_over_prior_sd"] = agg.apply(lambda r: r.post_sd / prior_sd[r.process], axis=1)
    pd.set_option("display.width", 200)
    print(agg.round(2).to_string(index=False))
    real = Parallel(n_jobs=2)(delayed(fit_real)(w) for w in (False, True))
    js = dict(aggregate=agg.to_dict("records"), prior_sd=prior_sd, nrep=NREP, sampler=SAMPLER,
              real=[r[0] for r in real])
    for r, _ in real:
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items() if k != "diagnostics"})
    json.dump(js, open(RESULTS / "m6_nu_evidence.json", "w"), indent=1, default=float)

    fig, axes = plt.subplots(1, 3, figsize=(9, 2.8), sharey=True)
    for ax, s, title in zip(axes, ["clean", "family", "eps"],
                            ["S1 generated by the model's own structure", "S2 CLM-family offset, no permafrost effect",
                             "S3 ε channel, no direct nitrogen effect on γ"]):
        g = R[R.scenario == s]
        for i, p in enumerate(PROC):
            h = g[g.process == p]
            ax.plot(h.post_mean, np.full(len(h), i) + np.linspace(-0.25, 0.25, len(h)), "o", ms=3, color=pl.SERIES[0])
            ax.plot([h.truth.iloc[0]] * 2, [i - 0.4, i + 0.4], color=pl.TEXT, lw=2)
        ax.set_yticks(range(4), [f"δ$_{{{p}}}$" for p in PROC])
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("posterior mean (PgC K$^{-1}$); bar = truth")
        ax.invert_yaxis()
    pl.save(fig, RESULTS / "fig_m6_identifiability.png")


if __name__ == "__main__":
    main()
