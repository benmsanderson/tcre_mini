"""M11: process completeness of the flag-based land model.

(1) Residual budget: share of between-model variance in log beta_L and gamma_L (COU-BGC, development-code
    inputs) explained by the four presence flags; ordinary and leave-one-out R2, with and without
    BCC-CSM2-MR, and repeated on the esm-deck gamma_L (11 models). Spread within each flag group.
(2) Real-data fit of the single-variable gamma_L process model: posterior of the shared-bias terms
    (sigma_struct, rho, sigma_common), the only part of the model that can absorb a process no model has.
(3) Synthetic tests on the real 13-model flag table with an unflagged mechanism (labelled 'demography'):
      U0  present in no model; the synthetic real world has an extra effect d_unrep on gamma_L;
      U1  present, unflagged, in one model (EC-Earth3-CC, the cohort-based vegetation model);
      U3  present, unflagged, in three of the four vegetation-flagged models.
    Coverage of the synthetic real-world value by the 90 % interval of the process model's
    real-world estimate (U * prod eta + sum delta), and where the unflagged effect is attributed.
(4) An explicit unrepresented-process term d_unrep on gamma_L added to the development configuration
    posterior (results/m6_prior_samples.npz) and its effect on TCRE.
Writes results/m11_completeness.json, results/m11_synthetic.csv, results/fig_m11_completeness.png
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc as pm
import xarray as xr
from joblib import Parallel, delayed

from tcre_evidence import plotting as pl
from tcre_evidence import upstream as up
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS
from tcre_evidence.upstream_patched import CMIP_process, model_builders as mb

SEED = 0
NREP = 20
PROC = ["nitrogen", "permafrost", "fire", "veg"]
TRUE_DELTA = {"nitrogen": 10.0, "permafrost": -30.0, "fire": -5.0, "veg": 5.0}
BASE = -15.0
D_DEM = -20.0  # illustrative size of the unflagged effect on gamma_L (PgC/K)
SUBSETS = {"U0": [], "U1": ["EC-Earth3-CC"], "U3": ["EC-Earth3-CC", "GFDL-ESM4", "MPI-ESM1-2-LR"]}
D_UNREP = [0.0, -10.0, -20.0, -30.0, -40.0, -60.0]
SAMPLER = dict(draws=1000, tune=1000, chains=4, cores=1, target_accept=0.99, progressbar=False)


# ---------------------------------------------------------------- (1) residual budget
def ols_budget(y, F):
    X = np.column_stack([np.ones(len(y)), F])
    n, p = X.shape
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    r2 = 1 - (r ** 2).sum() / ((y - y.mean()) ** 2).sum()
    loo = []
    for i in range(n):
        m = np.arange(n) != i
        bb, *_ = np.linalg.lstsq(X[m], y[m], rcond=None)
        loo.append(y[i] - X[i] @ bb)
    q2 = 1 - np.sum(np.square(loo)) / np.sum(np.square(y - y.mean()))
    return dict(n=n, sd=float(y.std(ddof=1)), r2=float(r2), adj_r2=float(1 - (1 - r2) * (n - 1) / (n - p)),
                loo_q2=float(q2), resid_sd=float(np.sqrt((r ** 2).sum() / (n - p))))


def residual_budget():
    ds = ud.cmip_sheet_dataset()
    D = ds[["beta_land", "gamma_land", "gamma_land_tropics", "gamma_land_extratropics"]].to_dataframe().dropna()
    D["log_beta_land"] = np.log(D.beta_land)
    F = ud.lookup_table().to_pandas().T.reindex(D.index)[PROC]
    rows = []
    for drop in ([], ["BCC-CSM2-MR"]):
        k = ~D.index.isin(drop)
        for q in ["log_beta_land", "gamma_land", "gamma_land_tropics", "gamma_land_extratropics"]:
            rows.append(dict(source="development inputs", quantity=q, dropped=",".join(drop) or "none",
                             **ols_budget(D.loc[k, q].values, F.loc[k].values)))
    E = pd.read_csv(RESULTS / "m7_nitrogen_historical.csv").set_index("model").dropna(subset=["gamma_land_coubgc"])
    Fe = ud.lookup_table().to_pandas().T.reindex(E.index)[PROC]
    for q, col in [("log_beta_land", "beta_land_coubgc"), ("gamma_land", "gamma_land_coubgc")]:
        y = np.log(E[col].values) if q.startswith("log") else E[col].values
        rows.append(dict(source="esm-deck 1pctCO2", quantity=q, dropped="none", **ols_budget(y, Fe.values)))
    groups = {p: dict(models=list(F.index[F[p] == 1]), gamma_land=D.gamma_land[F[p] == 1].round(1).tolist())
              for p in PROC}
    return pd.DataFrame(rows), groups


# ---------------------------------------------------------------- (2) real-data shared-bias terms
def build(y, lt):
    pri = up.priors.gamma_land
    return CMIP_process.build_vectorized_process_model(
        y, lt, additive_processes=PROC, noprocess_prior=mb._noprocess_prior(pri["baseline"], "g"),
        sigma_struct=pri["sigma_struct"], sigma_process_struct=pri["sigma_process_struct"],
        delta_priors=mb._delta_priors(pri["delta"]), likelihood_noise=pri["likelihood_noise"], var_name="g")


def real_fit():
    ds = ud.cmip_sheet_dataset()
    y = ds.gamma_land.dropna("model")
    lt = ud.lookup_table().reindex(model=y.model)
    with build(y, lt) as m:
        idata = pm.sample(random_seed=SEED + 11, **SAMPLER)
        prior = pm.sample_prior_predictive(draws=4000, random_seed=SEED)
    out = dict(diagnostics=up.diagnostics(idata, ["g", "g_delta", "g_sigma_struct", "g_rho"]))
    for v in ["g", "g_sigma_struct", "g_rho", "g_sigma_common", "g_sigma_m"]:
        s = up.flat(idata, v)
        p = prior.prior[v].values.reshape(-1)
        out[v] = dict(post_mean=float(s.mean()), post_sd=float(s.std()), post_p05=float(np.percentile(s, 5)),
                      post_p95=float(np.percentile(s, 95)), prior_mean=float(np.nanmean(p)), prior_sd=float(np.nanstd(p)))
    return out


# ---------------------------------------------------------------- (3) synthetic unflagged mechanism
def generate(subset, rng, lt):
    M = lt.sizes["model"]
    L = lt.sel(process=PROC).values.astype(float)
    s_proc = np.array([up.priors.gamma_land["sigma_process_struct"][p] for p in PROC]) * 0.25
    eps = rng.normal(0, s_proc[:, None], L.shape)
    delta = np.array([TRUE_DELTA[p] for p in PROC])
    dem = np.where(np.isin(lt.model.values, subset), D_DEM, 0.0)
    y = (BASE + (L * (delta[:, None] + eps)).sum(0) + rng.normal(0, 8.0) + rng.normal(0, 8.0, M)
         + rng.normal(0, 5.0, M) + dem)
    return xr.DataArray(y, coords=dict(model=lt.model.values))


def fit_synth(case, rep):
    rng = np.random.default_rng(5000 + 100 * list(SUBSETS).index(case) + rep)
    lt = ud.lookup_table()
    y = generate(SUBSETS[case], rng, lt)
    with build(y, lt):
        idata = pm.sample(random_seed=SEED + rep, **SAMPLER)
    g = up.flat(idata, "g")
    d = idata.posterior["g_delta"].values.reshape(-1, len(PROC))
    lo, hi = np.percentile(g, [5, 95])
    truth0 = BASE + sum(TRUE_DELTA.values())
    row = dict(case=case, rep=rep, g_mean=float(g.mean()), g_sd=float(g.std()), g_p05=float(lo), g_p95=float(hi),
               truth_without=truth0, **up.diagnostics(idata, ["g", "g_delta"]))
    for i, p in enumerate(PROC):
        row[f"delta_{p}"] = float(d[:, i].mean())
    for du in D_UNREP:  # U0: real world carries du, no model does
        row[f"cov90_unrep{int(du)}"] = bool(lo <= truth0 + du <= hi)
    row["cov90_with_dem"] = bool(lo <= truth0 + D_DEM <= hi)  # U1/U3: real world has the mechanism
    row["cov90_without_dem"] = bool(lo <= truth0 <= hi)
    return row


# ---------------------------------------------------------------- (4) unrepresented-process term
def unrep_term():
    S = np.load(RESULTS / "m6_prior_samples.npz")
    bl, gl, bo, go = (S[f"post_{v}"].astype(float) for v in ["beta_land", "gamma_land", "beta_ocean", "gamma_ocean"])
    rng = np.random.default_rng(SEED + 4)
    n = bl.size
    a = rng.normal(up.priors.alpha["TCR_mu"], up.priors.alpha["TCR_sigma"], n) / 280.0

    def summ(du):
        t = 1e3 * a / (ud.K + bl + bo + a * (gl + du + go))
        return dict(mean=float(t.mean()), p05=float(np.percentile(t, 5)), p50=float(np.median(t)),
                    p95=float(np.percentile(t, 95)))

    out = {"none": summ(0.0)}
    for s in (10.0, 20.0, 30.0):
        out[f"normal_sd{int(s)}"] = summ(rng.normal(0, s, n))
        out[f"neg_halfnormal_sd{int(s)}"] = summ(-np.abs(rng.normal(0, s, n)))
    return out


def main():
    pl.style()
    budget, groups = residual_budget()
    pd.set_option("display.width", 200)
    print(budget.round(2).to_string(index=False))
    tasks = [(c, r) for c in SUBSETS for r in range(NREP)]
    jobs = int(os.environ.get("TCRE_JOBS", 4))
    res = Parallel(n_jobs=jobs, verbose=5)(delayed(fit_synth)(c, r) for c, r in tasks)
    R = pd.DataFrame(res)
    R.to_csv(RESULTS / "m11_synthetic.csv", index=False, float_format="%.4g")
    agg = R.groupby("case").agg(**{f"delta_{p}": (f"delta_{p}", "mean") for p in PROC},
                                g_mean=("g_mean", "mean"), g_sd=("g_sd", "mean"),
                                cov90_with_dem=("cov90_with_dem", "mean"),
                                cov90_without_dem=("cov90_without_dem", "mean"),
                                rhat_max=("rhat_max", "max"), divergences=("divergences", "sum"))
    print(agg.round(2).to_string())
    u0 = R[R.case == "U0"]
    cov_unrep = {int(du): float(u0[f"cov90_unrep{int(du)}"].mean()) for du in D_UNREP}
    print("U0 coverage by unrepresented effect:", cov_unrep)
    real = real_fit()
    print({k: v for k, v in real.items() if k != "diagnostics"}, real["diagnostics"])
    term = unrep_term()
    print(json.dumps(term, indent=1))
    js = dict(budget=budget.to_dict("records"), flag_groups=groups, synthetic=agg.reset_index().to_dict("records"),
              u0_cov90_by_unrep=cov_unrep, truth_without=BASE + sum(TRUE_DELTA.values()), d_dem=D_DEM,
              subsets=SUBSETS, real_fit=real, unrep_term=term, nrep=NREP, sampler=SAMPLER, seed=SEED)
    json.dump(js, open(RESULTS / "m11_completeness.json", "w"), indent=1, default=float)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8, 2.8))
    a1.plot(list(cov_unrep), list(cov_unrep.values()), "o-", color=pl.SERIES[0])
    a1.axhline(0.9, color=pl.TEXT2, lw=0.8, ls="--")
    a1.set_xlabel("effect in the real world, absent from all models (PgC K$^{-1}$)")
    a1.set_ylabel("90 % interval coverage")
    a1.invert_xaxis()
    for j, c in enumerate(["U0", "U1", "U3"]):
        a2.bar(np.arange(len(PROC)) + (j - 1) * 0.27, [agg.loc[c, f"delta_{p}"] for p in PROC], 0.27,
               color=pl.SERIES[j], label={"U0": "none", "U1": "1 model", "U3": "3 veg-flagged"}[c])
    a2.scatter(range(len(PROC)), [TRUE_DELTA[p] for p in PROC], marker="_", s=300, color=pl.TEXT, zorder=3)
    a2.set_xticks(range(len(PROC)), PROC)
    a2.set_ylabel(r"posterior mean $\delta$ (PgC K$^{-1}$)")
    a2.legend(title="unflagged mechanism in", fontsize=7, title_fontsize=7)
    pl.save(fig, RESULTS / "fig_m11_completeness.png")


if __name__ == "__main__":
    main()
