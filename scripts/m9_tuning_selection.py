"""M9: does model development (calibration and conditional recalibration) bias flag-based process estimates?

Synthetic ensembles on the real nitrogen flags (13 models), with known truth:
  beta_L,m  = beta0 * eta^N_m * exp(u_m),  u ~ N(0, 0.25),  beta0 = 1.43, eta = 0.73 (nitrogen on in reality)
  gamma_L,m = g0 + v_m,  v ~ N(0, 25),  g0 = -35 PgC/K (no direct nitrogen effect on gamma)
  historical land sink (1960-2014, PgC) H_m = 80 beta_L,m + 0.6 gamma_L,m - LUC_m,  LUC ~ N(55, 15)
  target H_obs = the synthetic real world's own historical sink (80 beta_true + 0.6 g0 - mean LUC)
Development practices (applied in this order; each scenario switches some on):
  FACE   : beta_L,m <- beta_L,m * (beta_face / beta_L,m)^phi, phi = 0.5, beta_face = the real-world beta_L
           (models pulled part-way toward the observed CO2 response of productivity, flag-independent)
  RECAL_T: if |H_m - H_obs| > tau, rescale turnover kappa: beta, gamma <- kappa*beta, kappa*gamma so that H lands
           uniformly within +/- tau/2 of H_obs  (moves beta and gamma together)
  RECAL_B: as RECAL_T but rescales beta only (allocation-type adjustment)
  target error: recalibration towards a historical target that is itself biased (e.g. observational or
           land-use accounting error), here 25 PgC too high
  tau = 20 PgC (comparable to the GCB 1-sigma of the net land flux, 27 PgC; results/m5.json)
Estimators (as in the development configuration, least-squares stand-in for the Bayesian fit):
  eta_hat         = exp(slope) of log beta_L on the nitrogen flag
  beta_allon_hat  = exp(intercept + slope)  ("all processes on")
  R_allon         = carbon factor from beta_allon_hat, mean gamma_L, fixed ocean
  R_raw           = ensemble mean of the models' own R
  R_e2e           = flag regression on R itself, evaluated with the flag on
R = 1000 / E at 2xCO2 with dCa = 285 ppm, TCR = 1.8 K, k = 2.12, beta_O = 0.9, gamma_O = -9 (EgC^-1).
A check with the upstream Bayesian beta_L model (patched builder, upstream priors, no pseudo-observations)
is run on a few replicates per scenario.
Writes results/m9_tuning_selection.csv, results/m9_tuning_selection.json, results/fig_m9_tuning_selection.png
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed

from tcre_evidence import plotting as pl
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS

SEED = 0
NREP = 4000
NBAYES = 6
BETA0, ETA, G0 = 1.43, 0.73, -35.0
SD_U, SD_V = 0.25, 25.0
TAU, PHI, LUC_MEAN = 20.0, 0.5, 55.0
DCA, TCR, K, BO, GO = 285.0, 1.8, 2.12, 0.9, -9.0
SCEN = {
    "no development effects": dict(face=False, recal=None),
    "FACE calibration": dict(face=True, recal=None),
    "recalibration of turnover": dict(face=False, recal="T"),
    "recalibration of allocation": dict(face=False, recal="B"),
    "FACE + turnover recalibration": dict(face=True, recal="T"),
    "turnover recalibration, target 25 PgC too high": dict(face=False, recal="T", target_err=25.0),
}


def R_of(beta, gamma):
    return 1e3 / (K * DCA + (beta + BO) * DCA + (gamma + GO) * TCR)


BETA_TRUE = BETA0 * ETA
R_TRUE = R_of(BETA_TRUE, G0)
# the target is the synthetic real world's own historical sink (true parameters, mean land-use emissions),
# so recalibration pulls towards the truth and any bias comes from how it is achieved
H_OBS = 80 * BETA_TRUE + 0.6 * G0 - LUC_MEAN


def ensemble(N, rng, face, recal, target_err=0.0):
    n = len(N)
    beta = BETA0 * ETA ** N * np.exp(rng.normal(0, SD_U, n))
    gamma = G0 + rng.normal(0, SD_V, n)
    luc = rng.normal(LUC_MEAN, 15, n)
    if face:
        beta = beta * (BETA_TRUE / beta) ** PHI
    H = 80 * beta + 0.6 * gamma - luc
    nrecal = 0
    if recal:
        out = np.abs(H - (H_OBS + target_err)) > TAU
        nrecal = int(out.sum())
        tgt = H_OBS + target_err + rng.uniform(-TAU / 2, TAU / 2, n)
        if recal == "T":
            kap = (tgt + luc) / (80 * beta + 0.6 * gamma)
            beta = np.where(out, kap * beta, beta)
            gamma = np.where(out, kap * gamma, gamma)
        else:
            bnew = (tgt + luc - 0.6 * gamma) / 80
            beta = np.where(out, bnew, beta)
        beta = np.maximum(beta, 0.05)
    return beta, gamma, nrecal


def estimate(N, beta, gamma):
    X = np.c_[np.ones(len(N)), N]
    c = np.linalg.lstsq(X, np.log(beta), rcond=None)[0]
    eta_hat = np.exp(c[1])
    b_all = np.exp(c[0] + c[1])
    R_models = R_of(beta, gamma)
    ce = np.linalg.lstsq(X, R_models, rcond=None)[0]
    return dict(eta_hat=eta_hat, beta_allon=b_all, R_allon=R_of(b_all, gamma.mean()),
                R_raw=R_models.mean(), R_e2e=ce[0] + ce[1])


def bayes_eta(N, beta, models, seed):
    import pymc as pm
    from tcre_evidence import upstream as up
    from tcre_evidence.upstream_patched import model_builders as mb
    lt = ud.lookup_table().sel(model=models)
    pri = deepcopy(up.priors.beta_land)
    pri["nitrogen_obs"] = {"values": [1.0], "sigma": 1e3}  # no paired experiments: CMIP evidence only
    da = xr.DataArray(beta, coords=dict(model=models))
    with mb.build_beta_land_model(da, lt, pri):
        idata = pm.sample(draws=1000, tune=1000, chains=4, cores=1, target_accept=0.99, random_seed=seed,
                          progressbar=False)
    eta = idata.posterior["beta_land_eta"].values[..., procs_index(pri)].reshape(-1)
    d = up.diagnostics(idata, ["beta_land", "beta_land_eta"])
    return float(eta.mean()), float(idata.posterior["beta_land"].values.mean()), d["rhat_max"], d["divergences"]


def procs_index(pri):
    return list(pri["eta"].keys()).index("nitrogen")


def main():
    pl.style()
    lt = ud.lookup_table()
    models = [str(m) for m in lt.model.values]
    N = lt.sel(process="nitrogen").values.astype(float)
    rng = np.random.default_rng(SEED)
    rows = []
    for name, sc in SCEN.items():
        for r in range(NREP):
            beta, gamma, nrecal = ensemble(N, rng, **sc)
            rows.append(dict(scenario=name, rep=r, nrecal=nrecal, **estimate(N, beta, gamma)))
    D = pd.DataFrame(rows)
    D.to_csv(RESULTS / "m9_tuning_selection.csv", index=False, float_format="%.5g")
    agg = D.groupby("scenario", sort=False).agg(
        eta_hat=("eta_hat", "mean"), eta_sd=("eta_hat", "std"), beta_allon=("beta_allon", "mean"),
        R_allon=("R_allon", "mean"), R_raw=("R_raw", "mean"), R_e2e=("R_e2e", "mean"),
        R_e2e_sd=("R_e2e", "std"), nrecal=("nrecal", "mean"))
    for c in ["R_allon", "R_raw", "R_e2e"]:
        agg[f"{c}_bias_pct"] = 100 * (agg[c] / R_TRUE - 1)
    agg["eta_bias"] = agg.eta_hat - ETA
    agg["beta_allon_bias_pct"] = 100 * (agg.beta_allon / BETA_TRUE - 1)
    pd.set_option("display.width", 220)
    print(f"truth: eta {ETA}, beta_L {BETA_TRUE:.3f}, R {R_TRUE:.4f}, H {H_OBS:.1f} PgC")
    print(agg.round(3).to_string())

    # Bayesian check on a few replicates (upstream beta_L model, CMIP evidence only)
    rng_b = np.random.default_rng(SEED + 1)
    jobs = []
    for name in ["no development effects", "recalibration of turnover", "FACE + turnover recalibration"]:
        for r in range(NBAYES):
            beta, gamma, _ = ensemble(N, rng_b, **SCEN[name])
            jobs.append((name, r, beta, estimate(N, beta, gamma)["eta_hat"]))
    res = Parallel(n_jobs=int(os.environ.get("TCRE_JOBS", 4)))(
        delayed(bayes_eta)(N, b, models, SEED + 10 + i) for i, (_, _, b, _) in enumerate(jobs))
    B = pd.DataFrame([dict(scenario=n, rep=r, eta_ols=e, eta_bayes=rb[0], beta_land_bayes=rb[1], rhat=rb[2], div=rb[3])
                      for (n, r, _, e), rb in zip(jobs, res)])
    Bagg = B.groupby("scenario", sort=False).agg(eta_ols=("eta_ols", "mean"), eta_bayes=("eta_bayes", "mean"),
                                                 beta_land_bayes=("beta_land_bayes", "mean"), rhat_max=("rhat", "max"),
                                                 divergences=("div", "sum"))
    print(Bagg.round(3).to_string())
    json.dump(dict(truth=dict(eta=ETA, beta_L=BETA_TRUE, R=R_TRUE, gamma_L=G0),
                   settings=dict(nrep=NREP, nbayes=NBAYES, tau=TAU, phi=PHI, sd_u=SD_U, sd_v=SD_V, H_obs=H_OBS, luc_mean=LUC_MEAN),
                   aggregate=agg.reset_index().to_dict("records"), bayes_check=Bagg.reset_index().to_dict("records")),
              open(RESULTS / "m9_tuning_selection.json", "w"), indent=1, default=float)

    plot(agg)


def plot(agg):
    pl.style()
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.4), gridspec_kw=dict(width_ratios=[1, 1.4]))
    y = np.arange(len(agg))
    ax = axes[0]
    ax.barh(y, agg.eta_hat, color=pl.SERIES[0], height=0.5, zorder=3)
    for yi, v in zip(y, agg.eta_hat):
        ax.text(v + 0.01, yi, f"{v:.2f}", va="center", fontsize=8, color=pl.TEXT2)
    ax.axvline(ETA, color=pl.TEXT, lw=1.2, zorder=4)
    ax.text(ETA, -0.75, f" true η$_N$ = {ETA}", fontsize=8, color=pl.TEXT, va="center")
    ax.set_yticks(y, agg.index)
    ax.invert_yaxis()
    ax.set_xlim(0.5, 1.25)
    ax.set_xlabel("η$_N$ from flag regression (mean)")
    ax = axes[1]
    for j, (c, lab) in enumerate([("R_allon_bias_pct", "'all processes on'"), ("R_raw_bias_pct", "raw ensemble mean"),
                                  ("R_e2e_bias_pct", "flag regression on R")]):
        ax.barh(y + (j - 1) * 0.25, agg[c], 0.23, color=pl.SERIES[j], label=lab, zorder=3)
    ax.axvline(0, color=pl.TEXT, lw=1.2, zorder=4)
    ax.set_yticks(y, [])
    ax.invert_yaxis()
    ax.set_xlabel("bias in carbon factor R relative to truth (%)")
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=3, bbox_to_anchor=(0.6, -0.14), fontsize=8)
    pl.save(fig, RESULTS / "fig_m9_tuning_selection.png")


if __name__ == "__main__":
    import sys
    if "--plot-only" in sys.argv:
        plot(pd.DataFrame(json.load(open(RESULTS / "m9_tuning_selection.json"))["aggregate"]).set_index("scenario"))
    else:
        main()
