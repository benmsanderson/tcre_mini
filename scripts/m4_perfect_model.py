"""M4 / A4: perfect-model leave-one-out (loo) and leave-one-family-out (lfo) test of the upstream framework.

Targets: models with an actual 1pctCO2 TCRE (esm-deck C4MIP, COU, piControl baseline, year 61-80 mean)
and upstream ocean predictors. alpha = the held-out model's own dT/dCa, for the framework and baselines.
Baselines for the carbon factor R = TCRE/TCR (equivalently TCRE at fixed TCR):
  raw   Gaussian(mean, sd) of the training models' R (deck values), times TCR_m;
  re    random-effects predictive: Student-t, n-1 dof, scale sd*sqrt(1+1/n), times TCR_m.
Writes results/a4/<model>_<mode>.npz (posterior samples), results/m4_scores.csv, results/m4_summary.json,
results/fig_m4_pit.png
Usage: python scripts/m4_perfect_model.py [--jobs 4] [--modes loo lfo]
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

from tcre_evidence import feedbacks as fb
from tcre_evidence import perfect_model as pmod
from tcre_evidence import plotting as pl
from tcre_evidence import scoring as sc
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS

OUT = RESULTS / "a4"
SEED = 0
SAMPLER = dict(cores=1)


def run_one(m, mode, alpha):
    f = OUT / f"{m}_{mode}.npz"
    if f.exists():
        return
    r = pmod.fit_heldout(m, mode, alpha, seed=SEED, sampler=SAMPLER)
    np.savez_compressed(f, **{k: np.asarray(r[k], dtype="float32") for k in
                              ("beta_land", "gamma_land", "beta_ocean", "gamma_ocean", "af", "tcre")})
    json.dump(dict(model=m, mode=mode, heldout=r["heldout"], diagnostics=r["diagnostics"],
                   n_train_land=r["n_train_land"], n_train_ocean=r["n_train_ocean"]),
              open(OUT / f"{m}_{mode}.json", "w"), indent=1)


def truth_table():
    t = fb.table("2x", "picontrol")
    sheet = ud.cmip_sheet_dataset().drop_dims("process").to_dataframe()
    oc = ud.terhaar_ocean(69).to_dataframe()
    T = pd.DataFrame(dict(tcre=t.tcre, af=t.af, alpha=t.alpha, tcr=t.tcr, R=t.R))
    T["beta_land"], T["gamma_land"] = sheet.beta_land.reindex(T.index), sheet.gamma_land.reindex(T.index)
    T["beta_ocean"], T["gamma_ocean"] = oc.beta_ocean.reindex(T.index), oc.gamma_ocean.reindex(T.index)
    return T


def baselines(T, m, heldout):
    tr = T.drop(index=[x for x in heldout if x in T.index])
    tcr = T.loc[m, "tcr"]
    R = tr.R.values / 1e3 * 1e3  # R in 1/EgC (TCRE/TCR)
    n = len(R)
    mu, sd = R.mean(), R.std(ddof=1)
    y = T.loc[m, "tcre"]
    raw = dict(mean=tcr * mu, sd=tcr * sd)
    raw.update(crps=sc.crps_gaussian(raw["mean"], raw["sd"], y),
               logscore=float(-stats.norm.logpdf(y, raw["mean"], raw["sd"])),
               pit=float(stats.norm.cdf(y, raw["mean"], raw["sd"])))
    scale = tcr * sd * np.sqrt(1 + 1 / n)
    re = dict(mean=tcr * mu, sd=float(scale * np.sqrt((n - 1) / (n - 3))))
    draws = tcr * mu + scale * stats.t.rvs(n - 1, size=20000, random_state=SEED)
    re.update(crps=sc.crps_samples(draws, y), logscore=float(-stats.t.logpdf((y - tcr * mu) / scale, n - 1) + np.log(scale)),
              pit=float(stats.t.cdf((y - tcr * mu) / scale, n - 1)))
    for b in (raw, re):
        z = stats.norm.ppf if b is raw else (lambda q: stats.t.ppf(q, n - 1))
        s = raw["sd"] if b is raw else scale
        b["cov66"] = bool(abs(y - tcr * mu) <= z(0.83) * s)
        b["cov90"] = bool(abs(y - tcr * mu) <= z(0.95) * s)
        b["n_train"] = n
    return raw, re


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--modes", nargs="+", default=["loo", "lfo"])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    T = truth_table()
    targets = [m for m in T.index if m in set(ud.terhaar_predictors().model.values)]
    print("targets:", targets)
    tasks = [(m, mode, float(T.loc[m, "alpha"])) for mode in args.modes for m in targets]
    Parallel(n_jobs=args.jobs, verbose=5)(delayed(run_one)(*t) for t in tasks)

    rows = []
    for m, mode, _ in tasks:
        d = np.load(OUT / f"{m}_{mode}.npz")
        meta = json.load(open(OUT / f"{m}_{mode}.json"))
        dg = meta["diagnostics"]
        diag = dict(rhat_max=max(v["rhat_max"] for v in dg.values()), ess_bulk_min=min(v["ess_bulk_min"] for v in dg.values()),
                    divergences=sum(v["divergences"] for v in dg.values()))
        for q in ("tcre", "af", "beta_land", "gamma_land", "beta_ocean", "gamma_ocean"):
            y = T.loc[m, q]
            if not np.isfinite(y):
                continue
            rows.append(dict(model=m, mode=mode, quantity=q, method="framework", **sc.summarise(d[q], y), **diag,
                             heldout=";".join(meta["heldout"])))
        raw, re = baselines(T, m, meta["heldout"])
        for name, b in (("raw_ensemble", raw), ("random_effects", re)):
            rows.append(dict(model=m, mode=mode, quantity="tcre", method=name, truth=T.loc[m, "tcre"], **b,
                             heldout=";".join(meta["heldout"])))
        # component baselines: raw training ensemble of the same input values
        tr = T.drop(index=[x for x in meta["heldout"] if x in T.index])
        for q in ("beta_land", "gamma_land", "beta_ocean", "gamma_ocean"):
            src = (ud.cmip_sheet_dataset().drop_dims("process").to_dataframe()[q] if "land" in q
                   else ud.terhaar_ocean(69).to_dataframe()[q])
            v = src.drop(index=[x for x in meta["heldout"] if x in src.index]).dropna().values
            y = T.loc[m, q]
            if np.isfinite(y):
                rows.append(dict(model=m, mode=mode, quantity=q, method="raw_ensemble", truth=y, mean=v.mean(),
                                 sd=v.std(ddof=1), crps=sc.crps_gaussian(v.mean(), v.std(ddof=1), y),
                                 pit=float(stats.norm.cdf(y, v.mean(), v.std(ddof=1))),
                                 logscore=float(-stats.norm.logpdf(y, v.mean(), v.std(ddof=1))),
                                 cov66=bool(abs(y - v.mean()) <= 0.954 * v.std(ddof=1)),
                                 cov90=bool(abs(y - v.mean()) <= 1.645 * v.std(ddof=1)),
                                 heldout=";".join(meta["heldout"])))
    S = pd.DataFrame(rows)
    S.to_csv(RESULTS / "m4_scores.csv", index=False, float_format="%.4g")
    agg = (S.groupby(["mode", "quantity", "method"])
           .agg(n=("crps", "size"), crps=("crps", "mean"), logscore=("logscore", "mean"), cov66=("cov66", "mean"),
                cov90=("cov90", "mean"), mean_sd=("sd", "mean")).reset_index())
    skill = []
    for (mode, q), g in S.groupby(["mode", "quantity"]):
        f = g[g.method == "framework"].set_index("model").crps
        r = g[g.method == "raw_ensemble"].set_index("model").crps
        common = f.index.intersection(r.index)
        skill.append(dict(mode=mode, quantity=q, crpss_vs_raw=float(1 - f[common].mean() / r[common].mean()),
                          n_better=int((f[common] < r[common]).sum()), n=len(common)))
    pd.set_option("display.width", 200)
    print(agg.round(3).to_string(index=False))
    print(pd.DataFrame(skill).round(3).to_string(index=False))
    print(S[(S.quantity == "tcre")].pivot_table(index="model", columns=["mode", "method"], values="crps").round(3))
    json.dump(dict(aggregate=agg.to_dict("records"), skill=skill, targets=targets, sampler=SAMPLER),
              open(RESULTS / "m4_summary.json", "w"), indent=1, default=float)

    # figure: per-model TCRE predictive intervals vs truth (loo)
    pl.style()
    fig, axes = plt.subplots(1, len(args.modes), figsize=(3.6 * len(args.modes), 3.4), sharey=True, squeeze=False)
    for ax, mode in zip(axes[0], args.modes):
        g = S[(S.mode == mode) & (S.quantity == "tcre")]
        for j, (meth, lab) in enumerate([("framework", "upstream framework"), ("raw_ensemble", "raw ensemble")]):
            h = g[g.method == meth].set_index("model").reindex(targets)
            y = np.arange(len(targets)) + (j - 0.5) * 0.3
            lo = h["p05"] if meth == "framework" else h["mean"] - 1.645 * h["sd"]
            hi = h["p95"] if meth == "framework" else h["mean"] + 1.645 * h["sd"]
            ax.hlines(y, lo, hi, color=pl.SERIES[j], lw=2, label=f"{lab} 90 %")
            ax.plot(h["mean"], y, "o", ms=4, color=pl.SERIES[j])
        ax.plot(T.loc[targets, "tcre"], np.arange(len(targets)), "D", ms=5, color=pl.TEXT, label="actual")
        ax.set_yticks(np.arange(len(targets)), targets)
        ax.invert_yaxis()
        ax.set_title({"loo": "leave one model out", "lfo": "leave one family out"}[mode])
        ax.set_xlabel("TCRE (K EgC$^{-1}$), α = model's own")
    fig.legend(*axes[0][0].get_legend_handles_labels(), loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.08))
    pl.save(fig, RESULTS / "fig_m4_perfect_model.png")


if __name__ == "__main__":
    main()
