"""M5 / Workstream B: TCRE = TCR x R with R constrained by historical observables, and perfect-model
validation of the historical-to-1pctCO2 (and -to-flat10) transfer.

R_m = TCRE_m / TCR_m (EgC^-1) from esm-deck 1pctCO2 (COU, piControl baseline, years 61-80);
R_flat10_m = TCRE_flat10_m / TCR_m. H in {H_ocean, H_land, H_af} (historical 1960-2014, see
tcre_evidence.historical). Assessed TCR: N(1.7, 0.36) K (upstream prior) and N(1.8, 0.365) K (AR6 very
likely 1.2-2.4).
Writes results/m5_table.csv, results/m5_estimates.csv, results/m5_perfect_model.csv, results/m5.json,
results/fig_m5_relationships.png, results/fig_m5_estimates.png
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tcre_evidence import feedbacks as fb
from tcre_evidence import historical as hist
from tcre_evidence import perfect_model as pmod
from tcre_evidence import plotting as pl
from tcre_evidence import scoring as sc
from tcre_evidence import tcr_r
from tcre_evidence.paths import RESULTS
from tcre_evidence.rankcorr import fisher_ci

SEED = 0
N = 20000
TCR_PRIORS = {"upstream N(1.7, 0.36)": (1.7, 0.36), "AR6 N(1.8, 0.365)": (1.8, 0.365)}
HS = ["H_ocean", "H_land", "H_af"]


def build_table():
    t = fb.table("2x", "picontrol")
    H = hist.table()
    f = hist.flat10_tcre()
    T = pd.DataFrame(dict(tcr=t.tcr, tcre=t.tcre, R=t.R, af=t.af))
    T = T.join(H[HS]).join(f.tcre_flat10)
    T["R_flat10"] = T.tcre_flat10 / T.tcr
    T["pathway_diff"] = T.tcre - T.tcre_flat10
    T["family"] = [pmod.family(m) for m in T.index]
    return T


def estimates(T, o):
    rows, samples = [], {}
    rng = np.random.default_rng(SEED)
    for pname, (mu, sd) in TCR_PRIORS.items():
        tcr_star = rng.normal(mu, sd, N)
        for target in ["R", "R_flat10"]:
            for mode, hname in [("ensemble", None), ("unconstrained", None)] + [("substitute", h) for h in HS]:
                d = T.dropna(subset=[target] + ([hname] if hname else []))
                Rs = tcr_r.predict_R(d.tcr.values, d[target].values, d[hname].values if hname else None, tcr_star,
                                     o[hname] if hname else None, o[f"sigma_{hname}"] if hname else 0.0, mode, seed=SEED)
                te = tcr_star * Rs
                key = f"{pname}|{target}|{mode}|{hname or '-'}"
                samples[key] = te
                rows.append(dict(tcr_prior=pname, target=target, mode=mode, H=hname or "-", n=len(d),
                                 R_mean=Rs.mean(), R_sd=Rs.std(), tcre_mean=te.mean(), tcre_sd=te.std(),
                                 tcre_p05=np.percentile(te, 5), tcre_p17=np.percentile(te, 17),
                                 tcre_p50=np.percentile(te, 50), tcre_p83=np.percentile(te, 83),
                                 tcre_p95=np.percentile(te, 95)))
    return pd.DataFrame(rows), samples


def perfect_model(T, o):
    rows = []
    for target in ["R", "R_flat10"]:
        for hname in HS:
            d = T.dropna(subset=[target, hname])
            for lfo, obs_scale in ((False, 1.0), (True, 1.0), (False, 0.0)):
                for m in d.index:
                    out = [x for x in d.index if (pmod.family(x) == pmod.family(m) if lfo else x == m)]
                    tr = d.drop(index=out)
                    if len(tr) < 5:
                        continue
                    y = d.loc[m, target]
                    h_m = d.loc[m, hname]  # exact pseudo-observation; obs error enters via sigma_obs
                    for mode in ["ensemble", "unconstrained", "H_only", "substitute"]:
                        s = tcr_r.predict_R_given_tcr(tr.tcr.values, tr[target].values, tr[hname].values,
                                                      d.loc[m, "tcr"], h_m, obs_scale * o[f"sigma_{hname}"], mode, n=4000, seed=SEED)
                        r = sc.summarise(s * d.loc[m, "tcr"], y * d.loc[m, "tcr"])  # score on TCRE scale
                        rows.append(dict(target=target, H=hname, mode=mode, lfo=lfo, obs_err=obs_scale, model=m, n_train=len(tr), **r))
    return pd.DataFrame(rows)


def main():
    pl.style()
    T = build_table()
    o = hist.observed()
    T.to_csv(RESULTS / "m5_table.csv", float_format="%.4g")
    pd.set_option("display.width", 220)
    print(T.round(3).to_string())
    print("obs:", {k: round(v, 3) for k, v in o.items()})
    corr = {}
    for target in ["R", "R_flat10", "tcr"]:
        for h in HS + (["tcr"] if target != "tcr" else []):
            d = T.dropna(subset=[target, h])
            if target == h or len(d) < 4:
                continue
            r = float(np.corrcoef(d[target], d[h])[0, 1])
            corr[f"{target}~{h}"] = dict(r=r, n=len(d), ci90=fisher_ci(r, len(d)))
    for k, v in corr.items():
        print(f"  corr {k:16s} r={v['r']:+.2f} n={v['n']} 90% CI [{v['ci90'][0]:+.2f}, {v['ci90'][1]:+.2f}]")
    est, samples = estimates(T, o)
    est.to_csv(RESULTS / "m5_estimates.csv", index=False, float_format="%.4g")
    print(est[est.target == "R"].round(3).to_string(index=False))
    pmt = perfect_model(T, o)
    pmt.to_csv(RESULTS / "m5_perfect_model.csv", index=False, float_format="%.4g")
    agg = (pmt.groupby(["target", "H", "lfo", "obs_err", "mode"]).agg(n=("crps", "size"), crps=("crps", "mean"),
           logscore=("logscore", "mean"), cov66=("cov66", "mean"), cov90=("cov90", "mean"), sd=("sd", "mean"))
           .reset_index())
    base = agg[agg["mode"] == "ensemble"].set_index(["target", "H", "lfo", "obs_err"]).crps
    agg["crpss_vs_ensemble"] = 1 - agg.crps / agg.set_index(["target", "H", "lfo", "obs_err"]).index.map(base).values
    print(agg[(agg.lfo == False)].round(3).to_string(index=False))
    js = dict(observed=o, correlations=corr, perfect_model=agg.to_dict("records"),
              pathway_diff=dict(mean=float(T.pathway_diff.mean()), sd=float(T.pathway_diff.std()),
                                n=int(T.pathway_diff.notna().sum())),
              tcr_priors=TCR_PRIORS, n_samples=N, seed=SEED)
    json.dump(js, open(RESULTS / "m5.json", "w"), indent=1, default=float)
    np.savez_compressed(RESULTS / "m5_samples.npz", **{k.replace("|", "__").replace(" ", "_").replace("(", "").replace(")", "").replace(",", ""): v.astype("float32") for k, v in samples.items()})

    # relationships figure
    fig, axes = plt.subplots(1, 3, figsize=(9, 2.8), sharey=True)
    for ax, h, lab in zip(axes, HS, ["cumulative ocean sink 1960–2014 (PgC)", "cumulative net land flux 1960–2014 (PgC)",
                                     "airborne fraction of fossil emissions 1960–2014"]):
        d = T.dropna(subset=["R", h])
        ax.plot(d[h], d.R, "o", ms=5, color=pl.SERIES[0])
        for m, r in d.iterrows():
            ax.annotate(m.split("-")[0], (r[h], r.R), fontsize=6, color=pl.TEXT2, xytext=(3, 2), textcoords="offset points")
        ax.axvspan(o[h] - o[f"sigma_{h}"], o[h] + o[f"sigma_{h}"], color=pl.SERIES[1], alpha=0.15, lw=0)
        ax.axvline(o[h], color=pl.SERIES[1], lw=1.5)
        ax.set_xlabel(lab, fontsize=8)
        c = corr[f"R~{h}"]
        ax.set_title(f"r = {c['r']:+.2f} (n = {c['n']})", fontsize=8)
    axes[0].set_ylabel("R = TCRE / TCR (EgC$^{-1}$)")
    pl.save(fig, RESULTS / "fig_m5_relationships.png")


if __name__ == "__main__":
    main()
