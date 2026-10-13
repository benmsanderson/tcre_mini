"""M3 / A2: convention sensitivity (RAD reconstruction bias, eps) and variance attribution of the carbon
term under three allocations, at 2x and 4x, piControl and first-10-year baselines.

Writes results/m3_attribution.csv, results/m3_attribution.json, results/fig_m3_attribution.png,
results/fig_m3_eps.png, results/m3_eps_per_model.csv
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tcre_evidence import attribution as at
from tcre_evidence import feedbacks as fb
from tcre_evidence import identity as idt
from tcre_evidence import plotting as pl
from tcre_evidence.paths import RESULTS


def main():
    pl.style()
    rows, js = [], {}
    tabs = {}
    for level in ["2x", "4x"]:
        for base in ["picontrol", "forced_start"]:
            t = fb.table(level, base)
            tabs[(level, base)] = t
            key = f"{level}_{base}"
            js[key] = dict(tcre_log_shares=at.tcre_log_shares(t))
            for alloc in at.ALLOCATIONS:
                sh = at.shares(at.components(t, alloc))
                for comp, r in sh.iterrows():
                    rows.append(dict(level=level, baseline=base, allocation=alloc, component=comp, **r.to_dict()))
                beta_share = sh.loc[["B_land", "B_ocean"], "cov_share"].sum()
                js[key][alloc] = dict(beta_cov_share=float(beta_share), gamma_cov_share=float(1 - beta_share),
                                      land_cov_share=float(sh.loc[["B_land", "G_land"], "cov_share"].sum()),
                                      beta_indep_share=float(sh.loc[["B_land", "B_ocean"], "indep_share"].sum()))
            r = idt.reconstruct(t, "bgcrad")
            js[key]["rad_reconstruction"] = dict(
                af_err_mean=float(r.af_err.mean()), af_err_sd=float(r.af_err.std()),
                tcre_err_mean=float(r.tcre_err.mean()), tcre_err_sd=float(r.tcre_err.std()),
                n_within_0p05=int((r.tcre_err.abs() < 0.05).sum()), worst=r.tcre_err.abs().idxmax(),
                worst_err=float(r.tcre_err[r.tcre_err.abs().idxmax()]))
            js[key]["eps"] = {f"eps_{x}_{s}": float(getattr(t[f"eps_{x}"], s)()) for x in ("land", "ocean")
                              for s in ("mean", "std")}
            js[key]["eps"]["n_eps_land_pos"] = int((t.eps_land > 0).sum())
            js[key]["eps"]["n_eps_ocean_neg"] = int((t.eps_ocean < 0).sum())
            print(f"\n== {key}  ln-TCRE variance share: alpha {js[key]['tcre_log_shares']['alpha']:.2f}, "
                  f"carbon {js[key]['tcre_log_shares']['carbon']:.2f}")
            for alloc in at.ALLOCATIONS:
                d = js[key][alloc]
                print(f"   {alloc:10s} beta share {d['beta_cov_share']:+.2f} (indep {d['beta_indep_share']:.2f}), "
                      f"land share {d['land_cov_share']:+.2f}")
            print("   RAD recon:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in js[key]['rad_reconstruction'].items()})
            print("   eps:", {k: round(v, 1) for k, v in js[key]['eps'].items()})
    pd.DataFrame(rows).to_csv(RESULTS / "m3_attribution.csv", index=False, float_format="%.4g")
    # eps and gamma per model, both levels (state dependence, A7 preview)
    e = []
    for (level, base), t in tabs.items():
        e.append(t[["family", "eps_land", "eps_ocean", "gamma_land_coubgc", "gamma_land_rad", "gamma_ocean_coubgc",
                    "gamma_ocean_rad", "beta_land_coubgc", "beta_ocean_coubgc", "dT_cou", "dCa"]]
                 .assign(level=level, baseline=base).reset_index())
    pd.concat(e).to_csv(RESULTS / "m3_eps_per_model.csv", index=False, float_format="%.4g")
    json.dump(js, open(RESULTS / "m3_attribution.json", "w"), indent=2)

    # Figure 1: share of carbon-term variance from beta (land+ocean) under each allocation
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
    for ax, level in zip(axes, ["2x", "4x"]):
        x = np.arange(len(at.ALLOCATIONS))
        for j, (base, lab) in enumerate([("picontrol", "piControl baseline"), ("forced_start", "first-10-year baseline")]):
            v = [js[f"{level}_{base}"][a]["beta_cov_share"] for a in at.ALLOCATIONS]
            ax.bar(x + (j - 0.5) * 0.36, v, 0.34, color=pl.SERIES[j], label=lab)
        ax.axhline(0.5, color=pl.TEXT2, lw=0.8, ls=":")
        ax.set_xticks(x, ["COU−BGC\n(C4MIP)", "COU−RAD", "symmetric"])
        ax.set_title(f"{level}CO$_2$")
    axes[0].set_ylabel("share of var(β + αγ) from β terms")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.08))
    pl.save(fig, RESULTS / "fig_m3_attribution.png")

    # Figure 2: eps per model, land and ocean, 2x and 4x (piControl baseline)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for ax, x in zip(axes, ["land", "ocean"]):
        t2, t4 = tabs[("2x", "picontrol")], tabs[("4x", "picontrol")]
        y = np.arange(len(t2))
        ax.barh(y + 0.18, t2[f"eps_{x}"], 0.34, color=pl.SERIES[0], label="2xCO$_2$")
        ax.barh(y - 0.18, t4[f"eps_{x}"], 0.34, color=pl.SERIES[1], label="4xCO$_2$")
        ax.axvline(0, color=pl.TEXT2, lw=0.8)
        ax.set_yticks(y, t2.index)
        ax.set_title(f"ε$_{{{x}}}$ = ΔC$_{{COU}}$ − ΔC$_{{BGC}}$ − ΔC$_{{RAD}}$ (PgC)")
        ax.invert_yaxis()
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.06))
    pl.save(fig, RESULTS / "fig_m3_eps.png")


if __name__ == "__main__":
    main()
