"""Figures for the LaTeX technical note (note/). Reads saved results only; no sampling, so no seed is needed.

Outputs note/figures/{nitrogen,perfect_model,development_bias}.pdf from
results/m7_nitrogen_historical.csv, results/m4_scores.csv and results/m9_tuning_selection.json.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tcre_evidence import plotting as pl
from tcre_evidence.paths import RESULTS

OUT = Path(__file__).resolve().parents[1] / "note" / "figures"
Z90 = 1.645


def nitrogen():
    D = pd.read_csv(RESULTS / "m7_nitrogen_historical.csv")
    panels = [("beta_land_coubgc", r"$\beta_L$ (PgC ppm$^{-1}$, COU$-$BGC)"),
              ("R", r"$R$ = TCRE / TCR"),
              ("H_land", "net land flux 1960–2014 (PgC)")]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
    rng = np.random.default_rng(0)
    for ax, (col, lab) in zip(axs, panels):
        for i, (flag, name) in enumerate([(0, "no N limitation"), (1, "N limitation")]):
            v = D.loc[D.nitrogen == flag, col].dropna().to_numpy()
            x = i + rng.uniform(-0.08, 0.08, v.size)
            ax.plot(x, v, "o", ms=4, color=pl.SERIES[i], alpha=0.8)
            ax.hlines(v.mean(), i - 0.25, i + 0.25, color=pl.TEXT, lw=1.5)
            ns = ns + [v.size] if i else [v.size]
        ax.set_xticks([0, 1], [f"no N (n={ns[0]})", f"N-limited (n={ns[1]})"], fontsize=7.5)
        ax.set_xlim(-0.6, 1.6)
        ax.set_title(lab, fontsize=8.5)
        ax.grid(axis="x", visible=False)
    fig.tight_layout()
    pl.save(fig, OUT / "nitrogen.pdf")


def perfect_model():
    S = pd.read_csv(RESULTS / "m4_scores.csv")
    S = S[(S["mode"] == "loo") & (S.quantity == "tcre")]
    meths = [("raw_ensemble", "raw ensemble"), ("framework", "configuration (all processes on)"),
             ("framework_flags", "own-flag predictive")]
    models = sorted(S.model.unique())
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    for j, (meth, lab) in enumerate(meths):
        g = S[S.method == meth].set_index("model").reindex(models)
        err = g["mean"] - g["truth"]
        lo = np.where(g.p05.notna(), g.p05 - g.truth, err - Z90 * g["sd"])
        hi = np.where(g.p95.notna(), g.p95 - g.truth, err + Z90 * g["sd"])
        x = np.arange(len(models)) + (j - 1) * 0.25
        ax.vlines(x, lo, hi, color=pl.SERIES[j], lw=1.2, alpha=0.7)
        ax.plot(x, err, "o", ms=4, color=pl.SERIES[j], label=lab)
    ax.axhline(0, color=pl.TEXT2, lw=0.8)
    ax.set_xticks(np.arange(len(models)), models, rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel(r"estimate $-$ held-out TCRE (K EgC$^{-1}$)")
    ax.set_ylim(-1.2, 2.7)
    ax.legend(fontsize=7.5, loc="upper left", ncol=3)
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    pl.save(fig, OUT / "perfect_model.pdf")


def development_bias():
    d = json.load(open(RESULTS / "m9_tuning_selection.json"))
    A = pd.DataFrame(d["aggregate"])
    labels = ["none", "FACE", "recal. turnover", "recal. allocation", "FACE + recal.", "recal., target +25 PgC"]
    y = np.arange(len(A))[::-1]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
    a1.errorbar(A.eta_hat, y, xerr=A.eta_sd, fmt="o", color=pl.SERIES[0], ms=4, lw=1.2)
    a1.axvline(d["truth"]["eta"], color=pl.SERIES[1], lw=1, ls="--")
    a1.axvline(1.0, color=pl.TEXT2, lw=0.8, ls=":")
    a1.set_yticks(y, labels, fontsize=7.5)
    a1.set_xlabel(r"flag-estimated $\eta_N$ (truth 0.73)")
    a1.grid(axis="y", visible=False)
    est = [("R_raw_bias_pct", "raw ensemble"), ("R_allon_bias_pct", "all processes on"),
           ("R_e2e_bias_pct", "end to end on R")]
    for j, (col, lab) in enumerate(est):
        a2.barh(y + (1 - j) * 0.25, A[col], height=0.25, color=pl.SERIES[j], label=lab)
    a2.axvline(0, color=pl.TEXT2, lw=0.8)
    a2.set_xlabel(r"bias in $R$ (%)")
    a2.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, columnspacing=0.8)
    a2.grid(axis="y", visible=False)
    fig.tight_layout()
    pl.save(fig, OUT / "development_bias.pdf")


if __name__ == "__main__":
    pl.style()
    nitrogen()
    perfect_model()
    development_bias()
    print("wrote", *sorted(p.name for p in OUT.glob("*.pdf")))
