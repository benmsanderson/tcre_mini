"""Two figures on the treatment of paired factorials in the WIEMIP process model, from saved results.

1. paired_intervals.png: for each process, the 90 % interval and median of the prior, of the posterior with
   the development code's treatment of the factorials (one value per model, fixed error), and of the
   posterior with the between-model spread estimated (results/wiemip/paired_spread_wiemip.csv).
2. paired_coverage.png: in pseudo-ensembles whose factorials scatter between models as the real ones do,
   the share of replicates whose 90 % interval holds the true effect, for the same two treatments
   (results/wiemip/paired_tests_synthetic_summary.csv, scenario 'observed').

Writes results/wiemip/figures/. No sampling; seconds.
"""
import matplotlib.pyplot as plt
import pandas as pd

from tcre_evidence import plotting as pl
from tcre_evidence.paths import RESULTS

W = RESULTS / "wiemip"
FIG = W / "figures"
PROCS = ["nitrogen", "fire", "veg", "permafrost"]
LABELS = {"nitrogen": "Nitrogen", "fire": "Fire", "veg": "Dynamic vegetation", "permafrost": "Permafrost"}
QUANTS = {"gamma_land_coubgc": (r"$\gamma_L$", r"Change in $\gamma_L$ from the process (PgC/K)", 0.0, (-70, 60)),
          "beta_land_coubgc": (r"$\beta_L$", r"Factor on $\beta_L$ from the process ($\eta$)", 1.0, (0.45, 1.9))}
PRIOR = "#a9a8a2"                                    # neutral: the prior is reference, not a series
DEFAULT, ESTIMATE = pl.SERIES[0], pl.SERIES[1]
ROWS = [("prior", "Prior", PRIOR, "|"),
        ("upstream_models", "Factorials, fixed error (current default)", DEFAULT, "o"),
        ("estimate", "Factorials, between-model spread estimated", ESTIMATE, "s")]


def intervals():
    T = pd.read_csv(W / "paired_spread_wiemip.csv")
    T["process"] = T.term.str.extract(r"\[(.*)\]")[0]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
    for ax, (q, (_, xlabel, null, xlim)) in zip(axes, QUANTS.items()):
        for i, p in enumerate(PROCS):
            for j, (mode, _, colour, marker) in enumerate(ROWS):
                src = T[(T.quantity == q) & (T.process == p) & (T["mode"] == ("upstream_models" if mode == "prior"
                                                                                 else mode))].iloc[0]
                lo, mid, hi = (src[f"{'prior' if mode == 'prior' else 'post'}_{k}"] for k in ("p05", "median", "p95"))
                y = i + (j - 1) * 0.24
                ax.plot([max(lo, xlim[0]), min(hi, xlim[1])], [y, y], color=colour, lw=2.4, solid_capstyle="round")
                if lo < xlim[0]:                     # the permafrost prior runs off the axis
                    ax.plot(xlim[0], y, marker="<", color=colour, ms=5, clip_on=False)
                if mode != "prior":
                    ax.plot(mid, y, marker=marker, color=colour, ms=6, mec=pl.SURFACE, mew=1.2, zorder=3)
        ax.axvline(null, color=pl.TEXT2, lw=0.8, ls=(0, (3, 3)), zorder=0)
        ax.set_xlim(xlim)
        ax.set_xlabel(xlabel)
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(range(len(PROCS)), [LABELS[p] for p in PROCS])
    axes[0].invert_yaxis()
    handles = [plt.Line2D([], [], color=c, lw=2.4, marker=None if m == "|" else m, ms=6, mec=pl.SURFACE)
               for _, _, c, m in ROWS]
    fig.legend(handles, [lab for _, lab, _, _ in ROWS], loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("Process effects in WIEMIP: 90 % ranges (dashed line = no effect)", x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    pl.save(fig, FIG / "paired_intervals.png")


def coverage():
    S = pd.read_csv(W / "paired_tests_synthetic_summary.csv")
    S = S[(S.scenario == "observed") & S.term.isin(PROCS)]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 2.9), sharey=True)
    for ax, q in zip(axes, QUANTS):
        for i, p in enumerate(PROCS):
            for j, (mode, _, colour, marker) in enumerate(ROWS[1:]):
                v = S[(S.quantity == q) & (S.term == p) & (S.treatment == ("upstream" if mode == "upstream_models"
                                                                         else mode))].cov90.iloc[0]
                ax.plot(v, i + (j - 0.5) * 0.22, marker=marker, color=colour, ms=7, mec=pl.SURFACE, mew=1.2, ls="")
        ax.axvline(0.9, color=pl.TEXT2, lw=0.8, ls=(0, (3, 3)), zorder=0)
        ax.text(0.9, -0.75, "target", color=pl.TEXT2, ha="center", va="bottom", fontsize=8)
        ax.set_xlim(-0.03, 1.03)
        ax.set_xlabel("Share of 90 % ranges that contain the true effect")
        ax.set_title(QUANTS[q][0], loc="left")
        if q == "gamma_land_coubgc":
            i = PROCS.index("permafrost")
            ax.annotate("true value lies outside the prior", (0.02, i), xytext=(0.12, i), va="center",
                        color=pl.TEXT2, fontsize=8, arrowprops=dict(arrowstyle="-", color=pl.TEXT2, lw=0.6))
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(range(len(PROCS)), [LABELS[p] for p in PROCS])
    axes[0].set_ylim(len(PROCS) - 0.5, -0.9)
    handles = [plt.Line2D([], [], color=c, marker=m, ms=7, mec=pl.SURFACE, ls="") for _, _, c, m in ROWS[1:]]
    fig.legend(handles, [lab for _, lab, _, _ in ROWS[1:]], loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("Synthetic test: factorials that disagree between models as much as the real ones",
                 x=0.06, ha="left")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    pl.save(fig, FIG / "paired_coverage.png")


if __name__ == "__main__":
    pl.style()
    intervals()
    coverage()
