"""Redraw fig_m3_covariance.png from results/m3_samples.npz (written by m3_covariance_rerun.py)."""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from tcre_evidence import plotting as pl
from tcre_evidence.paths import RESULTS


def main():
    pl.style()
    s = np.load(RESULTS / "m3_samples.npz")
    x = np.linspace(0.2, 3.6, 400)
    fig, ax = plt.subplots(figsize=(6, 3))
    for i, (k, lab) in enumerate([("U", "framework as published"),
                                  ("J_a", "joint land + TCR–carbon rank correlation (deck)"),
                                  ("U_h", "TCR–carbon rank correlation (hard-coded α order)")]):
        ax.plot(x, stats.gaussian_kde(s[k])(x), color=pl.SERIES[i], label=lab)
        ax.axvline(np.percentile(s[k], 95), color=pl.SERIES[i], lw=1, ls=":")
    ax.set_xlabel("TCRE (K EgC$^{-1}$); dotted lines: 95th percentiles")
    ax.set_ylabel("density")
    ax.set_ylim(bottom=0)
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=1, bbox_to_anchor=(0.5, -0.32), fontsize=8)
    pl.save(fig, RESULTS / "fig_m3_covariance.png")


if __name__ == "__main__":
    main()
