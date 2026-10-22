"""Figures for the process-model notebooks. Each function draws one idea and returns the figure.

The cached results are read with load(); they are written by scripts/precompute.py into
results/precomputed/<name>/<quantity>/.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import plotting as pl
from .paths import RESULTS

METHOD_LABELS = {"raw_ensemble": "raw ensemble", "all_processes_on": "all processes on (real-world estimate)",
                 "own_flags": "own flags"}


def load(name: str, quantity: str) -> dict:
    """All cached results for one quantity as a dict of tables (missing files are skipped)."""
    d = RESULTS / "precomputed" / name / quantity
    if not d.exists():
        raise FileNotFoundError(f"{d} not found: run  python scripts/precompute.py --name {name} ...  first")
    out = {"dir": d, "run": json.load(open(d / "run.json"))}
    for f in d.glob("*.csv"):
        out[f.stem] = pd.read_csv(f)
    if (d / "budget.json").exists():
        out["budget"] = json.load(open(d / "budget.json"))
    if (d / "draws.npz").exists():
        z = np.load(d / "draws.npz")
        terms = dict.fromkeys(k.split("|")[0] for k in z.files)  # file order, without duplicates
        out["draws"] = {t: (z[f"{t}|prior"], z[f"{t}|posterior"]) for t in terms}
    return out


def _process(term: str) -> str:
    return term.split("[")[1].rstrip("]") if "[" in term else term


def pretty(term: str) -> str:
    """Plain-language name for a model term: y_delta[fire] -> 'δ fire', y_rho -> 'ρ (shared fraction)'."""
    names = {"y": "real-world value", "y_sigma_struct": "σ structural", "y_rho": "ρ (shared fraction)",
             "y_sigma_common": "σ shared bias"}
    if term.startswith("y_delta["):
        return f"δ {_process(term)}"
    if term.startswith("y_eta["):
        return f"η {_process(term)}"
    return names.get(term, term)


def by_flag(values: pd.Series, flags: pd.DataFrame, groups: pd.Series | None = None, label: str = ""):
    """Each model's value, split by whether it has each process. Horizontal bars are group means."""
    pl.style()
    procs = list(flags.columns)
    fig, axs = plt.subplots(1, len(procs), figsize=(2.1 * len(procs), 2.8), sharey=True)
    fams = sorted(groups.unique()) if groups is not None else []
    colour = {f: pl.SERIES[i % len(pl.SERIES)] for i, f in enumerate(fams)}
    rng = np.random.default_rng(0)
    for ax, p in zip(np.atleast_1d(axs), procs):
        for x in (0, 1):
            m = flags.index[flags[p] == x]
            jitter = rng.uniform(-0.12, 0.12, len(m))
            c = [colour[groups[k]] for k in m] if groups is not None else pl.SERIES[x]
            ax.scatter(x + jitter, values[m], s=18, c=c, alpha=0.85)
            if len(m):
                ax.hlines(values[m].mean(), x - 0.3, x + 0.3, color=pl.TEXT, lw=1.5)
        ax.set_xticks([0, 1], [f"without\n(n={int((flags[p] == 0).sum())})", f"with\n(n={int(flags[p].sum())})"],
                      fontsize=7.5)
        ax.set_xlim(-0.6, 1.6)
        ax.set_title(p, fontsize=9)
        ax.grid(axis="x", visible=False)
    np.atleast_1d(axs)[0].set_ylabel(label)
    if groups is not None:
        handles = [plt.Line2D([], [], marker="o", ls="", color=colour[f], label=f) for f in fams]
        fig.legend(handles=handles, fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5), title="family",
                   title_fontsize=7)
    fig.tight_layout()
    return fig


def prior_posterior(res: dict, terms: list[str] | None = None):
    """Prior (grey) and posterior (colour) of each term, with the ratio of their 90 % widths."""
    pl.style()
    draws = res["draws"]
    terms = terms or list(draws)
    table = res["prior_vs_posterior"].set_index("term")
    ncol = min(4, len(terms))
    nrow = int(np.ceil(len(terms) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.3 * ncol, 2.0 * nrow), squeeze=False)
    for ax, t in zip(axs.ravel(), terms):
        pr, po = draws[t]
        lo, hi = np.percentile(np.r_[pr, po], [1, 99])
        bins = np.linspace(lo, hi, 50)
        ax.hist(pr, bins=bins, density=True, color=pl.GRID, label="prior")
        ax.hist(po, bins=bins, density=True, color=pl.SERIES[0], alpha=0.7, label="posterior")
        ax.set_title(f"{pretty(t)}\nwidth ratio {table.loc[t, 'width_ratio']:.2f}", fontsize=8)
        ax.set_yticks([])
        ax.grid(False)
    for ax in axs.ravel()[len(terms):]:
        ax.set_visible(False)
    axs[0, 0].legend(fontsize=7)
    fig.tight_layout()
    return fig


def recovery(res: dict):
    """Synthetic test: estimated effect against the true effect for each process, in the two scenarios.

    An identified effect follows the dashed 1:1 line. One that is not identified stays near its prior
    median (dotted) whatever the truth.
    """
    pl.style()
    a, b = res["synthetic_effects_summary"].set_index("term"), res["synthetic_no_effects_summary"].set_index("term")
    prior = res["prior_vs_posterior"].assign(process=lambda t: t.term.map(_process)).set_index("process")
    procs = [p for p in a.index if p != "real_world"]
    fig, axs = plt.subplots(1, len(procs), figsize=(2.2 * len(procs), 2.6))
    for ax, p in zip(np.atleast_1d(axs), procs):
        truths = [b.loc[p, "truth"], a.loc[p, "truth"]]
        est = [b.loc[p, "mean_estimate"], a.loc[p, "mean_estimate"]]
        sd = [b.loc[p, "post_sd"], a.loc[p, "post_sd"]]
        lo, hi = min(truths + est) - max(sd), max(truths + est) + max(sd)
        ax.plot([lo, hi], [lo, hi], ls="--", color=pl.TEXT2, lw=0.8)
        if p in prior.index:
            ax.axhline(prior.loc[p, "prior_median"], ls=":", color=pl.SERIES[1], lw=1)
        ax.errorbar(truths, est, yerr=sd, fmt="o", color=pl.SERIES[0], ms=5, lw=1.2)
        ax.set_title(p, fontsize=9)
        ax.set_xlabel("true effect")
    np.atleast_1d(axs)[0].set_ylabel("estimated effect (± posterior sd)")
    fig.tight_layout()
    return fig


def loo(res: dict, label: str = ""):
    """Each model held out in turn: estimate minus its own value, with 5-95 % ranges, for three predictions."""
    pl.style()
    S = res["loo_scores"]
    models = list(dict.fromkeys(S.model))
    fig, ax = plt.subplots(figsize=(max(6, 0.55 * len(models)), 3))
    for j, meth in enumerate(["raw_ensemble", "all_processes_on", "own_flags"]):
        g = S[S.method == meth].set_index("model").reindex(models)
        x = np.arange(len(models)) + (j - 1) * 0.25
        ax.vlines(x, g.p05 - g.truth, g.p95 - g.truth, color=pl.SERIES[j], lw=1.2, alpha=0.7)
        ax.plot(x, g["mean"] - g.truth, "o", ms=4, color=pl.SERIES[j], label=METHOD_LABELS[meth])
    ax.axhline(0, color=pl.TEXT2, lw=0.8)
    ax.set_xticks(np.arange(len(models)), models, rotation=40, ha="right", fontsize=7.5)
    ax.set_ylabel(f"estimate − held-out value{(' (' + label + ')') if label else ''}")
    ax.legend(fontsize=7, ncol=3, loc="upper left")
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    return fig


def power(P: pd.DataFrame):
    """Models needed to detect each effect at two standard errors, against the number available now."""
    pl.style()
    fig, ax = plt.subplots(figsize=(5, 0.5 + 0.45 * len(P)))
    y = np.arange(len(P))[::-1]
    cap = lambda v: np.minimum(v.replace(np.inf, 5000), 5000)
    ax.barh(y + 0.18, cap(P.n_needed_same_mix), 0.36, color=pl.SERIES[0], label="current flag mix")
    ax.barh(y - 0.18, cap(P.n_needed_balanced), 0.36, color=pl.SERIES[2], label="balanced flags")
    ax.axvline(P.n_models_now.iloc[0], color=pl.SERIES[1], lw=1.5, label="models available")
    ax.set_xscale("log")
    ticks = [10, 30, 100, 300, 1000, 3000]
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_yticks(y, P.process)
    ax.set_xlabel("independent models needed (two-standard-error detection)")
    ax.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False)
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    return fig
