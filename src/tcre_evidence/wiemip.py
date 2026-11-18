"""WIEMIP 1pctCO2 land feedbacks from the global land carbon extracted by scripts/wiemip_extract.py.

Inputs: data/wiemip/cland_annual.csv (model, factorial, sim, forcing, year, cland_pgc) and
data/wiemip/forcing.csv (year, co2_ppm, gmst_k). The extraction is the only step that reads WIEMIP files
(through `wiemip_registry`); this module does not import the registry.

Differences from the C4MIP CMIP6 setting (feedbacks.diagnose):
- The runs are land only. The climate is prescribed from FaIR (median config) through pattern scaling,
  so every model sees the same GMST; dT is the FaIR GMST anomaly, identical in COU and RAD and for the
  three patterns (ukesm, gfdl, ipsl). The pattern changes only the spatial distribution.
- BGC (and ctrl) use a constant climate, so dT_bgc = 0 and beta_coubgc = beta_bgcrad = dC_bgc / dCa.
- No ocean and no interactive temperature: TCRE, R and the airborne fraction cannot be computed.
- Each factorial (noFire, noNitrogen, ...) has its own ctrl; anomalies are taken from it.
- Global land carbon is a stock (PgC), so no flux sign convention applies. Signs are checked by
  assert_signs (BGC land uptake at 2xCO2 must be positive).

Conventions follow feedbacks.py: windows from feedbacks.WINDOWS (2x: years 61-80, rows 60:80 from 1850),
beta in PgC/ppm, gamma in PgC/K, eps = dC_cou - dC_bgc - dC_rad in PgC.
"""
from __future__ import annotations

import re
import warnings

import numpy as np
import pandas as pd

from .feedbacks import WINDOWS
from .paths import DATA

WIEMIP_DATA = DATA / "wiemip"
PATTERNS = ("ukesm", "gfdl", "ipsl")
START_YEAR = 1850
RUN_KEYS = ["model", "factorial", "sim", "forcing"]


def load(path=WIEMIP_DATA) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(cland, forcing) as written by scripts/wiemip_extract.py; forcing indexed by year."""
    cl = pd.read_csv(path / "cland_annual.csv")
    fc = pd.read_csv(path / "forcing.csv").set_index("year")
    return cl, fc


def processes(factorial: str) -> frozenset[str]:
    """Processes switched off in a factorial: 'noFire_noNitrogen' -> {'fire', 'nitrogen'}; 'baseline' -> {}."""
    if factorial == "baseline":
        return frozenset()
    toks = factorial.split("_")
    if not all(re.fullmatch(r"no[A-Z]\w*", t) for t in toks):
        raise ValueError(f"unrecognised factorial {factorial!r}")
    return frozenset(t[2:].lower() for t in toks)


def screen(cl: pd.DataFrame, level: str = "2x") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop runs that cannot be used at this level; returns (kept, dropped with a reason).

    A run is dropped when it does not cover the window or has a non-finite value anywhere in it
    (an unmasked fill value, for example).
    """
    a, b = WINDOWS[level]
    years = range(START_YEAR + a, START_YEAR + b)
    bad = []
    for key, g in cl.groupby(RUN_KEYS):
        w = g.set_index("year").cland_pgc.reindex(years)
        if w.isna().any():
            bad.append((*key, "window not covered"))
        elif not np.isfinite(w).all():
            bad.append((*key, "non-finite values in window"))
    dropped = pd.DataFrame(bad, columns=RUN_KEYS + ["reason"])
    keep = cl.merge(dropped[RUN_KEYS], how="left", indicator=True).query("_merge == 'left_only'").drop(columns="_merge")
    return keep, dropped


def ctrl_drift(cl: pd.DataFrame) -> pd.Series:
    """Linear trend of each ctrl run's land carbon over its whole length (PgC per century)."""
    c = cl[cl.sim == "ctrl"].replace([np.inf, -np.inf], np.nan).dropna()
    return c.groupby(["model", "factorial"]).apply(lambda g: 100 * np.polyfit(g.year, g.cland_pgc, 1)[0])


def feedbacks(cl: pd.DataFrame, fc: pd.DataFrame, level: str = "2x") -> pd.DataFrame:
    """beta and gamma for every (model, factorial, pattern) with ctrl, bgc and cou; one row each.

    Runs with unusable windows are dropped first (screen); a factorial whose ctrl or bgc is dropped
    gives no rows, and a missing RAD leaves its columns NaN.
    """
    a, b = WINDOWS[level]
    years = list(range(START_YEAR + a, START_YEAR + b))
    cl, dropped = screen(cl, level)
    for _, r in dropped.iterrows():
        warnings.warn(f"WIEMIP {r.model} {r.factorial} {r.sim}/{r.forcing} dropped: {r.reason}")
    dCa = fc.co2_ppm.loc[years].mean() - fc.co2_ppm.loc[START_YEAR]
    dT = fc.gmst_k.loc[years].mean() - fc.gmst_k.loc[START_YEAR]
    w = cl[cl.year.isin(years)].groupby(RUN_KEYS).cland_pgc.mean()
    rows = []
    for (m, fac, sim, frc), c0 in w.items():
        if sim != "ctrl" or (m, fac, "bgc", "stable") not in w:
            continue
        dC_bgc = w[(m, fac, "bgc", "stable")] - c0
        for p in PATTERNS:
            if (m, fac, "cou", p) not in w:
                continue
            dC_cou = w[(m, fac, "cou", p)] - c0
            dC_rad = w[(m, fac, "rad", p)] - c0 if (m, fac, "rad", p) in w else np.nan
            g_cb = (dC_cou - dC_bgc) / dT  # dT_bgc = 0
            rows.append(dict(model=m, factorial=fac, pattern=p, level=level, dCa=dCa, dT=dT,
                             dC_land_bgc=dC_bgc, dC_land_cou=dC_cou, dC_land_rad=dC_rad,
                             beta_land_coubgc=dC_bgc / dCa, gamma_land_coubgc=g_cb,
                             gamma_land_rad=dC_rad / dT, eps_land=dC_cou - dC_bgc - dC_rad))
    out = pd.DataFrame(rows)
    assert_signs(out)
    return out


def assert_signs(fb: pd.DataFrame) -> None:
    """Stock sign check: BGC land uptake at the window must be positive for every row."""
    bad = fb[~(fb.dC_land_bgc > 0)]
    if len(bad):
        raise AssertionError(f"non-positive BGC land uptake:\n{bad[['model', 'factorial', 'dC_land_bgc']]}")


QUANTITIES = ["beta_land_coubgc", "gamma_land_coubgc", "gamma_land_rad", "eps_land"]


def values_table(fb: pd.DataFrame, pattern: str = "ukesm", factorial: str = "baseline") -> pd.DataFrame:
    """One row per model in the toolkit's values format (model, quantities...).

    pattern : one of PATTERNS, or 'mean' for the mean over the patterns a model ran (column n_patterns
              records how many; gamma_land_rad is averaged over the patterns with RAD).
    """
    f = fb[fb.factorial == factorial]
    if pattern == "mean":
        t = f.groupby("model")[QUANTITIES].mean()
        t["n_patterns"] = f.groupby("model").size()
    else:
        t = f[f.pattern == pattern].set_index("model")[QUANTITIES]
    return t.reset_index()


def factorial_flags(cl: pd.DataFrame) -> pd.DataFrame:
    """Process flags implied by the factorials each model ran: 1 where a model switched a process off.

    Absence of a factorial says nothing (a model may have the process and not have run it off), so the
    other cells are NaN, not 0. This is a cross-check for a full flag table, not a substitute.
    """
    procs = cl.groupby("model").factorial.unique().map(lambda fs: set().union(*map(processes, fs)))
    cols = sorted(set().union(*procs))
    return pd.DataFrame({p: procs.map(lambda s: 1.0 if p in s else np.nan) for p in cols}).rename_axis("model")


def paired_effects(fb: pd.DataFrame) -> pd.DataFrame:
    """Process effects from pairs of factorials of the same model that differ by one process.

    process on vs off, with the other switched-off processes as `background` ('' for baseline vs noX):
    eta_beta_land = beta_on / beta_off (the development code's eta, e.g. nitrogen 0.69, 0.77),
    delta_gamma_land_coubgc and delta_gamma_land_rad = gamma_on - gamma_off (PgC/K).
    """
    rows = []
    idx = fb.set_index(["model", "pattern", "factorial"])
    for (m, p), g in fb.groupby(["model", "pattern"]):
        facs = {f: processes(f) for f in g.factorial}
        for f_on, off_on in facs.items():
            for f_off, off_off in facs.items():
                extra = off_off - off_on
                if not (off_on < off_off and len(extra) == 1):
                    continue
                on, off = idx.loc[(m, p, f_on)], idx.loc[(m, p, f_off)]
                rows.append(dict(model=m, pattern=p, process=next(iter(extra)),
                                 background="_".join(sorted(off_on)), factorial_on=f_on, factorial_off=f_off,
                                 eta_beta_land=on.beta_land_coubgc / off.beta_land_coubgc,
                                 delta_gamma_land_coubgc=on.gamma_land_coubgc - off.gamma_land_coubgc,
                                 delta_gamma_land_rad=on.gamma_land_rad - off.gamma_land_rad))
    return pd.DataFrame(rows).sort_values(["process", "model", "background", "pattern"]).reset_index(drop=True)


def paired_by_model(pe: pd.DataFrame, column: str, background: str = "",
                    as_series: bool = False) -> dict[str, list[float]] | dict[str, pd.Series]:
    """One paired value per model and process, for process_tests.build(paired=...).

    The three climate patterns of one model are not independent experiments (beta does not depend on
    the pattern at all), so they are averaged: the arithmetic mean for delta, the geometric mean for eta.
    Only pairs with the given `background` are used ('' = baseline against noX). Non-finite values are
    dropped. Processes are named as flag columns (FLAG_NAMES). Returns {process: [value per model]} in
    model order, or {process: pd.Series indexed by model} if as_series (for process_tests.leave_one_out).
    """
    d = pe[(pe.background.fillna("") == background) & np.isfinite(pe[column])]
    if column.startswith("eta_"):
        agg = d.groupby(["process", "model"])[column].agg(lambda x: float(np.exp(np.log(x).mean())))
    else:
        agg = d.groupby(["process", "model"])[column].mean()
    procs = agg.index.get_level_values("process").unique()
    return {FLAG_NAMES.get(p, p): (agg.loc[p] if as_series else agg.loc[p].tolist()) for p in procs}


FLAG_SHEET = WIEMIP_DATA / "process_inclusion_sheet.csv"
SHEET_NAMES = {"LPX-Bem": "LPX_Bern"}  # spelling in the WIEMIP spreadsheet -> registry name
FLAG_NAMES = {"dynveg": "veg"}          # process switched off by a factorial (noDynVeg) -> flag column


def read_flag_sheet(path=FLAG_SHEET, factorial_flags: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Process flags from the WIEMIP 'process inclusion' sheet, with registry model names.

    Cells marked '?' are unknown (NaN). With `factorial_flags` (from factorial_flags()), an unknown cell is
    set to 1 where the model ran the noX factorial, since a model can only switch off a process it has.
    Returns (flags, conflicts); conflicts lists cells the sheet sets to 0 although the model ran noX.
    """
    F = pd.read_csv(path, dtype=str).set_index("model")
    F.index = [SHEET_NAMES.get(m, m.replace("-", "_")) for m in F.index]
    F = F.apply(pd.to_numeric, errors="coerce").rename_axis("model")
    conflicts = []
    if factorial_flags is not None:
        factorial_flags = factorial_flags.rename(columns=FLAG_NAMES)
        for p in F.columns.intersection(factorial_flags.columns):
            ran = factorial_flags[p].reindex(F.index) == 1
            conflicts += [dict(model=m, process=p, sheet=0) for m in F.index[ran & (F[p] == 0)]]
            F.loc[ran & F[p].isna(), p] = 1.0
    return F, pd.DataFrame(conflicts, columns=["model", "process", "sheet"])
