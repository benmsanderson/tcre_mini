"""Historical observables H for Workstream B (concentration-driven `historical`, 1960-2014) and their
observational counterparts from GCB 2025.

Model side (PgC, uptake positive):
  H_ocean = cumulative fgco2, minus a pre-industrial baseline: the piControl mean (first 200 years) when
            it agrees with the 1850-1869 historical mean to within 0.25 PgC/yr, else the 1850-1869 mean;
  H_land  = cumulative nbp (net land flux incl. land-use change), minus the esm-deck piControl mean land
            flux where available (Pangeo has no piControl nbp for most models), else no baseline;
  H_af    = G / (G + H_land + H_ocean), with G the observed atmospheric growth (GCB) over the window:
            prescribed CO2 follows observations, so models differ only through their sinks.
Observations (GCB 2025), defined to match the model quantities: model nbp includes land-use change,
so G + H_land + H_ocean is the model's compatible fossil emission and H_af is the airborne fraction of
fossil emissions. Hence
  H_ocean_obs = ocean_sink; H_land_obs = fossil - G - ocean_sink (net land flux as a budget residual,
  avoiding the bookkeeping LUC estimate and the budget imbalance); H_af_obs = G / fossil.
1-sigma observational errors, fully correlated in time over the window: fossil 5 %, ocean 0.4 PgC/yr
(GCB stated 1-sigma); G is treated as exact.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import c4mip, obs
from .paths import DERIVED

CMIP = DERIVED / "cmip_globalmean"
WINDOW = (1960, 2014)


def _series(exp, model):
    f = CMIP / exp / f"{model}.csv"
    return pd.read_csv(f).set_index("year") if f.exists() else None


def _deck_ctrl_land(model):
    deck = c4mip.CMIP_TO_DECK.get(model)
    if deck is None:
        return None
    return float(c4mip.read_run(deck, "ctrl").fland.mean())


def model_observables(model: str, G: float, window=WINDOW) -> dict:
    """Cumulative historical land and ocean sinks and airborne fraction for one model over `window`."""
    h = _series("historical", model)
    out = dict(model=model)
    if h is None:
        return out
    a, b = window
    n = b - a + 1
    sel = (h.index >= a) & (h.index <= b)
    early = h.index < 1870
    if "fgco2" in h and h.fgco2[sel].notna().all():
        c = _series("piControl", model)
        e = float(h.fgco2[early].mean())
        base, src = e, "historical 1850-1869"
        if c is not None and "fgco2" in c and abs(float(c.fgco2.mean()) - e) < 0.25:
            base, src = float(c.fgco2.mean()), "piControl mean"
        out.update(H_ocean=float(h.fgco2[sel].sum() - n * base), ocean_baseline=base, ocean_baseline_src=src)
    if "nbp" in h and h.nbp[sel].notna().all():
        cl = _deck_ctrl_land(model)
        out.update(H_land=float(h.nbp[sel].sum() - n * (cl or 0.0)), land_baseline=cl or 0.0,
                   land_baseline_src="esm-deck piControl mean" if cl is not None else "none")
    if "H_ocean" in out and "H_land" in out:
        out["H_af"] = G / (G + out["H_land"] + out["H_ocean"])
    return out


def observed(window=WINDOW) -> dict:
    """The same observables from the Global Carbon Budget over `window`."""
    g = obs.gcb()
    a, b = window
    s = g.loc[a:b]
    n = len(s)
    G, FF, O = float(s.atm_growth.sum()), float(s.fossil.sum()), float(s.ocean_sink.sum())
    sFF, sO = 0.05 * FF, 0.4 * n
    return dict(G=G, fossil=FF, n_years=n,
                H_ocean=O, sigma_H_ocean=sO,
                H_land=FF - G - O, sigma_H_land=float(np.hypot(sFF, sO)),
                H_af=G / FF, sigma_H_af=G / FF * 0.05,
                gcb_land_sink_minus_luc=float((s.land_sink - s.luc).sum()),
                gcb_budget_imbalance=float(s.budget_imbalance.sum()))


def table(window=WINDOW) -> pd.DataFrame:
    """model_observables() for every model with a downloaded historical run."""
    o = observed(window)
    models = sorted(p.stem for p in (CMIP / "historical").glob("*.csv"))
    return pd.DataFrame([model_observables(m, o["G"], window) for m in models]).set_index("model")


def flat10_tcre() -> pd.DataFrame:
    """flat10 TCRE at 1000 PgC from esm-deck output (temperature anomaly, mean of years 91-110)."""
    from .paths import ESM_DECK
    T = pd.read_csv(ESM_DECK / "output" / "temperatures_flat10.csv")
    name = {"UKESM": "UKESM1-0-LL", "GISS": "GISS-E2-1-G"}
    t = T.iloc[90:110].mean().rename(index=lambda c: name.get(c, c))
    met = pd.read_csv(ESM_DECK / "esm_metrics_flat10.csv", index_col=0).rename(index=lambda c: name.get(c, c))
    return pd.DataFrame(dict(tcre_flat10=t, tcre_flat10_esmdeck=met.tcre.reindex(t.index)))
