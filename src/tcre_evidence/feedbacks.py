"""Carbon-cycle feedback parameters from 1pctCO2 BGC/COU/RAD global means.

Conventions (names carry the convention):
- *_coubgc: Arora et al. (2020) 2x2 solve from COU and BGC (gamma from COU-BGC, accounting for BGC warming);
- *_bgcrad: 2x2 solve from BGC and RAD (gamma_rad = dC_RAD / dT_RAD, beta from BGC corrected for its warming);
- eps = dC_COU - (dC_BGC + dC_RAD): non-additive interaction (PgC), per reservoir.

Anomalies: temperature relative to a linear fit of piControl (or the mean of the first 10 forced years);
carbon fluxes drift-corrected by a linear fit of piControl (or not); CO2 relative to the forced run's
first year. All quantities are means over a 20-year window (default simulation years 61-80, centred on
year 70.5 for 2xCO2; years 131-140 for 4xCO2, a 10-year window because three models stop at year 140).
Units: beta PgC/ppm, gamma PgC/K, alpha K/ppm, TCRE K/EgC.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import c4mip
from .upstream_data import K

WINDOWS = {"2x": (60, 80), "4x": (130, 140)}


def _linfit(y: pd.Series, n: int) -> np.ndarray:
    t = np.arange(len(y))
    p = np.polyfit(t, y.values, 1)
    return np.polyval(p, np.arange(n))


def anomalies(runs: dict, baseline: str = "picontrol") -> dict[str, pd.DataFrame]:
    """Per forced run: dT (K), dCa (ppm), cumulative dC_land, dC_ocean (PgC, uptake positive)."""
    ctrl = runs["ctrl"]
    out = {}
    for r in ("bgc", "cou", "rad"):
        d = runs[r]
        n = len(d)
        if baseline == "picontrol":
            t0 = _linfit(ctrl.tas, n)
            fl0, fo0 = _linfit(ctrl.fland, n), _linfit(ctrl.focean, n)
        elif baseline == "forced_start":  # HANDOVER §2.4 quick computation
            t0 = runs[r].tas.iloc[:10].mean()
            fl0 = fo0 = 0.0
        else:
            raise ValueError(baseline)
        out[r] = pd.DataFrame(dict(
            dT=d.tas.values - t0,
            dCa=d.co2.values - runs["cou"].co2.iloc[0],
            dCL=np.cumsum(d.fland.values - fl0),
            dCO=np.cumsum(d.focean.values - fo0),
        ))
    return out


def window_means(an: dict, level: str = "2x") -> dict[str, pd.Series]:
    a, b = WINDOWS[level]
    res = {}
    for r, df in an.items():
        if len(df) < b:
            raise ValueError(f"run {r} has {len(df)} years, window needs {b}")
        res[r] = df.iloc[a:b].mean()
    return res


def diagnose(model: str, level: str = "2x", baseline: str = "picontrol") -> dict:
    """All per-model quantities for one deck model at one forcing level."""
    runs = c4mip.read_model(model)
    c4mip.assert_signs(runs, model)
    w = window_means(anomalies(runs, baseline), level)
    cou, bgc, rad = w["cou"], w["bgc"], w["rad"]
    dCa = cou.dCa
    out = dict(model=c4mip.MODELS[model].cmip_name, deck_name=model, family=c4mip.MODELS[model].family,
               level=level, baseline=baseline, dCa=dCa, dT_cou=cou.dT, dT_bgc=bgc.dT, dT_rad=rad.dT)
    for x, col in (("land", "dCL"), ("ocean", "dCO")):
        Cc, Cb, Cr = cou[col], bgc[col], rad[col]
        out[f"dC_{x}_cou"], out[f"dC_{x}_bgc"], out[f"dC_{x}_rad"] = Cc, Cb, Cr
        g_cb = (Cc - Cb) / (cou.dT - bgc.dT)
        out[f"gamma_{x}_coubgc"] = g_cb
        out[f"beta_{x}_coubgc"] = (Cb - g_cb * bgc.dT) / dCa
        g_r = Cr / rad.dT
        out[f"gamma_{x}_rad"] = g_r
        out[f"beta_{x}_bgcrad"] = (Cb - g_r * bgc.dT) / dCa
        out[f"beta_{x}_bgconly"] = Cb / dCa
        out[f"eps_{x}"] = Cc - (Cb + Cr)
        # COU-RAD allocation: gamma from RAD, beta absorbs eps (beta_courad = (C_cou - gamma_rad dT_cou)/dCa)
        out[f"beta_{x}_courad"] = (Cc - g_r * cou.dT) / dCa
    E = K * dCa + out["dC_land_cou"] + out["dC_ocean_cou"]
    out.update(alpha=cou.dT / dCa, tcr=cou.dT, E=E, tcre=1e3 * cou.dT / E, af=K * dCa / E, R=1e3 / E)
    return out


def table(level: str = "2x", baseline: str = "picontrol", models=None) -> pd.DataFrame:
    models = models or list(c4mip.MODELS)
    return pd.DataFrame([diagnose(m, level, baseline) for m in models]).set_index("model")
