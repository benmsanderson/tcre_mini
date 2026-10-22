"""Variance attribution of AF and TCRE spread under three allocations of the interaction term (A2).

Exact per-model decomposition of the carbon denominator, for x in {land, ocean}:
    S = (dC_land + dC_ocean)_COU / dCa = sum_x (B_x + G_x),   G_x = gamma_x * dT_COU / dCa = alpha * gamma_x,
with B_x = S_x - G_x so every allocation closes exactly. Allocations differ only in gamma_x:
    'coubgc'    gamma_coubgc (C4MIP; eps sits in gamma)
    'courad'    gamma_rad (eps sits in beta)
    'symmetric' (gamma_coubgc + gamma_rad) / 2 (eps shared equally; Stein and Alpert 1993 / Shapley split
                of the two-factor interaction)
Variance shares use the covariance (Shapley for a sum) decomposition: share_i = cov(X_i, S) / var(S),
which sums to one and assigns covariance terms half to each partner. The JF2020-style independent
shares var(X_i) / sum var(X_j) are reported alongside.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .upstream_data import K

ALLOCATIONS = ("coubgc", "courad", "symmetric")


def components(tab: pd.DataFrame, allocation: str) -> pd.DataFrame:
    """Split each model's carbon term into CO2 (B_x) and warming (G_x = alpha*gamma) parts, x = land, ocean,
    with gamma taken from COU-BGC, RAD (courad) or their mean (symmetric)."""
    out = pd.DataFrame(index=tab.index)
    for x in ("land", "ocean"):
        S_x = tab[f"dC_{x}_cou"] / tab.dCa
        if allocation == "coubgc":
            g = tab[f"gamma_{x}_coubgc"]
        elif allocation == "courad":
            g = tab[f"gamma_{x}_rad"]
        elif allocation == "symmetric":
            g = 0.5 * (tab[f"gamma_{x}_coubgc"] + tab[f"gamma_{x}_rad"])
        else:
            raise ValueError(allocation)
        out[f"G_{x}"] = tab.alpha * g
        out[f"B_{x}"] = S_x - out[f"G_{x}"]
    return out


def shares(comp: pd.DataFrame) -> pd.DataFrame:
    """Share of the variance of the summed carbon term from each component: covariance share and independent share."""
    S = comp.sum(axis=1)
    cov = pd.Series({c: np.cov(comp[c], S)[0, 1] for c in comp}) / S.var()
    ind = comp.var() / comp.var().sum()
    return pd.DataFrame({"cov_share": cov, "indep_share": ind, "mean": comp.mean(), "sd": comp.std()})


def tcre_log_shares(tab: pd.DataFrame) -> dict:
    """ln TCRE = ln alpha - ln(k + S): covariance shares of alpha versus the carbon term."""
    la = np.log(tab.alpha)
    lc = -np.log(K + tab.dC_land_cou / tab.dCa + tab.dC_ocean_cou / tab.dCa)
    lt = la + lc
    v = lt.var()
    return dict(alpha=float(np.cov(la, lt)[0, 1] / v), carbon=float(np.cov(lc, lt)[0, 1] / v),
                corr_alpha_carbon=float(np.corrcoef(la, lc)[0, 1]),
                sd_ln_alpha=float(la.std()), sd_ln_carbon=float(lc.std()))
