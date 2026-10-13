"""JF2020 identity: AF and TCRE reconstruction from (alpha, beta, gamma), and closure checks."""
from __future__ import annotations

import pandas as pd

from .upstream_data import K


def af(alpha, beta, gamma, k=K):
    """Dimensionless airborne fraction k / (k + beta + alpha*gamma)."""
    return k / (k + beta + alpha * gamma)


def tcre(alpha, beta, gamma, k=K):
    """TCRE in K/EgC."""
    return 1e3 * alpha / (k + beta + alpha * gamma)


def reconstruct(tab: pd.DataFrame, convention: str) -> pd.DataFrame:
    """Reconstruct AF and TCRE from a feedbacks.table() under a convention.

    convention: 'coubgc' (beta_coubgc, gamma_coubgc), 'bgcrad' (beta_bgcrad, gamma_rad),
    'courad' (beta_courad, gamma_rad), 'bgconly_rad' (beta from BGC alone, gamma_rad).
    """
    bmap = {"coubgc": "coubgc", "bgcrad": "bgcrad", "courad": "courad", "bgconly_rad": "bgconly"}
    gmap = {"coubgc": "coubgc", "bgcrad": "rad", "courad": "rad", "bgconly_rad": "rad"}
    b = tab[f"beta_land_{bmap[convention]}"] + tab[f"beta_ocean_{bmap[convention]}"]
    g = tab[f"gamma_land_{gmap[convention]}"] + tab[f"gamma_ocean_{gmap[convention]}"]
    out = pd.DataFrame(index=tab.index)
    out["af_rec"] = af(tab.alpha, b, g)
    out["tcre_rec"] = tcre(tab.alpha, b, g)
    out["af_err"] = out.af_rec - tab.af
    out["tcre_err"] = out.tcre_rec - tab.tcre
    return out
