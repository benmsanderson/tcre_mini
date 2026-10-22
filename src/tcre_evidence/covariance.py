"""Recombination of per-model feedback components under different sampling schemes (HANDOVER §2.3, A3).

Each scheme draws model indices with replacement; components that share an index
keep their cross-model covariance, components with separate indices lose it.
"""
from __future__ import annotations

import numpy as np

from .upstream_data import K

COMPONENTS = ("alpha", "beta_land", "gamma_land", "beta_ocean", "gamma_ocean")

# Groups of components drawn with a common model index.
SCHEMES = {
    "paired": [COMPONENTS],
    # upstream notebook: beta_L, gamma_L, alpha sampled independently; ocean pair joint
    "upstream": [("alpha",), ("beta_land",), ("gamma_land",), ("beta_ocean", "gamma_ocean")],
    "independent": [(c,) for c in COMPONENTS],
}


def tcre_from_components(alpha, beta_land, gamma_land, beta_ocean, gamma_ocean, k=K):
    """JF2020 identity. Returns (AF, TCRE in K/EgC); AF is the dimensionless airborne fraction."""
    denom = k + beta_land + beta_ocean + alpha * (gamma_land + gamma_ocean)
    return k / denom, 1e3 * alpha / denom


def resample(values: dict, scheme: str, n: int = 400_000, seed: int = 0) -> dict:
    """Bootstrap-recombine per-model component arrays under a sampling scheme."""
    rng = np.random.default_rng(seed)
    m = len(values["alpha"])
    out = {}
    for group in SCHEMES[scheme]:
        idx = rng.integers(0, m, n)
        for c in group:
            out[c] = np.asarray(values[c])[idx]
    out["af"], out["tcre"] = tcre_from_components(*(out[c] for c in COMPONENTS))
    return out


def summary(draws: dict) -> dict:
    """Spread of the airborne fraction and TCRE summary statistics from resampled draws."""
    t = draws["tcre"]
    return dict(sd_af=float(np.std(draws["af"])), tcre_mean=float(t.mean()), tcre_sd=float(t.std()),
                tcre_p05=float(np.percentile(t, 5)), tcre_p95=float(np.percentile(t, 95)))


def pair_correlations(v: dict) -> dict:
    """Correlations across models between the pairs of components that matter for TCRE."""
    a, bL, gL, bO, gO = (np.asarray(v[c]) for c in COMPONENTS)
    c = lambda x, y: float(np.corrcoef(x, y)[0, 1])
    return {"beta_L,gamma_L": c(bL, gL), "beta,alpha*gamma": c(bL + bO, a * (gL + gO)),
            "alpha,alpha*gamma": c(a, a * (gL + gO)), "beta_O,gamma_O": c(bO, gO), "alpha,beta_O": c(a, bO)}
