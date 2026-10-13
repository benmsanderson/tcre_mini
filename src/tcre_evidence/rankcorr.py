"""Iman and Conover (1982) rank-correlation induction: impose a target rank correlation on independent
samples while keeping each marginal exactly."""
from __future__ import annotations

import numpy as np
from scipy import stats


def iman_conover(X: np.ndarray, C: np.ndarray, seed: int = 0) -> np.ndarray:
    """Reorder columns of X (n x d) so their rank correlation approximates C (d x d)."""
    rng = np.random.default_rng(seed)
    n, d = X.shape
    scores = stats.norm.ppf(np.arange(1, n + 1) / (n + 1))
    S = np.column_stack([rng.permutation(scores) for _ in range(d)])
    E = np.corrcoef(S, rowvar=False)
    T = S @ np.linalg.inv(np.linalg.cholesky(E)).T @ np.linalg.cholesky(C).T
    out = np.empty_like(X)
    for j in range(d):
        # column j takes the sorted X values in the rank order of T[:, j]
        out[:, j] = np.sort(X[:, j])[np.argsort(np.argsort(T[:, j]))]
    return out


def fisher_ci(r: float, n: int, level: float = 0.90) -> tuple[float, float]:
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    q = stats.norm.ppf(0.5 + level / 2)
    return float(np.tanh(z - q * se)), float(np.tanh(z + q * se))
