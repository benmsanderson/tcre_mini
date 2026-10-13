"""Probabilistic scores for perfect-model tests: central-interval coverage, CRPS, log score."""
from __future__ import annotations

import numpy as np
from scipy import stats


def covered(samples: np.ndarray, truth: float, level: float) -> bool:
    lo, hi = np.percentile(samples, [50 * (1 - level), 50 * (1 + level)])
    return bool(lo <= truth <= hi)


def pit(samples: np.ndarray, truth: float) -> float:
    return float(np.mean(samples < truth))


def crps_samples(samples: np.ndarray, truth: float) -> float:
    """CRPS = E|X - y| - 0.5 E|X - X'| (sorted-sample form, O(n log n))."""
    x = np.sort(np.asarray(samples, float))
    n = x.size
    e1 = np.mean(np.abs(x - truth))
    i = np.arange(1, n + 1)
    e2 = 2 * np.sum((2 * i - n - 1) * x) / n ** 2
    return float(e1 - 0.5 * e2)


def crps_gaussian(mu: float, sd: float, truth: float) -> float:
    z = (truth - mu) / sd
    return float(sd * (z * (2 * stats.norm.cdf(z) - 1) + 2 * stats.norm.pdf(z) - 1 / np.sqrt(np.pi)))


def logscore_samples(samples: np.ndarray, truth: float) -> float:
    """Negative log predictive density at the truth (Gaussian KDE, Scott bandwidth). Lower is better."""
    return float(-np.log(max(stats.gaussian_kde(samples)(truth)[0], 1e-300)))


def summarise(samples: np.ndarray, truth: float) -> dict:
    return dict(truth=float(truth), mean=float(np.mean(samples)), sd=float(np.std(samples)),
                p05=float(np.percentile(samples, 5)), p95=float(np.percentile(samples, 95)),
                pit=pit(samples, truth), cov66=covered(samples, truth, 0.66), cov90=covered(samples, truth, 0.90),
                crps=crps_samples(samples, truth), logscore=logscore_samples(samples, truth))
