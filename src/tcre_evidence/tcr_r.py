"""Workstream B: TCRE = TCR x R, with R constrained by a historical observable H (HANDOVER §5).

Cross-model model: z_m = (TCR_m, R_m, H_m) ~ MVN(mu, Sigma). Small-ensemble uncertainty in (mu, Sigma)
is carried by the Normal-inverse-Wishart (Jeffreys) predictive, a multivariate Student-t with n - d dof
(Ribes et al. 2021 style kriging, closed form). The real world z* is a new draw from that predictive.

Two ways to use the assessed TCR:
  'substitute'  TCR* drawn from the assessed distribution; R* drawn from p(R* | TCR*, H_obs) under the
                ensemble predictive (the TCR-R covariance is kept, the ensemble's TCR marginal is not used);
  'unconstrained'  as 'substitute' without H (R* | TCR* only);
  'ensemble'    R* from the ensemble marginal, independent of TCR (no covariance).
H_obs carries observational error sigma_obs, added to the H-H block.
R is in EgC^-1 (TCRE / TCR), so TCRE = TCR x R in K/EgC.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def mvt_conditional(Z: np.ndarray, idx_t: list[int], idx_c: list[int], x_c: np.ndarray, extra_cov_c=None):
    """Predictive of target block given conditioning block for a new draw from the NIW predictive.

    Z: (n, d) ensemble. Returns (loc, scale_matrix, dof) of the conditional multivariate t.
    """
    n, d = Z.shape
    mu = Z.mean(0)
    S = np.cov(Z, rowvar=False) * (n - 1)            # scatter matrix
    nu = n - d                                       # posterior dof (Jeffreys prior)
    Psi = S * (n + 1) / (n * nu)                     # predictive scale (marginal t with nu dof)
    if extra_cov_c is not None:
        Psi = Psi.copy()
        Psi[np.ix_(idx_c, idx_c)] += np.asarray(extra_cov_c) * 1.0
    A = Psi[np.ix_(idx_t, idx_t)]
    B = Psi[np.ix_(idx_t, idx_c)]
    C = Psi[np.ix_(idx_c, idx_c)]
    Ci = np.linalg.inv(C)
    dx = np.asarray(x_c) - mu[idx_c]
    loc = mu[idx_t] + B @ Ci @ dx
    nu_c = nu + len(idx_c)
    scale = (nu + dx @ Ci @ dx) / nu_c * (A - B @ Ci @ B.T)
    return loc, scale, nu_c


def sample_t(loc, scale, dof, size, rng):
    loc, scale = np.atleast_1d(loc), np.atleast_2d(scale)
    z = rng.multivariate_normal(np.zeros(len(loc)), scale, size)
    w = rng.chisquare(dof, size) / dof
    return loc + z / np.sqrt(w)[:, None]


def predict_R(tcr: np.ndarray, R: np.ndarray, H: np.ndarray | None, tcr_star: np.ndarray, h_obs: float | None = None,
              sigma_obs: float = 0.0, mode: str = "substitute", seed: int = 0) -> np.ndarray:
    """Samples of R* for each TCR* sample (vectorised by binning TCR* draws is unnecessary: closed form per draw)."""
    rng = np.random.default_rng(seed)
    tcr_star = np.asarray(tcr_star)
    N = tcr_star.size
    if mode == "ensemble":
        n = R.size
        return R.mean() + R.std(ddof=1) * np.sqrt(1 + 1 / n) * stats.t.rvs(n - 1, size=N, random_state=rng)
    if mode == "unconstrained" or H is None:
        Z = np.c_[tcr, R]
        idx_c = [0]
        extra = None
    else:
        Z = np.c_[tcr, R, H]
        idx_c = [0, 2]
        extra = np.diag([0.0, sigma_obs ** 2])
    out = np.empty(N)
    # the conditional loc is linear in TCR*, the scale depends on TCR* only through the Mahalanobis term
    for i, t in enumerate(tcr_star):
        xc = [t] if len(idx_c) == 1 else [t, h_obs]
        loc, scale, dof = mvt_conditional(Z, [1], idx_c, np.array(xc), extra)
        out[i] = loc[0] + np.sqrt(scale[0, 0]) * rng.standard_t(dof)
    return out


def predict_R_given_tcr(tcr, R, H, tcr_m, h_m, sigma_obs, mode, n=4000, seed=0):
    """Perfect-model predictive of R for a held-out model with known TCR_m and observable H_m."""
    rng = np.random.default_rng(seed)
    if mode == "ensemble":
        k = R.size
        return R.mean() + R.std(ddof=1) * np.sqrt(1 + 1 / k) * stats.t.rvs(k - 1, size=n, random_state=rng)
    if mode == "unconstrained" or H is None:
        loc, scale, dof = mvt_conditional(np.c_[tcr, R], [1], [0], np.array([tcr_m]))
    elif mode == "H_only":
        loc, scale, dof = mvt_conditional(np.c_[R, H], [0], [1], np.array([h_m]), np.array([[sigma_obs ** 2]]))
    else:
        loc, scale, dof = mvt_conditional(np.c_[tcr, R, H], [1], [0, 2], np.array([tcr_m, h_m]),
                                          np.diag([0.0, sigma_obs ** 2]))
    return sample_t(loc, scale, dof, n, rng)[:, 0]
