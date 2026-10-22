"""Ensemble-agnostic tests of flag-based process models.

The process model of Marvel et al. (2026, sections 3.5 and 4.3), as used in the TCRE development code,
explains each model's value of a feedback by the processes it represents. With L_mp = 1 when model m
represents process p:

    additive (e.g. gamma_L in PgC/K):
        y_m = U + sum_p L_mp (delta_p + e_mp) + b_common + b_m
    multiplicative (e.g. beta_L in PgC/ppm):
        y_m = U * exp(sum_p L_mp (log eta_p + e_mp)) + b_common + b_m

b_common is a bias shared by all models (its share of the structural variance is rho) and b_m a bias
specific to model m. The estimate for the real world has every process switched on and no bias:
U + sum_p delta_p, or U * prod_p eta_p.

This module builds that model with the development code's own builder (via upstream_patched) and runs
five tests on any ensemble. Inputs are two tables indexed by model name:

    values : pd.Series, one value of the quantity per model
    flags  : pd.DataFrame, one 0/1 column per process

load_ensemble() reads them from CSV; data/examples/ holds the CMIP6 tables used in the technical note.

Tests, cheapest first:

    flag_budget         share of the spread between models that the flags explain, in sample and
                        for each model predicted from the others (seconds)
    detection_power     standard error of each flag effect, and the number of models needed to detect
                        an effect of a given size at two standard errors (seconds)
    prior_vs_posterior  how far the ensemble moves each term away from its prior (one fit, minutes)
    synthetic_recovery  are known effects recovered when pseudo-data are generated on the real flag
                        table? (one fit per replicate)
    leave_one_out       each model in turn plays the real world; is it predicted better than by the
                        raw ensemble? (one fit per model)

Sampling uses PyMC with the settings in upstream.SAMPLER_DEFAULTS. When fits run in parallel (n_jobs > 1),
set OMP_NUM_THREADS=1 (and the OpenBLAS/MKL equivalents) before Python starts; forked workers with
multithreaded BLAS can stall. scripts/process_tests.py does this for you.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pymc as pm
import xarray as xr
from joblib import Parallel, delayed
from scipy import stats

from . import scoring
from . import upstream as up
from .upstream_patched import CMIP_process
from .upstream_patched import model_builders as mb

ADDITIVE, MULTIPLICATIVE = "additive", "multiplicative"
VAR = "y"  # name of the quantity inside the PyMC model; the real-world estimate is the variable "y"


# ----------------------------------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------------------------------

def load_ensemble(values_csv: str | Path, flags_csv: str | Path, quantity: str,
                  group_column: str = "family") -> tuple[pd.Series, pd.DataFrame, pd.Series | None]:
    """Read one quantity and the process flags, keeping models present in both.

    values_csv : CSV with a 'model' column and one column per quantity.
    flags_csv  : CSV with a 'model' column, one 0/1 column per process and, optionally, a column
                 `group_column` naming each model's family (used for leave-one-family-out).
    Returns (values, flags, groups); groups is None when the column is absent.
    """
    V = pd.read_csv(values_csv).set_index("model")
    F = pd.read_csv(flags_csv).set_index("model")
    if quantity not in V:
        raise KeyError(f"{quantity!r} not in {values_csv}; columns are {list(V.columns)}")
    groups = F.pop(group_column) if group_column in F else None
    bad = [c for c in F if not F[c].dropna().isin([0, 1]).all()]
    if bad:
        raise ValueError(f"flag columns must be 0/1; check {bad}")
    y = V[quantity].dropna()
    models = [m for m in y.index if m in F.index and F.loc[m].notna().all()]
    dropped = sorted((set(V.index) | set(F.index)) - set(models))
    if dropped:
        print(f"load_ensemble: using {len(models)} models; not in both tables or incomplete: {dropped}")
    flags = F.loc[models].astype(int)
    unused = [c for c in flags if flags[c].nunique() < 2]
    if unused:
        print(f"load_ensemble: flags {unused} are the same for every model and cannot be estimated")
    return y.loc[models].astype(float), flags, (groups.loc[models] if groups is not None else None)


def default_priors(kind: str, processes: list[str]) -> dict:
    """The development code's priors (config/priors.py), restricted to `processes`.

    additive uses the gamma_land priors and multiplicative the beta_land priors. A process the
    development code does not have (for example 'demography') gets a wide prior centred on no effect:
    delta ~ N(0, 10) PgC/K, or log eta ~ N(0, 0.3).
    """
    if kind == ADDITIVE:
        src = up.priors.gamma_land
        return dict(kind=kind, baseline=src["baseline"],
                    delta={p: src["delta"].get(p, dict(dist="Normal", mu=0.0, sigma=10.0)) for p in processes},
                    sigma_process_struct={p: src["sigma_process_struct"].get(p, 10.0) for p in processes},
                    sigma_struct=src["sigma_struct"], likelihood_noise=src["likelihood_noise"])
    if kind == MULTIPLICATIVE:
        src = up.priors.beta_land
        return dict(kind=kind, baseline=src["baseline"],
                    eta={p: src["eta"].get(p, dict(mu=0.0, sigma=0.3)) for p in processes},
                    sigma_mult_struct={p: src["sigma_mult_struct"].get(p, 0.01) for p in processes},
                    sigma_struct=src["sigma_struct"], likelihood_noise=src["likelihood_noise"])
    raise ValueError(f"kind must be {ADDITIVE!r} or {MULTIPLICATIVE!r}")


# ----------------------------------------------------------------------------------------------------
# The Bayesian process model
# ----------------------------------------------------------------------------------------------------

def build(values: pd.Series, flags: pd.DataFrame, priors: dict,
          paired: dict[str, list[float]] | None = None, paired_sigma: float | None = None) -> pm.Model:
    """PyMC process model for one quantity, built with the development code's builder.

    paired : optional evidence from paired process experiments (the same model with a process on and
             off), as {process: [values]}; values are delta (additive) or eta (multiplicative). In the
             development code these are, for example, beta_land nitrogen eta = [0.69, 0.77].
    """
    procs = list(flags.columns)
    y = xr.DataArray(values.values, coords=dict(model=values.index.values), dims="model")
    lt = xr.DataArray(flags.T.values.astype(float), coords=dict(process=procs, model=values.index.values),
                      dims=("process", "model"))
    common = dict(noprocess_prior=mb._noprocess_prior(priors["baseline"], VAR), sigma_struct=priors["sigma_struct"],
                  likelihood_noise=priors["likelihood_noise"], var_name=VAR)
    if priors["kind"] == ADDITIVE:
        model = CMIP_process.build_vectorized_process_model(
            y, lt, additive_processes=procs, delta_priors=mb._delta_priors(priors["delta"]),
            sigma_process_struct=priors["sigma_process_struct"], **common)
    else:
        model = CMIP_process.build_vectorized_process_model(
            y, lt, multiplicative_processes=procs, eta_priors=priors["eta"],
            sigma_mult_struct=priors["sigma_mult_struct"], **common)
    if paired:
        sigma = paired_sigma if paired_sigma is not None else (5.0 if priors["kind"] == ADDITIVE else 0.1)
        with model:
            for p, vals in paired.items():
                CMIP_process.add_process_information(model, p, np.asarray(vals, float), process_sigma=sigma,
                                                     process_iteration=1, var_name=VAR)
    return model


def fit(model: pm.Model, seed: int = 0, sampler: dict | None = None):
    """Sample the posterior. Returns (idata, diagnostics); check diagnostics['rhat_max'] < 1.01 and divergences.

    Chains run one after another (cores=1): parallelism comes from n_jobs in the tests, and PyMC's own
    multiprocess sampling in the parent process breaks later joblib workers.
    """
    idata = up._sample(model, seed, {"cores": 1, **(sampler or {})}, progressbar=False)
    names = [VAR] + [v for v in (f"{VAR}_delta", f"{VAR}_eta") if v in idata.posterior]
    return idata, up.diagnostics(idata, names)


def _draws(idata, var: str) -> np.ndarray:
    x = idata.posterior[var].values
    return x.reshape(x.shape[0] * x.shape[1], *x.shape[2:])


def effect_name(kind: str) -> str:
    """Name of the process-effect variable in the PyMC model: y_delta (additive) or y_eta (multiplicative)."""
    return f"{VAR}_delta" if kind == ADDITIVE else f"{VAR}_eta"


def own_flags_predictive(idata, flag_row: pd.Series, priors: dict, seed: int = 0) -> np.ndarray:
    """Posterior predictive for a new model with the given flags (the model's own generative structure):
    fresh process scatter e, model bias b_m ~ N(0, sigma_m) and likelihood noise, plus the shared bias."""
    rng = np.random.default_rng(seed)
    L = flag_row.values.astype(float)
    procs = list(flag_row.index)
    U = _draws(idata, f"{VAR}_unscaled")
    n = U.shape[0]
    if priors["kind"] == ADDITIVE:
        d = _draws(idata, f"{VAR}_delta")
        e = rng.normal(0, [priors["sigma_process_struct"][p] for p in procs], d.shape)
        core = U + ((d + e) * L).sum(1)
    else:
        log_eta = _draws(idata, f"{VAR}_log_eta")
        e = rng.normal(0, [priors["sigma_mult_struct"][p] for p in procs], log_eta.shape)
        core = U * np.exp(((log_eta + e) * L).sum(1))
    return (core + _draws(idata, f"{VAR}_bias_common") + rng.normal(0, 1, n) * _draws(idata, f"{VAR}_sigma_m")
            + rng.normal(0, priors["likelihood_noise"], n))


# ----------------------------------------------------------------------------------------------------
# Test 1: how much of the spread do the flags explain?
# ----------------------------------------------------------------------------------------------------

def _design(flags: pd.DataFrame) -> np.ndarray:
    return np.column_stack([np.ones(len(flags)), flags.values.astype(float)])


def flag_budget(values: pd.Series, flags: pd.DataFrame) -> dict:
    """Least-squares regression of the quantity on the flags.

    r2       share of the variance between models explained in sample
    adj_r2   the same, adjusted for the number of flags (negative: worse than chance)
    loo_q2   share explained when each model is predicted from the others' regression;
             negative means the flags predict worse than the ensemble mean
    resid_sd scatter the flags leave unexplained (in the quantity's units)
    """
    y, X = values.values.astype(float), _design(flags)
    n, p = X.shape
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    ss = ((y - y.mean()) ** 2).sum()
    r2 = 1 - (r ** 2).sum() / ss
    loo = []
    for i in range(n):
        keep = np.arange(n) != i
        bi, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
        loo.append(y[i] - X[i] @ bi)
    return dict(n=n, sd=float(y.std(ddof=1)), r2=float(r2), adj_r2=float(1 - (1 - r2) * (n - 1) / (n - p)),
                loo_q2=float(1 - np.sum(np.square(loo)) / ss), resid_sd=float(np.sqrt((r ** 2).sum() / (n - p))),
                coef=dict(zip(flags.columns, b[1:].round(4).tolist())))


# ----------------------------------------------------------------------------------------------------
# Test 2: how many models would be needed?
# ----------------------------------------------------------------------------------------------------

def se_factor(X: np.ndarray) -> np.ndarray:
    """sqrt(diag((X'X)^-1)) for the flag columns; times the residual sd gives each effect's standard error."""
    XtX = X.T @ X
    if np.linalg.matrix_rank(XtX) < XtX.shape[0]:
        return np.full(X.shape[1] - 1, np.inf)
    return np.sqrt(np.diag(np.linalg.inv(XtX)))[1:]


def detection_power(values: pd.Series, flags: pd.DataFrame, effects: dict[str, float],
                    seed: int = 0, max_models: int = 2000) -> pd.DataFrame:
    """Standard error of each flag effect now, and models needed for a two-standard-error detection.

    A least-squares stand-in for the Bayesian model (which behaves similarly with weak priors). effects
    are in the quantity's units (use log values and log effects for a multiplicative quantity). Two
    ways of growing the ensemble: replicating the current flag mix, and a balanced design in which
    each flag is set independently with probability 0.5 (median over 50 random designs per size, so the
    balanced numbers depend on `seed` and on the order of `effects`). Assumes the new models are independent.
    """
    X = _design(flags)
    s = flag_budget(values, flags)["resid_sd"]
    f0 = se_factor(X)
    rng = np.random.default_rng(seed)
    rows = []
    for p, eff in effects.items():
        i = list(flags.columns).index(p)
        n_same = next((len(X) * k for k in range(1, max_models // len(X) + 1)
                       if s * se_factor(np.tile(X, (k, 1)))[i] <= abs(eff) / 2), np.inf)
        n_bal = np.inf
        for n in range(8, max_models, 2):
            se = np.median([s * se_factor(np.c_[np.ones(n), rng.integers(0, 2, (n, X.shape[1] - 1))])[i]
                            for _ in range(50)])
            if se <= abs(eff) / 2:
                n_bal = n
                break
        rows.append(dict(process=p, effect=eff, residual_sd=s, se_now=float(s * f0[i]), n_models_now=len(X),
                         n_needed_same_mix=n_same, n_needed_balanced=n_bal))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------------
# Test 3: has the ensemble updated the priors?
# ----------------------------------------------------------------------------------------------------

def prior_vs_posterior(values: pd.Series, flags: pd.DataFrame, priors: dict, seed: int = 0,
                       sampler: dict | None = None, paired: dict | None = None, return_draws: bool = False):
    """Prior and posterior of each process effect, the real-world value and the shared-bias terms.

    width_ratio = width of the posterior 90 % interval / width of the prior's: near 1 means the ensemble
                  has not informed the term.
    shift       = (posterior median - prior median) / width of the prior 90 % interval.
    Quantiles are used because several priors are skewed (permafrost, nitrogen) and their standard
    deviations are set by long tails. A narrower posterior can also mean only that the data rule out
    the prior's tail; synthetic_recovery is the direct test of whether a term is identified.

    Returns (table, diagnostics), plus {term: (prior draws, posterior draws)} if return_draws.
    """
    model = build(values, flags, priors, paired=paired)
    with model:
        prior = pm.sample_prior_predictive(draws=4000, random_seed=seed).prior
    idata, diag = fit(model, seed=seed, sampler=sampler)
    eff = effect_name(priors["kind"])
    draws = {f"{eff}[{p}]": (prior[eff].values[..., i].ravel(), _draws(idata, eff)[:, i])
             for i, p in enumerate(flags.columns)}
    for v in (VAR, f"{VAR}_sigma_struct", f"{VAR}_rho", f"{VAR}_sigma_common"):
        if v in idata.posterior:
            draws[v] = (prior[v].values.ravel(), _draws(idata, v))
    table = pd.DataFrame([_compare(term, pr, po) for term, (pr, po) in draws.items()])
    return (table, diag, draws) if return_draws else (table, diag)


def _compare(term, pr, po) -> dict:
    q = lambda x: np.percentile(x[np.isfinite(x)], [5, 50, 95])
    (a05, a50, a95), (b05, b50, b95) = q(pr), q(po)
    return dict(term=term, prior_p05=a05, prior_median=a50, prior_p95=a95, post_p05=b05, post_median=b50,
                post_p95=b95, width_ratio=(b95 - b05) / (a95 - a05), shift=(b50 - a50) / (a95 - a05))


# ----------------------------------------------------------------------------------------------------
# Test 4: are known effects recovered?
# ----------------------------------------------------------------------------------------------------

def synthetic_ensemble(flags: pd.DataFrame, kind: str, baseline: float, effects: dict[str, float],
                       scatter_sd: float, common_sd: float = 0.0, rng=None) -> pd.Series:
    """Pseudo-data on the real flag table with known effects.

    additive:       y_m = baseline + sum_p L_mp effects[p] + c + e_m
    multiplicative: y_m = baseline * prod_p effects[p]^L_mp * exp(c + e_m)   (effects are eta factors)
    with c ~ N(0, common_sd) shared by all models and e_m ~ N(0, scatter_sd); for multiplicative,
    scatter_sd and common_sd are in log units.
    """
    rng = rng or np.random.default_rng(0)
    L = flags.values.astype(float)
    eff = np.array([effects.get(p, 0.0 if kind == ADDITIVE else 1.0) for p in flags.columns])
    noise = rng.normal(0, common_sd) + rng.normal(0, scatter_sd, len(flags))
    if kind == ADDITIVE:
        y = baseline + L @ eff + noise
    else:
        y = baseline * np.exp(L @ np.log(eff) + noise)
    return pd.Series(y, index=flags.index)


def _recover_one(flags, priors, baseline, effects, scatter_sd, common_sd, rep, seed, sampler):
    kind = priors["kind"]
    y = synthetic_ensemble(flags, kind, baseline, effects, scatter_sd, common_sd,
                           rng=np.random.default_rng(10_000 * (seed + 1) + rep))
    idata, diag = fit(build(y, flags, priors), seed=seed + rep, sampler=sampler)
    est = _draws(idata, effect_name(kind))
    truth_rw = (baseline + sum(effects.values())) if kind == ADDITIVE else baseline * np.prod(list(effects.values()))
    rows = []
    for i, p in enumerate(flags.columns):
        t = effects.get(p, 0.0 if kind == ADDITIVE else 1.0)
        rows.append(dict(rep=rep, term=p, truth=t, **_estimate(est[:, i], t), **diag))
    rows.append(dict(rep=rep, term="real_world", truth=truth_rw, **_estimate(_draws(idata, VAR), truth_rw), **diag))
    return rows


def _estimate(s, truth) -> dict:
    lo, hi = np.percentile(s, [5, 95])
    return dict(post_mean=float(s.mean()), post_sd=float(s.std()), cov90=bool(lo <= truth <= hi))


def synthetic_recovery(flags: pd.DataFrame, priors: dict, baseline: float, effects: dict[str, float],
                       scatter_sd: float, common_sd: float = 0.0, nrep: int = 12, seed: int = 0,
                       sampler: dict | None = None, n_jobs: int = 1) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit the process model to `nrep` synthetic ensembles with known effects.

    Choose scatter_sd from flag_budget(...)['resid_sd'] on the real data so that the pseudo-data are as
    noisy as the real ensemble. Returns (per-replicate rows, summary by term). In the summary,
    post_sd close to the prior sd and an estimate near the prior mean, whatever the truth, mean the
    term is not identified.
    """
    res = Parallel(n_jobs=n_jobs)(delayed(_recover_one)(flags, priors, baseline, effects, scatter_sd, common_sd,
                                                        r, seed, sampler) for r in range(nrep))
    R = pd.DataFrame([x for r in res for x in r])
    S = (R.groupby("term", sort=False)
         .agg(truth=("truth", "first"), mean_estimate=("post_mean", "mean"), spread_between_reps=("post_mean", "std"),
              post_sd=("post_sd", "mean"), cov90=("cov90", "mean"), rhat_max=("rhat_max", "max"),
              divergences=("divergences", "sum")).reset_index())
    return R, S


# ----------------------------------------------------------------------------------------------------
# Test 5: is each held-out model predicted?
# ----------------------------------------------------------------------------------------------------

def _loo_one(values, flags, groups, priors, m, seed, sampler):
    out = groups.index[groups == groups[m]] if groups is not None else [m]
    train = [x for x in values.index if x not in out]
    idata, diag = fit(build(values[train], flags.loc[train], priors), seed=seed, sampler=sampler)
    truth = float(values[m])
    yt = values[train].values
    rng = np.random.default_rng(seed + 1)
    raw = yt.mean() + yt.std(ddof=1) * np.sqrt(1 + 1 / len(yt)) * stats.t.rvs(len(yt) - 1, size=4000,
                                                                               random_state=rng)
    preds = {"all_processes_on": _draws(idata, VAR),
             "own_flags": own_flags_predictive(idata, flags.loc[m], priors, seed=seed + 2),
             "raw_ensemble": raw}
    return [dict(model=m, method=k, n_train=len(train), **scoring.summarise(v, truth), **diag) for k, v in preds.items()]


def leave_one_out(values: pd.Series, flags: pd.DataFrame, priors: dict, groups: pd.Series | None = None,
                  seed: int = 0, sampler: dict | None = None, n_jobs: int = 1) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Each model in turn plays the real world; it (or its whole family, if `groups` is given) is removed
    from training, and three predictions of its value are scored:

    all_processes_on  the model's real-world estimate (what the configuration reports for the real world)
    own_flags         predictive for a new model with the held-out model's own flags
    raw_ensemble      Student-t predictive from the training models' mean and spread

    Returns (per-model scores, summary). In the summary, crps_skill = 1 - CRPS / CRPS(raw_ensemble):
    negative means worse than the raw ensemble. cov66 and cov90 are the shares of held-out values inside
    the central 66 % and 90 % intervals.
    """
    res = Parallel(n_jobs=n_jobs)(delayed(_loo_one)(values, flags, groups, priors, m, seed, sampler)
                                  for m in values.index)
    R = pd.DataFrame([x for r in res for x in r])
    R["error"] = R["mean"] - R["truth"]
    S = R.groupby("method").agg(crps=("crps", "mean"), bias=("error", "mean"), cov66=("cov66", "mean"),
                                cov90=("cov90", "mean"), mean_sd=("sd", "mean")).reset_index()
    S["crps_skill"] = 1 - S.crps / float(S.loc[S.method == "raw_ensemble", "crps"].iloc[0])
    return R, S


# ----------------------------------------------------------------------------------------------------
# Summary across tests
# ----------------------------------------------------------------------------------------------------

def evidence_table(prior_table: pd.DataFrame, synth_effects: pd.DataFrame, synth_none: pd.DataFrame,
                   power: pd.DataFrame) -> pd.DataFrame:
    """One row per process, combining the tests:

    models_needed  two-standard-error detection with the current flag mix (power test)
    tracking       how far the estimate follows the truth between the two synthetic scenarios, averaged
                   over replicates: (estimate with effects - estimate without) / (true difference).
                   1 = follows the truth, 0 = ignores it (the answer comes from the prior)
    sd_over_effect posterior sd of one fit / the true difference. Tracking is an average over many
                   pseudo-ensembles; this says whether a single ensemble resolves the effect. Above about
                   0.5 the effect is not detected at two standard deviations
    width_ratio    posterior over prior 90 % width on the real ensemble (1 = not informed)
    """
    a, b = synth_effects.set_index("term"), synth_none.set_index("term")
    procs = [p for p in a.index if p != "real_world"]
    width = prior_table.assign(process=prior_table.term.str.extract(r"\[(.*)\]")[0]).set_index("process")
    P = power.set_index("process")
    rows = []
    for p in procs:
        dt = a.loc[p, "truth"] - b.loc[p, "truth"]
        rows.append(dict(process=p, models_available=int(P.loc[p, "n_models_now"]),
                         models_needed=P.loc[p, "n_needed_same_mix"],
                         tracking=(a.loc[p, "mean_estimate"] - b.loc[p, "mean_estimate"]) / dt if dt else np.nan,
                         sd_over_effect=0.5 * (a.loc[p, "post_sd"] + b.loc[p, "post_sd"]) / abs(dt) if dt else np.nan,
                         width_ratio=width.loc[p, "width_ratio"]))
    return pd.DataFrame(rows)
