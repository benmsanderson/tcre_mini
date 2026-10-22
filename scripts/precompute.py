"""Run the slow (sampling) parts of the process-model tests once and cache them for the notebooks.

Default: the CMIP6 example tables, land gamma (additive model) and land beta (multiplicative model).
About an hour with --jobs 4; --quick (short chains) takes a few minutes and is fine for a first look,
but not for reported numbers.

    python scripts/precompute.py --jobs 4

Another ensemble (e.g. WIEMIP), with a name for its cache folder:

    python scripts/precompute.py --name wiemip --values my_values.csv --flags my_flags.csv \
        --quantities gamma_land_coubgc:additive beta_land_coubgc:multiplicative --jobs 4

Writes results/precomputed/<name>/<quantity>/:
    budget.json, power.csv             least-squares tests (also quick to run live)
    prior_vs_posterior.csv, draws.npz  one fit of the real ensemble; prior and posterior draws
    synthetic_<scenario>_{reps,summary}.csv
                                       pseudo-ensembles on the real flags with known effects;
                                       scenario 'effects' has every process active, 'no_effects' none
    loo_{scores,summary}.csv           each model (or family, if the flags table has a 'family'
                                       column) held out in turn
    run.json                           settings, priors, effect sizes, diagnostics, package versions
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # parallel PyMC fits stall with multithreaded BLAS
import argparse
import json
from datetime import date
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

from tcre_evidence import process_tests as pt
from tcre_evidence.paths import DATA, RESULTS

SEED = 0
NREP = 12


def scenario_effects(values: pd.Series, priors: dict) -> dict[str, dict[str, float]]:
    """Known effects for the synthetic test.

    'effects': every process active. Additive: half the spread between models, with the sign of the
    prior (so negative for permafrost). Multiplicative: eta = 0.75.
    'no_effects': no process has any effect.
    """
    procs = list(priors["delta"] if priors["kind"] == pt.ADDITIVE else priors["eta"])
    if priors["kind"] == pt.ADDITIVE:
        size = 0.5 * float(values.std(ddof=1))
        negative = lambda spec: spec["dist"] == "NegLogNormal" or (spec["dist"] == "Normal" and spec["mu"] < 0)
        sign = {p: -1.0 if negative(spec) else 1.0 for p, spec in priors["delta"].items()}
        return {"effects": {p: round(sign[p] * size, 2) for p in procs}, "no_effects": {p: 0.0 for p in procs}}
    return {"effects": {p: 0.75 for p in procs}, "no_effects": {p: 1.0 for p in procs}}


def run_quantity(values_csv, flags_csv, quantity, kind, out: Path, jobs: int, sampler: dict | None, nrep: int):
    out.mkdir(parents=True, exist_ok=True)
    y, flags, groups = pt.load_ensemble(values_csv, flags_csv, quantity)
    priors = pt.default_priors(kind, list(flags.columns))
    mult = kind == pt.MULTIPLICATIVE
    logy = np.log(y) if mult else y
    scen = scenario_effects(y, priors)
    run = dict(quantity=quantity, kind=kind, values=str(values_csv), flags=str(flags_csv), n_models=len(y),
               processes=list(flags.columns), seed=SEED, nrep=nrep, sampler=sampler or "default", priors=priors,
               scenarios=scen, date=str(date.today()),
               versions=dict(pymc=pm.__version__, arviz=az.__version__, numpy=np.__version__))
    print(f"\n=== {quantity} ({kind}, {len(y)} models) -> {out}")

    budget = pt.flag_budget(logy, flags)
    json.dump(budget, open(out / "budget.json", "w"), indent=1)
    size = {p: abs(np.log(e)) if mult else abs(e) for p, e in scen["effects"].items()}
    pt.detection_power(logy, flags, size, seed=SEED).to_csv(out / "power.csv", index=False, float_format="%.4g")

    table, diag, draws = pt.prior_vs_posterior(y, flags, priors, seed=SEED, sampler=sampler, return_draws=True)
    table.to_csv(out / "prior_vs_posterior.csv", index=False, float_format="%.4g")
    rng = np.random.default_rng(SEED)
    thin = lambda a: rng.choice(a[np.isfinite(a)], size=min(2000, np.isfinite(a).sum()), replace=False)
    np.savez_compressed(out / "draws.npz", **{f"{k}|prior": thin(pr) for k, (pr, _) in draws.items()},
                        **{f"{k}|posterior": thin(po) for k, (_, po) in draws.items()})
    run["diagnostics"] = {"prior_vs_posterior": diag}
    print(table.round(3).to_string(index=False))

    baseline = float(np.exp(logy.mean())) if mult else float(y.mean())
    for name, effects in scen.items():
        R, S = pt.synthetic_recovery(flags, priors, baseline, effects, scatter_sd=budget["resid_sd"], nrep=nrep,
                                     seed=SEED, sampler=sampler, n_jobs=jobs)
        R.to_csv(out / f"synthetic_{name}_reps.csv", index=False, float_format="%.4g")
        S.to_csv(out / f"synthetic_{name}_summary.csv", index=False, float_format="%.4g")
        run["diagnostics"][f"synthetic_{name}"] = dict(rhat_max=float(R.rhat_max.max()), divergences=int(R.divergences.sum()))
        print(f"synthetic, {name}:\n{S.round(3).to_string(index=False)}")

    R, S = pt.leave_one_out(y, flags, priors, groups=groups, seed=SEED, sampler=sampler, n_jobs=jobs)
    R.to_csv(out / "loo_scores.csv", index=False, float_format="%.4g")
    S.to_csv(out / "loo_summary.csv", index=False, float_format="%.4g")
    run["diagnostics"]["loo"] = dict(rhat_max=float(R.rhat_max.max()), divergences=int(R.divergences.sum()),
                                     held_out="family" if groups is not None else "model")
    print(f"leave one out:\n{S.round(3).to_string(index=False)}")
    json.dump(run, open(out / "run.json", "w"), indent=1, default=str)


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("--name", default="cmip6")
    a.add_argument("--values", default=DATA / "examples" / "cmip6_land_feedbacks.csv")
    a.add_argument("--flags", default=DATA / "examples" / "cmip6_process_flags.csv")
    a.add_argument("--quantities", nargs="+", metavar="COLUMN:KIND",
                   default=["gamma_land_coubgc:additive", "beta_land_coubgc:multiplicative"])
    a.add_argument("--jobs", type=int, default=1)
    a.add_argument("--nrep", type=int, default=NREP, help="synthetic ensembles per scenario")
    a.add_argument("--quick", action="store_true", help="short chains (2 x 500) for a first look")
    args = a.parse_args()
    sampler = dict(draws=500, tune=500, chains=2) if args.quick else None
    for qk in args.quantities:
        quantity, kind = qk.split(":")
        run_quantity(args.values, args.flags, quantity, kind, RESULTS / "precomputed" / args.name / quantity,
                     args.jobs, sampler, args.nrep)


if __name__ == "__main__":
    main()
