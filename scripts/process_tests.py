"""Run the process-model tests on any ensemble: a table of values per model and a table of 0/1 flags.

Example (CMIP6 land gamma, the cheap tests only, a few seconds):

    python scripts/process_tests.py --quantity gamma_land_coubgc --kind additive --tests budget power

All five tests on CMIP6 beta_L, with family-wise hold-out in the leave-one-out test:

    python scripts/process_tests.py --quantity beta_land_coubgc --kind multiplicative \
        --tests budget power prior synthetic loo --effects nitrogen=0.73 --groups --jobs 4

For another ensemble (e.g. WIEMIP) pass --values and --flags; see docs/new_ensemble.md for the format.

--effects sets the effect sizes used by the power and synthetic tests, in natural units: delta in the
quantity's units for --kind additive, eta factors (e.g. 0.73) for --kind multiplicative. Processes not
given get the regression estimate (power) or no effect (synthetic).

Outputs go to --out (default results/process_tests/<quantity>/): budget.json, power.csv,
prior_vs_posterior.csv, synthetic_reps.csv, synthetic_summary.csv, loo_scores.csv, loo_summary.csv,
and run.json with the settings and package versions.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # parallel PyMC fits stall with multithreaded BLAS
import argparse
import json
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

from tcre_evidence import process_tests as pt
from tcre_evidence.paths import DATA, RESULTS

TESTS = ["budget", "power", "prior", "synthetic", "loo"]


def parse():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("--values", default=DATA / "examples" / "cmip6_land_feedbacks.csv")
    a.add_argument("--flags", default=DATA / "examples" / "cmip6_process_flags.csv")
    a.add_argument("--quantity", required=True, help="column of --values, e.g. gamma_land_coubgc")
    a.add_argument("--kind", required=True, choices=[pt.ADDITIVE, pt.MULTIPLICATIVE])
    a.add_argument("--tests", nargs="+", default=["budget", "power"], choices=TESTS)
    a.add_argument("--effects", nargs="*", default=[], metavar="PROCESS=SIZE")
    a.add_argument("--nrep", type=int, default=12, help="synthetic replicates")
    a.add_argument("--groups", action="store_true", help="leave out whole families (needs a 'family' column)")
    a.add_argument("--jobs", type=int, default=1, help="parallel fits")
    a.add_argument("--quick", action="store_true", help="short chains (2 x 500) for trying things out")
    a.add_argument("--seed", type=int, default=0)
    a.add_argument("--out", type=Path, default=None)
    return a.parse_args()


def main():
    args = parse()
    out = args.out or RESULTS / "process_tests" / args.quantity
    out.mkdir(parents=True, exist_ok=True)
    y, flags, groups = pt.load_ensemble(args.values, args.flags, args.quantity)
    priors = pt.default_priors(args.kind, list(flags.columns))
    sampler = dict(draws=500, tune=500, chains=2) if args.quick else None
    effects = {k: float(v) for k, v in (e.split("=") for e in args.effects)}
    mult = args.kind == pt.MULTIPLICATIVE
    run = dict(vars(args), out=str(out), n_models=len(y), processes=list(flags.columns), priors=priors,
               versions=dict(pymc=pm.__version__, arviz=az.__version__, numpy=np.__version__))
    print(f"{args.quantity}: {len(y)} models, processes {list(flags.columns)}, {args.kind} model\n")

    # the budget is cheap and sets defaults for the other tests
    budget = pt.flag_budget(np.log(y) if mult else y, flags)
    if "budget" in args.tests:
        json.dump(budget, open(out / "budget.json", "w"), indent=1)
        print(f"[budget] flags explain {budget['r2']:.0%} of the spread in sample (adjusted {budget['adj_r2']:.2f}); "
              f"leave-one-out {budget['loo_q2']:.2f} (negative = worse than the ensemble mean). "
              f"Unexplained scatter {budget['resid_sd']:.3g}{' in log units' if mult else ''}.\n")

    if "power" in args.tests:
        size = {p: abs(np.log(effects[p]) if mult else effects[p]) if p in effects else abs(budget["coef"][p])
                for p in flags.columns}
        P = pt.detection_power(np.log(y) if mult else y, flags, size, seed=args.seed)
        P.to_csv(out / "power.csv", index=False, float_format="%.4g")
        print("[power] models needed for a two-standard-error detection" + (" (effects in log units)" if mult else ""))
        print(P.round(3).to_string(index=False), "\n")

    if "prior" in args.tests:
        T, diag = pt.prior_vs_posterior(y, flags, priors, seed=args.seed, sampler=sampler)
        T.to_csv(out / "prior_vs_posterior.csv", index=False, float_format="%.4g")
        run["prior_diagnostics"] = diag
        print("[prior] width_ratio near 1 and shift near 0: the ensemble has not informed the term "
              "(a narrower posterior can also mean only that the data rule out the prior's tail)")
        print(T.round(3).to_string(index=False), f"\n{diag}\n")

    if "synthetic" in args.tests:
        truth = {p: effects.get(p, 1.0 if mult else 0.0) for p in flags.columns}
        baseline = float(np.exp(np.log(y).mean())) if mult else float(y.mean())
        R, S = pt.synthetic_recovery(flags, priors, baseline, truth, scatter_sd=budget["resid_sd"], nrep=args.nrep,
                                     seed=args.seed, sampler=sampler, n_jobs=args.jobs)
        R.to_csv(out / "synthetic_reps.csv", index=False, float_format="%.4g")
        S.to_csv(out / "synthetic_summary.csv", index=False, float_format="%.4g")
        print(f"[synthetic] {args.nrep} pseudo-ensembles, scatter {budget['resid_sd']:.3g}; "
              "a term whose estimate sits at its prior whatever the truth is not identified")
        print(S.round(3).to_string(index=False), "\n")

    if "loo" in args.tests:
        if args.groups and groups is None:
            raise SystemExit("--groups needs a 'family' column in the flags table")
        R, S = pt.leave_one_out(y, flags, priors, groups=groups if args.groups else None, seed=args.seed,
                                sampler=sampler, n_jobs=args.jobs)
        R.to_csv(out / "loo_scores.csv", index=False, float_format="%.4g")
        S.to_csv(out / "loo_summary.csv", index=False, float_format="%.4g")
        print("[loo] crps_skill < 0: worse than the raw ensemble; bias = mean(estimate - held-out value)")
        print(S.round(3).to_string(index=False), "\n")
        print(f"max r-hat {R.rhat_max.max():.3f}, divergences {int(R.divergences.sum())}")

    json.dump(run, open(out / "run.json", "w"), indent=1, default=str)
    print(f"written to {out}")


if __name__ == "__main__":
    main()
