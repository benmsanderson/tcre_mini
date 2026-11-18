"""Synthetic recovery and leave one out on WIEMIP, with the paired factorials treated in different ways.

The ensemble is WIEMIP (baseline values, pattern mean; flags from results/wiemip/flags_resolved.csv, written
by wiemip_paired_spread.py --ensemble wiemip). Paired factorials are one value per model
(wiemip.paired_by_model). Treatments (process_tests.build):
- no_paired   flags only;
- upstream    the development code's: each factorial observes the mean effect with error 5 PgC/K (0.1 on eta);
- estimate    each factorial is one model's own effect; the between-model scatter tau_p is estimated.

synthetic  Pseudo-ensembles on the real flag table with known effects (the WIEMIP paired means), residual
           scatter from the flag budget of the real values, and pseudo factorials with the real number of
           models per process. Two scatter scenarios for the factorials:
             observed  tau_p = the between-model sd of the real factorial effects;
             small     tau_p = 0 (the factorials agree up to measurement error).
           A treatment is calibrated if cov90 (share of replicates whose 90 % interval holds the truth) is
           near 0.9.
loo        Each WIEMIP model in turn is held out with its own factorials; its baseline value is predicted
           with the held-out model's own flags (own_flags) and the all-processes estimate, and scored
           against the raw ensemble (crps_skill > 0 is better than the raw ensemble).

Writes results/wiemip/paired_tests_{synthetic,loo}[_quick].csv and *_summary[_quick].csv.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse

import numpy as np
import pandas as pd

from tcre_evidence import process_tests as pt
from tcre_evidence import wiemip as wm
from tcre_evidence.paths import RESULTS

SEED = 0
OUT = RESULTS / "wiemip"
QUANTS = {"gamma_land_coubgc": (pt.ADDITIVE, "delta_gamma_land_coubgc"),
          "beta_land_coubgc": (pt.MULTIPLICATIVE, "eta_beta_land")}
MEAS_SD = {pt.ADDITIVE: 2.0, pt.MULTIPLICATIVE: 0.02}  # build()'s default measurement error in "estimate" mode
TREATMENTS = {"no_paired": None, "upstream": {}, "estimate": {"paired_spread": "estimate"}}


def ensemble(quantity):
    kind, column = QUANTS[quantity]
    y, flags, _ = pt.load_ensemble(OUT / "values_patternmean.csv", OUT / "flags_resolved.csv", quantity)
    pe = pd.read_csv(OUT / "paired_effects.csv")
    paired = {p: s for p, s in wm.paired_by_model(pe, column, as_series=True).items() if p in flags.columns}
    return kind, y, flags, paired


def synthetic(quantity, nrep, sampler, jobs):
    kind, y, flags, paired = ensemble(quantity)
    to_eff = (lambda s: s.mean()) if kind == pt.ADDITIVE else (lambda s: float(np.exp(np.log(s).mean())))
    effects = {p: to_eff(s) for p, s in paired.items()}
    spread = {p: (s.std(ddof=1) if kind == pt.ADDITIVE else np.log(s).std(ddof=1)) if len(s) > 1 else 0.0
              for p, s in paired.items()}
    logy = np.log(y) if kind == pt.MULTIPLICATIVE else y   # flag_budget's scatter in the synthetic units
    scatter = pt.flag_budget(logy, flags)["resid_sd"]
    baseline = float(y[(flags == 0).all(axis=1)].mean()) if (flags == 0).all(axis=1).any() else float(y.median())
    priors = pt.default_priors(kind, list(flags.columns))
    rows = []
    for scen, taus in {"observed": spread, "small": {p: 0.0 for p in spread}}.items():
        design = {p: (len(paired[p]), taus[p]) for p in paired}
        for name, kw in TREATMENTS.items():
            if name == "no_paired" and scen == "small":
                continue                                      # no factorials, so the scenario does not apply
            R, S = pt.synthetic_recovery(flags, priors, baseline, effects, scatter_sd=scatter, nrep=nrep, seed=SEED,
                                         sampler=sampler, n_jobs=jobs,
                                         paired_design=None if kw is None else design,
                                         paired_meas_sd=MEAS_SD[kind], paired_kw=kw)
            print(f"{quantity} {scen} {name}: max r-hat {R.rhat_max.max():.3f}, divergences {R.divergences.sum()}")
            extra = dict(quantity=quantity, scenario=scen, treatment=name, scatter_sd=scatter, baseline=baseline)
            rows.append((R.assign(**extra), S.assign(**extra, tau=S.term.map(taus), n_paired=S.term.map(
                {p: len(v) for p, v in paired.items()}))))
    return rows


def loo(quantity, sampler, jobs):
    kind, y, flags, paired = ensemble(quantity)
    priors = pt.default_priors(kind, list(flags.columns))
    rows = []
    for name, kw in TREATMENTS.items():
        R, S = pt.leave_one_out(y, flags, priors, seed=SEED, sampler=sampler, n_jobs=jobs,
                                paired=None if kw is None else paired, paired_kw=kw)
        print(f"{quantity} loo {name}: max r-hat {R.rhat_max.max():.3f}, divergences {R.divergences.sum()}")
        rows.append((R.assign(quantity=quantity, treatment=name), S.assign(quantity=quantity, treatment=name)))
    return rows


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("test", choices=["synthetic", "loo"])
    a.add_argument("--jobs", type=int, default=1)
    a.add_argument("--nrep", type=int, default=16)
    a.add_argument("--quick", action="store_true", help="short chains (2 x 500) and 4 replicates for a first look")
    args = a.parse_args()
    sampler = dict(draws=500, tune=500, chains=2) if args.quick else None
    nrep = 4 if args.quick else args.nrep
    suffix = "_quick" if args.quick else ""
    rows = []
    for q in QUANTS:  # written after each quantity, so a run stopped part-way keeps what it finished
        rows += synthetic(q, nrep, sampler, args.jobs) if args.test == "synthetic" else loo(q, sampler, args.jobs)
        R = pd.concat([r for r, _ in rows], ignore_index=True)
        S = pd.concat([s for _, s in rows], ignore_index=True)
        R.to_csv(OUT / f"paired_tests_{args.test}{suffix}.csv", index=False, float_format="%.5g")
        S.to_csv(OUT / f"paired_tests_{args.test}_summary{suffix}.csv", index=False, float_format="%.4g")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
