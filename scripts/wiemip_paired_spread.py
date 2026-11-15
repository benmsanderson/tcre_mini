"""How much the WIEMIP paired factorials constrain the process effects, under four treatments of them.

The process model's paired evidence (process_tests.build, paired=...) can treat a factorial difference as:
- upstream_rows    the development code's treatment, with every pattern row (UKESM, GFDL, IPSL) entered as a
                   separate value: each observes delta_p (eta_p) directly with error 5 PgC/K (0.1 on eta);
- upstream_models  the same, with the three patterns of a model averaged to one value;
- struct           one value per model, each that model's own effect delta_p + e_mp, with e_mp from the
                   process scatter the regression already uses (sigma_process_struct / sigma_mult_struct);
- estimate         as struct, with the between-model scatter tau_p ~ HalfNormal(30 PgC/K, or 0.5 in log eta)
                   estimated from the factorials themselves.
plus no_paired (flags only) as the reference.

The regression part uses one of two ensembles (--ensemble):
- cmip6   the CMIP6 example tables (data/examples/);
- wiemip  WIEMIP baseline values (pattern mean, results/wiemip/values_patternmean.csv) with the flags of the
          WIEMIP 'process inclusion' sheet (data/wiemip/process_inclusion_sheet.csv, resolved with
          wiemip.read_flag_sheet and written to results/wiemip/flags_resolved.csv). Here the same models give
          both the regression data and the factorials, so the two are not independent evidence.
In both, the flags are nitrogen, permafrost, fire and veg, and the paired values are the WIEMIP
baseline-vs-noX pairs at 2xCO2 (results/wiemip/paired_effects.csv) for those processes (bvoc has no flag).

Writes results/wiemip/paired_spread[_wiemip].csv (one row per quantity, mode and term) and
paired_spread_diag[_wiemip].csv.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse

import pandas as pd
from joblib import Parallel, delayed

from tcre_evidence import process_tests as pt
from tcre_evidence import wiemip as wm
from tcre_evidence.paths import DATA, RESULTS

SEED = 0
OUT = RESULTS / "wiemip"
QUANTS = {"gamma_land_coubgc": (pt.ADDITIVE, "delta_gamma_land_coubgc"),
          "beta_land_coubgc": (pt.MULTIPLICATIVE, "eta_beta_land")}
MODES = {"no_paired": None, "upstream_rows": None, "upstream_models": None, "struct": "struct",
         "estimate": "estimate"}


def paired_rows(pe: pd.DataFrame, column: str) -> dict[str, list[float]]:
    """Every pattern row as its own value (the treatment the per-model average replaces)."""
    d = pe[(pe.background.fillna("") == "") & pe[column].notna()]
    return {p: g[column].tolist() for p, g in d.groupby("process")}


ENSEMBLES = {"cmip6": (DATA / "examples" / "cmip6_land_feedbacks.csv", DATA / "examples" / "cmip6_process_flags.csv"),
             "wiemip": (OUT / "values_patternmean.csv", OUT / "flags_resolved.csv")}


def write_wiemip_flags():
    cl, _ = wm.load()
    flags, conflicts = wm.read_flag_sheet(factorial_flags=wm.factorial_flags(cl))
    if len(conflicts):
        raise ValueError(f"flag sheet contradicts the factorials:\n{conflicts}")
    flags.to_csv(ENSEMBLES["wiemip"][1], float_format="%.0f")


def run_one(ensemble, quantity, mode, sampler):
    kind, column = QUANTS[quantity]
    y, flags, _ = pt.load_ensemble(*ENSEMBLES[ensemble], quantity)
    pe = pd.read_csv(OUT / "paired_effects.csv")
    vals = paired_rows(pe, column) if mode == "upstream_rows" else wm.paired_by_model(pe, column)
    vals = {p: v for p, v in vals.items() if p in flags.columns}
    paired = None if mode == "no_paired" else vals
    table, diag = pt.prior_vs_posterior(y, flags, pt.default_priors(kind, list(flags.columns)), seed=SEED,
                                        sampler=sampler, paired=paired, paired_spread=MODES[mode])
    n = {f"{pt.effect_name(kind)}[{p}]": len(v) for p, v in vals.items()}
    table = table.assign(quantity=quantity, mode=mode, n_paired=table.term.map(n).fillna(0).astype(int)
                         if paired else 0)
    return table, dict(quantity=quantity, mode=mode, **diag)


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("--ensemble", choices=list(ENSEMBLES), default="cmip6")
    a.add_argument("--jobs", type=int, default=1)
    a.add_argument("--quick", action="store_true", help="short chains (2 x 500) for a first look")
    args = a.parse_args()
    sampler = dict(draws=500, tune=500, chains=2) if args.quick else None
    if args.ensemble == "wiemip":
        write_wiemip_flags()
    res = Parallel(n_jobs=args.jobs)(delayed(run_one)(args.ensemble, q, m, sampler) for q in QUANTS for m in MODES)
    T = pd.concat([t for t, _ in res], ignore_index=True)
    D = pd.DataFrame([d for _, d in res])
    suffix = ("" if args.ensemble == "cmip6" else f"_{args.ensemble}") + ("_quick" if args.quick else "")
    cols = ["quantity", "mode", "term", "n_paired", "prior_p05", "prior_median", "prior_p95", "post_p05",
            "post_median", "post_p95", "width_ratio", "shift"]
    T[cols].to_csv(OUT / f"paired_spread{suffix}.csv", index=False, float_format="%.4g")
    D.to_csv(OUT / f"paired_spread_diag{suffix}.csv", index=False)
    show = T[T.term.str.contains(r"\[")]
    print(show[["quantity", "mode", "term", "n_paired", "post_p05", "post_median", "post_p95", "width_ratio"]]
          .round(3).to_string(index=False))
    print(D.to_string(index=False))


if __name__ == "__main__":
    main()
