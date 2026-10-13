"""M0: reproduce the upstream TCRE headline (notebook: mean 1.766, 95 % HDI [0.927, 2.737] K/EgC).

Runs the patched upstream builders on the frozen snapshot for several seeds and
writes results/m0_upstream_reproduction.json.
Usage: python scripts/m0_reproduce_upstream.py [--seeds 0 100 200]
"""
import argparse
import json
import platform

import arviz as az
import numpy as np
import pymc as pm
import pytensor

from tcre_evidence import upstream
from tcre_evidence.paths import RESULTS

TARGET = dict(mean=1.766, hdi95_lo=0.927, hdi95_hi=2.737)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 100, 200])
    args = ap.parse_args()
    runs = []
    for s in args.seeds:
        r = upstream.run(seed=s)
        diag = {k: f.diagnostics for k, f in r["fits"].items()}
        runs.append(dict(seed=s, tcre=r["tcre"], diagnostics=diag))
        print(f"seed {s}: TCRE mean {r['tcre']['mean']:.3f}, 95% HDI [{r['tcre']['hdi95_lo']:.3f}, {r['tcre']['hdi95_hi']:.3f}]")
        for k, d in diag.items():
            print(f"   {k:10s} rhat_max {d['rhat_max']:.4f} ess_bulk_min {d['ess_bulk_min']:.0f} "
                  f"ess_tail_min {d['ess_tail_min']:.0f} divergences {d['divergences']}")
    agg = {k: dict(mean=float(np.mean([r["tcre"][k] for r in runs])),
                   sd_across_seeds=float(np.std([r["tcre"][k] for r in runs], ddof=1)) if len(runs) > 1 else None)
           for k in runs[0]["tcre"]}
    out = dict(
        target=TARGET, runs=runs, across_seeds=agg,
        sampler=upstream.SAMPLER_DEFAULTS,
        versions=dict(python=platform.python_version(), pymc=pm.__version__, arviz=az.__version__,
                      pytensor=pytensor.__version__, numpy=np.__version__),
        notes="Patched upstream builders (src/tcre_evidence/upstream_patched); data from data/upstream_snapshot (2026-10-13).",
    )
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "m0_upstream_reproduction.json").write_text(json.dumps(out, indent=2))
    print("across seeds:", {k: round(v["mean"], 3) for k, v in agg.items()})


if __name__ == "__main__":
    main()
