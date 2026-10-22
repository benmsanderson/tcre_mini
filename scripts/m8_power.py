"""M8: how large an ensemble would flag-based process regression need?

Ordinary least squares stand-in for the flag regression (the Bayesian model with weak priors behaves
similarly): y_m = a + sum_p L_mp d_p + e_m, e_m ~ N(0, s^2). The standard error of each d_p is
s * sqrt([(X'X)^-1]_pp), so it depends on the residual scatter s and on the design (how many models carry
each flag, and how flags overlap). Residual s is estimated from the CMIP6 inputs (sheet gamma_L and
log beta_L on the lookup-table flags). Ensembles are enlarged by
  (a) replicating the current design k times (new models with the same flag combinations), and
  (b) a balanced design (flags assigned independently with probability 0.5).
Reports the number of independent models needed for a 2-sigma detection of each effect size.
The calculations are tcre_evidence.process_tests.flag_budget and detection_power.
Writes results/m8_power.csv, results/m8_power.json
"""
import json

import numpy as np
import pandas as pd

from tcre_evidence import process_tests as pt
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS

PROCS = ["nitrogen", "permafrost", "fire", "veg"]
EFFECTS = {"gamma_land": {"nitrogen": 10.0, "permafrost": 30.0, "fire": 5.0, "veg": 5.0},
           "log_beta_land": {"nitrogen": abs(np.log(0.73)), "fire": 0.1, "veg": 0.1}}


def main():
    lt = ud.lookup_table().sel(process=PROCS).to_pandas().T  # models x processes
    sh = ud.cmip_sheet_dataset().drop_dims("process").to_dataframe()
    D = lt.join(sh[["gamma_land", "beta_land"]]).dropna()
    F = D[PROCS]
    out = {"n_models": len(D), "flag_counts": {p: int(F[p].sum()) for p in PROCS},
           "flag_correlation": F.corr().round(2).to_dict()}
    # variance inflation: how well each flag is predicted by the others
    out["vif"] = {p: 1 / (1 - pt.flag_budget(F[p], F.drop(columns=p))["r2"]) for p in PROCS}
    rows = []
    for tgt, effects in EFFECTS.items():
        y = D.gamma_land if tgt == "gamma_land" else np.log(D.beta_land)
        budget = pt.flag_budget(y, F)
        P = pt.detection_power(y, F, effects)
        out[tgt] = dict(residual_sd=budget["resid_sd"], ols_coef=budget["coef"],
                        se_current=dict(zip(PROCS, (budget["resid_sd"] * pt.se_factor(pt._design(F))).round(3).tolist())))
        rows.append(P.assign(target=tgt))
    P = pd.concat(rows).rename(columns={"n_needed_same_mix": "n_needed_same_design"})
    P = P[["target", "process", "effect", "residual_sd", "se_now", "n_needed_same_design", "n_needed_balanced"]]
    P.to_csv(RESULTS / "m8_power.csv", index=False, float_format="%.3g")
    out["power"] = P.to_dict("records")
    json.dump(out, open(RESULTS / "m8_power.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "power"}, indent=1, default=float))
    print(P.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
