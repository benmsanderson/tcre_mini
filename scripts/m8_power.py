"""M8: how large an ensemble would flag-based process regression need?

Ordinary least squares stand-in for the flag regression (the Bayesian model with weak priors behaves
similarly): y_m = a + sum_p L_mp d_p + e_m, e_m ~ N(0, s^2). The standard error of each d_p is
s * sqrt([(X'X)^-1]_pp), so it depends on the residual scatter s and on the design (how many models carry
each flag, and how flags overlap). Residual s is estimated from the CMIP6 inputs (sheet gamma_L and
log beta_L on the lookup-table flags). Ensembles are enlarged by
  (a) replicating the current design k times (new models with the same flag combinations), and
  (b) a balanced design (flags assigned independently with probability 0.5).
Reports the number of independent models needed for a 2-sigma detection of each effect size.
Writes results/m8_power.csv, results/m8_power.json
"""
import json

import numpy as np
import pandas as pd

from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS

PROCS = ["nitrogen", "permafrost", "fire", "veg"]


def se_factor(X):
    """sqrt of diag((X'X)^-1) for the process columns (column 0 is the intercept)."""
    XtX = X.T @ X
    if np.linalg.matrix_rank(XtX) < XtX.shape[0]:
        return np.full(X.shape[1] - 1, np.inf)
    return np.sqrt(np.diag(np.linalg.inv(XtX)))[1:]


def main():
    lt = ud.lookup_table().sel(process=PROCS).to_pandas().T  # models x processes
    sh = ud.cmip_sheet_dataset().drop_dims("process").to_dataframe()
    D = lt.join(sh[["gamma_land", "beta_land"]]).dropna()
    X = np.c_[np.ones(len(D)), D[PROCS].values]
    out = {"n_models": len(D), "flag_counts": {p: int(D[p].sum()) for p in PROCS}}
    # flag overlap: correlation between flags and condition number
    out["flag_correlation"] = D[PROCS].corr().round(2).to_dict()
    out["vif"] = {p: float(1 / (1 - np.linalg.lstsq(np.c_[np.ones(len(D)), D[[q for q in PROCS if q != p]].values],
                                                   D[p].values, rcond=None)[1][0] / ((D[p] - D[p].mean()) ** 2).sum()))
                  for p in PROCS}
    targets = {"gamma_land": dict(effects={"nitrogen": 10.0, "permafrost": 30.0, "fire": 5.0, "veg": 5.0}),
               "log_beta_land": dict(effects={"nitrogen": abs(np.log(0.73)), "fire": 0.1, "veg": 0.1})}
    rows = []
    for tgt, spec in targets.items():
        y = D.gamma_land.values if tgt == "gamma_land" else np.log(D.beta_land.values)
        coef, res, *_ = np.linalg.lstsq(X, y, rcond=None)
        s = float(np.sqrt(res[0] / (len(y) - X.shape[1])))
        f0 = se_factor(X)
        out[tgt] = dict(residual_sd=s, ols_coef=dict(zip(PROCS, coef[1:].round(3).tolist())),
                        se_current=dict(zip(PROCS, (s * f0).round(3).tolist())))
        rng = np.random.default_rng(0)
        for p, eff in spec["effects"].items():
            i = PROCS.index(p)
            # (a) replicate current design
            n_rep = next((len(D) * k for k in range(1, 200) if s * se_factor(np.tile(X, (k, 1)))[i] <= eff / 2), np.inf)
            # (b) balanced random design, median over draws
            n_bal = np.inf
            for n in range(8, 2000, 2):
                se = np.median([s * se_factor(np.c_[np.ones(n), rng.integers(0, 2, (n, 4))])[i] for _ in range(50)])
                if se <= eff / 2:
                    n_bal = n
                    break
            rows.append(dict(target=tgt, process=p, effect=eff, residual_sd=s, se_now=float(s * f0[i]),
                             n_needed_same_design=n_rep, n_needed_balanced=n_bal))
    P = pd.DataFrame(rows)
    P.to_csv(RESULTS / "m8_power.csv", index=False, float_format="%.3g")
    out["power"] = P.to_dict("records")
    json.dump(out, open(RESULTS / "m8_power.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "power"}, indent=1, default=float))
    print(P.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
