"""M7: does the nitrogen-flag difference in beta_L carry over to the historical land sink?

Compares models with and without the upstream nitrogen flag (lookup_table.nc) on:
  beta_land_coubgc at 2xCO2 (esm-deck, piControl baseline), upstream sheet beta_land, and the
  historical cumulative net land flux 1960-2014 (H_land, tcre_evidence.historical) and ocean sink.
If models without N limitation compensate through tuning, their 1pctCO2 beta_L is higher while their
historical land sink is not. Small samples: reported as group means, Welch t and a permutation p-value.
Writes results/m7_nitrogen_historical.csv, results/m7_nitrogen_historical.json
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from tcre_evidence import feedbacks as fb
from tcre_evidence import historical as hist
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS


def perm_p(a, b, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.r_[a, b]
    d0 = abs(a.mean() - b.mean())
    k = len(a)
    c = 0
    for _ in range(n):
        p = rng.permutation(x)
        c += abs(p[:k].mean() - p[k:].mean()) >= d0
    return (c + 1) / (n + 1)


def main():
    lt = ud.lookup_table().sel(process="nitrogen").to_series()
    t = fb.table("2x", "picontrol")
    H = hist.table()
    sh = ud.cmip_sheet_dataset().drop_dims("process").to_dataframe()
    D = pd.DataFrame(dict(nitrogen=lt)).join(t[["beta_land_coubgc", "gamma_land_coubgc", "R"]], how="left")
    D = D.join(H[["H_land", "H_ocean"]], how="left").join(sh[["beta_land"]].rename(columns={"beta_land": "beta_land_sheet"}), how="left")
    D.to_csv(RESULTS / "m7_nitrogen_historical.csv", float_format="%.4g")
    out = {}
    for v in ["beta_land_coubgc", "beta_land_sheet", "H_land", "H_ocean", "gamma_land_coubgc", "R"]:
        d = D.dropna(subset=[v])
        a, b = d[d.nitrogen == 1][v].values, d[d.nitrogen == 0][v].values
        out[v] = dict(mean_with_N=float(a.mean()), mean_without_N=float(b.mean()), n_with=len(a), n_without=len(b),
                      ratio_with_over_without=float(a.mean() / b.mean()) if b.mean() else None,
                      welch_p=float(stats.ttest_ind(a, b, equal_var=False).pvalue), perm_p=float(perm_p(a, b)))
    both = D.dropna(subset=["beta_land_coubgc", "H_land"])
    out["corr_beta_land_H_land"] = dict(r=float(np.corrcoef(both.beta_land_coubgc, both.H_land)[0, 1]), n=len(both))
    json.dump(out, open(RESULTS / "m7_nitrogen_historical.json", "w"), indent=1)
    print(D.round(3).to_string())
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
