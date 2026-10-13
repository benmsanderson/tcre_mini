"""M6 / A5: decompose the difference between the upstream posterior TCRE and the TCRE implied by the CMIP6
input means into contributions of each component to the denominator k + beta + alpha*gamma.

Evaluated at TCR = 1.7 K (alpha = 1.7/280), using component posterior means and CMIP input means from
results/m6_prior_info.csv (written by m6_prior_sensitivity.py). First-order, exact for the denominator.
Also the first-order TCRE effect of the gamma_land change from the nu evidence (results/m6_nu_evidence.json).
Writes results/m6_shift_decomposition.json
"""
import json

import pandas as pd

from tcre_evidence.paths import RESULTS
from tcre_evidence.upstream_data import K


def main():
    info = pd.read_csv(RESULTS / "m6_prior_info.csv").set_index("quantity")
    a = 1.7 / 280.0
    cm, po = info.cmip_mean, info.post_mean
    terms = {"beta_land": 1.0, "beta_ocean": 1.0, "gamma_land": a, "gamma_ocean": a}
    D_c = K + sum(cm[q] * w for q, w in terms.items())
    D_p = K + sum(po[q] * w for q, w in terms.items())
    d = {q: float((po[q] - cm[q]) * w) for q, w in terms.items()}
    tot = D_p - D_c
    out = dict(alpha=a, tcre_cmip_means=1e3 * a / D_c, tcre_posterior_means=1e3 * a / D_p,
               denominator_change=tot, contribution=d, share={q: v / tot for q, v in d.items()})
    # first-order effect of adding the paired-experiment nu evidence (A6) on TCRE, via gamma_land
    nu = json.load(open(RESULTS / "m6_nu_evidence.json"))["real"]
    g0 = [r for r in nu if not r["with_nu"]][0]["gamma_land_mean"]
    g1 = [r for r in nu if r["with_nu"]][0]["gamma_land_mean"]
    D_nu = D_p + a * (g1 - g0)
    out["nu_evidence"] = dict(gamma_land_without=g0, gamma_land_with=g1,
                              tcre_change_first_order=1e3 * a / D_nu - 1e3 * a / D_p)
    json.dump(out, open(RESULTS / "m6_shift_decomposition.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
