"""WIEMIP detection power before the flag table exists: cross-model flag regression against paired factorials.

The full flag table (which models include each process) is being compiled separately, so the real flag
budget cannot be run yet. Detection power depends on the flag table only through how many models are
flagged and how the flags overlap, so it is computed here on hypothetical designs:
- the baseline models (UKESM pattern), one flag set for k < n of them, plus 3 further flags set at random;
- residual scatter = the WIEMIP ensemble sd times sqrt(1 - r2), for r2 = 0 (flags explain nothing) and
  0.2 (about the share the four flags explained for gamma_L in CMIP6); also with a robust sd (1.4826 MAD),
  which limits the weight of single outlying models.
The effect sizes are the paired factorial effects (process on minus off in the same model; baseline
background only, UKESM pattern), whose mean has its own standard error sd/sqrt(n_pairs).
beta is treated on the log scale (multiplicative effects, eta = beta_on / beta_off).

Writes results/wiemip/flags_from_factorials.csv (data/wiemip is input only), results/wiemip/power_paired.csv
and results/wiemip/power_designs.csv.
"""
import numpy as np
import pandas as pd

from tcre_evidence import process_tests as pt
from tcre_evidence import wiemip as wm
from tcre_evidence.paths import RESULTS

SEED = 0
OUT = RESULTS / "wiemip"
N_OTHER_FLAGS = 3
KS = [3, 6, 9, 12, 15]
R2S = [0.0, 0.2]
QUANTS = {"gamma_land_coubgc": ("delta_gamma_land_coubgc", False), "beta_land_coubgc": ("eta_beta_land", True)}


def n_needed(s: float, eff: float, max_models: int = 4000) -> float:
    """Smallest balanced ensemble (each flag p = 0.5) with standard error <= |eff| / 2."""
    n_max = max_models - 2 if max_models % 2 == 0 else max_models - 1  # last n the scan below reaches
    if s * pt.design_se_factor(n_max, n_max // 2, N_OTHER_FLAGS, seed=SEED, ndraw=50) > abs(eff) / 2:
        return np.inf  # not reached within max_models (se falls with n), so skip the scan
    for n in range(8, max_models, 2):
        if s * pt.design_se_factor(n, n // 2, N_OTHER_FLAGS, seed=SEED, ndraw=50) <= abs(eff) / 2:
            return n
    return np.inf


def main():
    cl, fc = wm.load()
    fb = wm.feedbacks(cl, fc, "2x")
    wm.factorial_flags(cl).to_csv(OUT / "flags_from_factorials.csv")
    vals = wm.values_table(fb, "ukesm").set_index("model")
    pe = wm.paired_effects(fb).query("pattern == 'ukesm' and background == ''")
    n = len(vals)

    paired, designs = [], []
    for q, (ecol, log) in QUANTS.items():
        y = np.log(vals[q]) if log else vals[q]
        sds = {"sd": y.std(ddof=1), "robust_sd": 1.4826 * (y - y.median()).abs().median()}
        for proc, g in pe.groupby("process"):
            e = np.log(g[ecol]) if log else g[ecol]
            eff, e_sd = e.mean(), e.std(ddof=1)
            row = dict(quantity=q, scale="log" if log else "linear", process=proc, n_pairs=len(e),
                       mean_effect=eff, sd_between_models=e_sd, se_paired=e_sd / np.sqrt(len(e)),
                       models=" ".join(g.model))
            for sname, sd in sds.items():
                for r2 in R2S:
                    s = sd * np.sqrt(1 - r2)
                    row[f"n_needed_{sname}_r2_{r2}"] = n_needed(s, eff) if len(e) > 1 else np.nan
            paired.append(row)
        for sname, sd in sds.items():
            for r2 in R2S:
                for k in (k for k in KS if k < n):  # a flag on every model is the intercept
                    designs.append(dict(quantity=q, scale="log" if log else "linear", scatter=sname, r2=r2,
                                        n_models=n, n_flagged=k, resid_sd=sd * np.sqrt(1 - r2),
                                        se_cross_model=sd * np.sqrt(1 - r2) *
                                        pt.design_se_factor(n, k, N_OTHER_FLAGS, seed=SEED)))
    paired, designs = pd.DataFrame(paired), pd.DataFrame(designs)
    paired.to_csv(OUT / "power_paired.csv", index=False, float_format="%.4g")
    designs.to_csv(OUT / "power_designs.csv", index=False, float_format="%.4g")

    pd.set_option("display.width", 200)
    print(f"{n} baseline models (UKESM pattern); hypothetical designs with {N_OTHER_FLAGS} further flags\n")
    print(paired.drop(columns="models").round(3).to_string(index=False))
    print()
    print(designs.pivot_table(index=["quantity", "scatter", "r2"], columns="n_flagged",
                              values="se_cross_model").round(3).to_string())


if __name__ == "__main__":
    main()
