"""M2 / A1: per-model feedbacks from esm-deck C4MIP files, closure test, comparison with upstream inputs.

Writes:
  data/derived/c4mip_feedbacks_{2x,4x}_{picontrol,forced_start}.csv
  results/m2_closure.csv, results/m2_upstream_comparison.csv, results/m2_summary.json
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import json

import numpy as np
import pandas as pd

from tcre_evidence import feedbacks as fb
from tcre_evidence import identity as idt
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import DERIVED, RESULTS

CONV = ["coubgc", "bgcrad", "courad", "bgconly_rad"]


def main():
    RESULTS.mkdir(exist_ok=True)
    summary, closure = {}, []
    tabs = {}
    for level in ["2x", "4x"]:
        for base in ["picontrol", "forced_start"]:
            t = fb.table(level, base)
            tabs[(level, base)] = t
            t.to_csv(DERIVED / f"c4mip_feedbacks_{level}_{base}.csv", float_format="%.5g")
            s = dict(n=len(t), tcre_mean=t.tcre.mean(), tcre_sd=t.tcre.std(), af_mean=t.af.mean(),
                     tcr_mean=t.tcr.mean(), R_mean=t.R.mean(), R_cv=t.R.std() / t.R.mean(),
                     tcr_cv=t.tcr.std() / t.tcr.mean())
            for c in CONV:
                r = idt.reconstruct(t, c)
                s[f"{c}_af_err_mean"], s[f"{c}_af_err_sd"] = r.af_err.mean(), r.af_err.std()
                s[f"{c}_tcre_err_mean"], s[f"{c}_tcre_err_sd"] = r.tcre_err.mean(), r.tcre_err.std()
                closure.append(r.assign(level=level, baseline=base, convention=c).reset_index())
            for v in ["beta_land_coubgc", "beta_ocean_coubgc", "gamma_land_coubgc", "gamma_ocean_coubgc",
                      "gamma_land_rad", "gamma_ocean_rad", "eps_land", "eps_ocean"]:
                s[f"{v}_mean"], s[f"{v}_sd"] = t[v].mean(), t[v].std()
            summary[f"{level}_{base}"] = {k: float(v) for k, v in s.items()}
            print(f"\n== {level} {base}: TCRE {s['tcre_mean']:.2f} ± {s['tcre_sd']:.2f}; "
                  f"CV TCR {s['tcr_cv']:.2f}, CV R {s['R_cv']:.2f}")
            for c in CONV:
                print(f"   {c:12s} AF err {s[c+'_af_err_mean']:+.4f} ± {s[c+'_af_err_sd']:.4f}   "
                      f"TCRE err {s[c+'_tcre_err_mean']:+.3f} ± {s[c+'_tcre_err_sd']:.3f}")
    pd.concat(closure).to_csv(RESULTS / "m2_closure.csv", index=False, float_format="%.5g")

    # A1 comparison: our COU-BGC values vs upstream CSV (land), Terhaar netCDF (ocean), sheet
    t2 = tabs[("2x", "picontrol")]
    t4 = tabs[("4x", "picontrol")]
    land2, land4 = ud.land_csv("2xCO2"), ud.land_csv("4xCO2")
    oc = ud.terhaar_ocean(69)
    sheet = ud.cmip_sheet_dataset()
    cmp = pd.DataFrame(index=t2.index)
    cmp["beta_land_ours"] = t2.beta_land_coubgc
    cmp["beta_land_upcsv"] = land2.beta_land.reindex(model=t2.index).values
    cmp["beta_land_sheet"] = sheet.beta_land.reindex(model=t2.index).values
    cmp["gamma_land_ours"] = t2.gamma_land_coubgc
    cmp["gamma_land_ours_rad"] = t2.gamma_land_rad
    cmp["gamma_land_upcsv"] = land2.gamma_land.reindex(model=t2.index).values
    cmp["gamma_land_sheet"] = sheet.gamma_land.reindex(model=t2.index).values
    cmp["beta_ocean_ours"] = t2.beta_ocean_coubgc
    cmp["beta_ocean_terhaar"] = oc.beta_ocean.reindex(model=t2.index).values
    cmp["gamma_ocean_ours"] = t2.gamma_ocean_coubgc
    cmp["gamma_ocean_ours_rad"] = t2.gamma_ocean_rad
    cmp["gamma_ocean_terhaar"] = oc.gamma_ocean.reindex(model=t2.index).values
    cmp["gamma_ocean_sheet"] = sheet.gamma_ocean.reindex(model=t2.index).values
    cmp["alpha_ours"] = t2.alpha
    cmp["alpha_sheet"] = sheet.alpha.reindex(model=t2.index).values
    cmp["gamma_land4x_ours"] = t4.gamma_land_coubgc
    cmp["gamma_land4x_upcsv"] = land4.gamma_land.reindex(model=t4.index).values
    cmp.to_csv(RESULTS / "m2_upstream_comparison.csv", float_format="%.4g")
    pd.set_option("display.width", 250)
    print("\nA1 comparison (2x, piControl baseline):\n", cmp.round(3).to_string())
    corr = {c: float(cmp[f"{c}_ours"].corr(cmp[o])) for c, o in
            [("beta_land", "beta_land_upcsv"), ("gamma_land", "gamma_land_upcsv"),
             ("beta_ocean", "beta_ocean_terhaar"), ("gamma_ocean", "gamma_ocean_terhaar")]}
    rad_corr = float(cmp.gamma_ocean_ours_rad.corr(cmp.gamma_ocean_terhaar))
    summary["upstream_comparison"] = dict(corr_ours_vs_upstream=corr, corr_gamma_ocean_rad_vs_terhaar=rad_corr,
                                          mean_diff={c: float((cmp[c + "_ours"] - cmp[o]).mean()) for c, o in
                                                     [("beta_land", "beta_land_upcsv"), ("gamma_land", "gamma_land_upcsv"),
                                                      ("beta_ocean", "beta_ocean_terhaar"), ("gamma_ocean", "gamma_ocean_terhaar")]})
    print("\ncorrelations ours vs upstream:", {k: round(v, 3) for k, v in corr.items()}, " gamma_O_rad vs Terhaar:", round(rad_corr, 3))
    # Closure of upstream inputs: TCRE reconstructed from the sheet (posterior inputs) and from
    # CSV/Terhaar + hard-coded alpha (CMIP reference) against actual 1pctCO2 TCRE from the deck files.
    sh = sheet.drop_dims("process").to_dataframe()
    up = pd.DataFrame(index=t2.index)
    up["tcre_actual"] = t2.tcre
    up["tcre_rec_sheet"] = idt.tcre(sh.alpha, sh.beta_land + sh.beta_ocean, sh.gamma_land + sh.gamma_ocean).reindex(t2.index)
    up["tcre_rec_sheet_ouralpha"] = idt.tcre(t2.alpha, (sh.beta_land + sh.beta_ocean).reindex(t2.index),
                                              (sh.gamma_land + sh.gamma_ocean).reindex(t2.index))
    cp = ud.carbon_parameters("2xCO2").to_dataframe()
    cp["alpha_hard"] = ud.HARDCODED_ALPHA
    up["tcre_rec_reference_hardalpha"] = idt.tcre(cp.alpha_hard, cp.beta_land + cp.beta_ocean,
                                                  cp.gamma_land + cp.gamma_ocean).reindex(t2.index)
    up["err_sheet"] = up.tcre_rec_sheet - up.tcre_actual
    up["err_sheet_ouralpha"] = up.tcre_rec_sheet_ouralpha - up.tcre_actual
    up["err_reference_hardalpha"] = up.tcre_rec_reference_hardalpha - up.tcre_actual
    up.to_csv(RESULTS / "m2_upstream_closure.csv", float_format="%.4g")
    print("\nUpstream-input closure vs actual 1pctCO2 TCRE (K/EgC):\n", up.round(3).to_string())
    summary["upstream_closure"] = {c: dict(mean=float(up[c].mean()), sd=float(up[c].std()),
                                           max_abs=float(up[c].abs().max()), n=int(up[c].notna().sum()))
                                   for c in ["err_sheet", "err_sheet_ouralpha", "err_reference_hardalpha"]}
    print(summary["upstream_closure"])
    json.dump(summary, open(RESULTS / "m2_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
