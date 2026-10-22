"""WIEMIP values tables and paired process effects (2xCO2 window, years 61-80).

Reads data/wiemip/ (written by scripts/wiemip_extract.py) through tcre_evidence.wiemip. No sampling, so
no seed is needed.

Writes results/wiemip/:
- feedbacks_long.csv       every (model, factorial, pattern) with beta, gamma (coubgc, rad) and eps
- values_ukesm.csv         baseline, UKESM pattern (all models), in the toolkit's values format
- values_patternmean.csv   baseline, mean over the patterns each model ran
- paired_effects.csv       one-process factorial pairs: eta for beta, delta for gamma
- screening.csv            runs dropped (window not covered or non-finite) and ctrl drift
"""
import warnings

import pandas as pd

from tcre_evidence import wiemip as wm
from tcre_evidence.paths import RESULTS

OUT = RESULTS / "wiemip"
LEVEL = "2x"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cl, fc = wm.load()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fb = wm.feedbacks(cl, fc, LEVEL)
    for w in caught:
        print(w.message)
    fb.to_csv(OUT / "feedbacks_long.csv", index=False, float_format="%.6g")
    wm.values_table(fb, "ukesm").to_csv(OUT / "values_ukesm.csv", index=False, float_format="%.6g")
    wm.values_table(fb, "mean").to_csv(OUT / "values_patternmean.csv", index=False, float_format="%.6g")
    pe = wm.paired_effects(fb)
    pe.to_csv(OUT / "paired_effects.csv", index=False, float_format="%.6g")
    _, dropped = wm.screen(cl, LEVEL)
    drift = wm.ctrl_drift(cl).rename("ctrl_drift_pgc_per_century").reset_index()
    scr = drift.merge(dropped.assign(sim_forcing=dropped.sim + "/" + dropped.forcing)
                      [["model", "factorial", "sim_forcing", "reason"]], how="outer", on=["model", "factorial"])
    scr.to_csv(OUT / "screening.csv", index=False, float_format="%.4g")

    pd.set_option("display.width", 160)
    dT, dCa = fb.dT.iloc[0], fb.dCa.iloc[0]
    print(f"\n{LEVEL}: dT = {dT:.3f} K (FaIR median), dCa = {dCa:.1f} ppm; {len(fb)} rows, "
          f"{fb.model.nunique()} models")
    print("\nbaseline, UKESM pattern:")
    print(wm.values_table(fb, "ukesm").round(2).to_string(index=False))
    print("\npaired effects (UKESM pattern):")
    print(pe[pe.pattern == "ukesm"].drop(columns=["pattern", "factorial_on", "factorial_off"])
          .round(2).to_string(index=False))
    print("\nctrl drift beyond 20 PgC/century:")
    print(drift[drift.ctrl_drift_pgc_per_century.abs() > 20].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
