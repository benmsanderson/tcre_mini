"""M1: reproduce HANDOVER §2.3 (covariance and recombination) and check upstream's hard-coded alpha alignment.

Writes results/m1_covariance.json and results/m1_alpha_alignment.csv,
results/m1_input_sources.csv.
"""
import json

import numpy as np
import pandas as pd

from tcre_evidence import covariance as cv
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS

N, SEED = 400_000, 0


def table(values, label):
    rows = {}
    for s in cv.SCHEMES:
        rows[s] = cv.summary(cv.resample(values, s, N, SEED))
    print(f"\n{label}")
    print(f"  {'scheme':12s} {'sd AF':>6s} {'mean':>5s} {'sd':>6s} {'5-95 %':>10s}")
    for s, r in rows.items():
        print(f"  {s:12s} {r['sd_af']:6.3f} {r['tcre_mean']:5.2f} {r['tcre_sd']:6.3f} {r['tcre_p05']:.2f}-{r['tcre_p95']:.2f}")
    return rows


def main():
    cp = ud.carbon_parameters("2xCO2")
    models = [str(m) for m in cp.model.values]
    carbon = {c: cp[c].values for c in cv.COMPONENTS[1:]}

    # (a) as upstream notebook: hard-coded alpha assigned in intersect1d order
    v_hard = {"alpha": ud.HARDCODED_ALPHA, **carbon}
    corr_hard = cv.pair_correlations(v_hard)
    print("models (intersect1d order):", models)
    print("correlations, hard-coded alpha:", {k: round(x, 2) for k, x in corr_hard.items()})
    tab_hard = table(v_hard, "Hard-coded alpha in upstream order (HANDOVER §2.3 targets)")

    # (b) alpha matched by model name from the sheet 'alpha' column
    sheet = ud.cmip_sheet_dataset()
    a_sheet = sheet.alpha.reindex(model=cp.model).values
    align = pd.DataFrame(dict(model=models, alpha_hardcoded_as_assigned=ud.HARDCODED_ALPHA,
                              alpha_sheet_by_name=a_sheet))
    # nearest sheet model for each hard-coded value (identifies the order actually used)
    sheet_alpha = sheet.alpha.dropna("model")
    align["nearest_sheet_model_to_hardcoded"] = [
        str(sheet_alpha.model.values[np.argmin(np.abs(sheet_alpha.values - a))]) for a in ud.HARDCODED_ALPHA]
    align["nearest_sheet_alpha"] = [float(sheet_alpha.values[np.argmin(np.abs(sheet_alpha.values - a))])
                                    for a in ud.HARDCODED_ALPHA]
    print("\nalpha alignment:\n", align.to_string(index=False))
    ok = ~np.isnan(a_sheet)
    v_name = {"alpha": a_sheet[ok], **{c: x[ok] for c, x in carbon.items()}}
    corr_name = cv.pair_correlations(v_name)
    print(f"\ncorrelations, alpha matched by name ({ok.sum()} models):", {k: round(x, 2) for k, x in corr_name.items()})
    tab_name = table(v_name, f"Alpha matched by model name from sheet ({ok.sum()} models)")

    # (c) sources: sheet values (used for posteriors) vs CSV + Terhaar (used for CMIP reference)
    src = pd.DataFrame(dict(model=models))
    for c in carbon:
        src[f"{c}_csv_terhaar"] = carbon[c]
        src[f"{c}_sheet"] = sheet[c].reindex(model=cp.model).values
    src.to_csv(RESULTS / "m1_input_sources.csv", index=False, float_format="%.4f")
    print("\nupstream inputs, CSV/Terhaar (reference) vs sheet (posteriors):\n", src.round(3).to_string(index=False))

    RESULTS.mkdir(exist_ok=True)
    align.to_csv(RESULTS / "m1_alpha_alignment.csv", index=False)
    json.dump(dict(n=N, seed=SEED, models_intersect_order=models,
                   hardcoded_alpha=dict(correlations=corr_hard, schemes=tab_hard),
                   alpha_by_name=dict(models=[m for m, o in zip(models, ok) if o], correlations=corr_name,
                                      schemes=tab_name)),
              open(RESULTS / "m1_covariance.json", "w"), indent=2)


if __name__ == "__main__":
    main()
