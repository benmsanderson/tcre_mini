# Upstream data snapshot

Retrieved 2026-10-13.

- `gsheet_*.csv`: Google Sheet `17_GeL1M1g0ilXK8w2j8803ewMeJxapWcf3mCFx4LGBE` (title "CMIP6"), tabs `beta_gamma`, `gamma_iav`, `Ocean predictors`, read through the Google Sheets API (the gviz CSV endpoint used by the upstream notebook is not reachable from the build environment). Blank spacer rows in `Ocean predictors` and `gamma_iav` are dropped; values are otherwise verbatim strings.
- `DATA/`: copy of `external/TCRE_Assessment/DATA` at commit 98c74cb8d4ccea8327cf8d325c9fe6266a2ff122.

Notes recorded at retrieval:
- `gamma_iav` sheet now includes UKESM1-0-LL (-6.79 ± 0.72), which is absent from the notebook's saved output; the OBS row is -5.20 ± 0.74 (HANDOVER §3.1 quotes -4.3 ± 0.67, the value hard-coded in `data_utils.load_emergent_constraint_evidence`).
- The sheet is live; the upstream headline (1.766 [0.927, 2.737]) was produced from whatever the sheet held at the time of the notebook run.
