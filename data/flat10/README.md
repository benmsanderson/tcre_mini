# flat10 temperatures and metrics

- `temperatures_flat10.csv`: annual global-mean temperature anomaly (K) for the flat10 experiment
  (10 PgC/yr for 100 years), one column per model.
- `esm_metrics_flat10.csv`: derived metrics per model (TCRE, ZEC50, ZEC100 and others).

Copied unchanged from `output/temperatures_flat10.csv` and `esm_metrics_flat10.csv` in the esm-deck
repository (`benmsanderson/esm-deck`), commit eef04dc63ee128d3fa0b84ceabd101d8849c0ea4. They are read by
`tcre_evidence.historical.flat10_tcre()` (used in `scripts/m5_tcr_r.py`).
