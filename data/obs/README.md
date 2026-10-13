# Observational data

Raw files are in `raw/` exactly as downloaded; tidy CSVs in this directory are
built from them by `python -m tcre_evidence.obs` (which also re-downloads).
Loaders: `tcre_evidence.obs.gcb()`, `noaa_co2()`, `hadcrut5()`.
Retrieval date for all files: 2026-10-13.

## Global Carbon Budget 2025 (v1.0)

- Raw: `raw/Global_Carbon_Budget_2025_v1.0.xlsx`
  (sha256 `a928cf06c57576b66761d1fec8224c9a41a781eb63442afc03ba902858dd64a9`).
- Source: ICOS Carbon Portal object `qSjPBsV1drZnYdH-yCJMmkGn`,
  https://meta.icos-cp.eu/objects/qSjPBsV1drZnYdH-yCJMmkGn, data citation
  "Global Carbon Project, 2025. Supplemental data of Global Carbon Budget 2025",
  https://doi.org/10.18160/GCP-2025 (linked from https://globalcarbonbudget.org/gcb-2025/).
  Paper: Friedlingstein et al. 2025 (see the workbook's Summary tab).
- Tidy: `gcb_global.csv` (sheet "Global Carbon Budget", 1959-2024) and
  `gcb_historical.csv` (sheet "Historical Budget", 1750-2024).
- Units: GtC/yr = PgC/yr. Columns: `fossil` (fossil emissions excluding
  carbonation), `luc` (land-use change emissions), `atm_growth`,
  `ocean_sink`, `land_sink`, `cement_carbonation_sink`, `budget_imbalance`.
  Sinks are positive for uptake. Uncertainties (1 sigma) are given as text in
  the workbook header rows and are not tabulated.

## NOAA GML CO2

- Raw: `raw/co2_annmean_gl.txt`, `raw/co2_gr_gl.txt`, `raw/co2_annmean_mlo.txt`
  from https://gml.noaa.gov/webdata/ccgg/trends/co2/ (file creation stamps
  9 Oct 2026; data through 2025).
- Tidy: `noaa_co2.csv`, columns `co2_gl` (global marine boundary layer annual
  mean, ppm, 1979-2025), `growth_gl` (global annual growth rate, ppm/yr,
  1959-2025; pre-1980 values approximated from Mauna Loa and South Pole),
  `co2_mlo` (Mauna Loa annual mean, ppm, 1959-2025), each with `*_unc`
  (1 sigma). The last year is preliminary per NOAA.
- Conversion used elsewhere: 1 ppm = 2.124 PgC.

## HadCRUT5

- Version HadCRUT.5.2.0.0 (updated by the Met Office on 18/09/2026), analysis
  (infilled) ensemble, summary series, global annual.
- Raw: `raw/HadCRUT.5.2.0.0.analysis.summary_series.global.annual.csv` from
  https://www.metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.2.0.0/analysis/diagnostics/HadCRUT.5.2.0.0.analysis.summary_series.global.annual.csv
- Tidy: `hadcrut5_global.csv`, columns `anomaly` (K, relative to 1961-1990),
  `lower_2p5`, `upper_97p5` (95% range). 1850-2026; the 2026 row is a partial
  year, dropped by default in `hadcrut5()`.
