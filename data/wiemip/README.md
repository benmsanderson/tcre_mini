# WIEMIP 1pctCO2 global land carbon

Written by `scripts/wiemip_extract.py` on the WIEMIP JupyterHub (data root `/srv/wiemip`) through the
hub's `wiemip_registry`; settings and registry code hash in `provenance.json`.

- `cland_annual.csv`: annual-mean global land carbon (PgC, `land_carbon_stock`) for every model, factorial
  and 1pctCO2 run the registry reaches (ctrl and bgc with forcing `stable`; cou and rad per GCM pattern).
- `forcing.csv`: CO2 (ppm) from the WIEMIP 1pctCO2 input file and the FaIR median GMST anomaly (K, config
  1051936) from which all three patterns' climate drivers were built.

Unpublished WIEMIP output: the repository is private, so share it only with WIEMIP authors.
