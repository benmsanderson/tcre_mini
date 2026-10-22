# C4MIP global means (11 CMIP6 models)

`c4mip_cmip6_ascii_files_11_models/<model>/<model>_{1pc_BGC,1pc_COU,1pc_RAD,PI_CONTROL}.txt`: annual
global means from the CMIP6 C4MIP experiments 1pctCO2-bgc, 1pctCO2, 1pctCO2-rad and piControl (CO2,
temperature, land and ocean carbon fluxes, and for some models GPP, NPP, heterotrophic respiration and
carbon pools).

- **Source.** Copied unchanged from `C4MIP_data/c4mip_cmip6_ascii_files_11_models/` in the esm-deck
  repository (`benmsanderson/esm-deck`), commit eef04dc63ee128d3fa0b84ceabd101d8849c0ea4. Only these
  files were used from that repository, so it is no longer a submodule.
- **Formats differ between models.** Column names, flux sign conventions, temperature units and year
  stamps all vary. Do not read the files directly: use `tcre_evidence.c4mip.read_run()`, which
  harmonises them with an explicit per-model table and asserts the flux signs.
- **Original data.** CMIP6 output from the modelling groups, distributed under the CMIP6 terms of use.
