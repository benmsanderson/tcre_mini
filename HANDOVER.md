# HANDOVER: TCRE lines-of-evidence testbed

Handover for Claude Code. Read this file in full before writing any code. It records the scientific context, what has already been established, the two workstreams, the data needed, and the order of work.

Owner: Ben Sanderson (CICERO). Context: IPCC AR7 WGI Chapter 5 (TCRE, remaining carbon budgets, ZEC).

Writing conventions for all docs, comments and figure text in this repo: UK English; no em-dashes (use commas, semicolons, colons or parentheses); neutral, technical register. This work concerns a community method developed by colleagues, so describe that method's behaviour factually ("the posterior is prior-dominated for δ_permafrost") and never characterise people or intentions.

---

## 1. Purpose

Two linked aims:

- **Workstream A, test the community framework.** Stress-test the Bayesian process-based TCRE framework in `netzeroasap/TCRE_Assessment`, which is likely to be a primary input to the WCRP TCRE assessment and AR7. Establish, with reproducible diagnostics, whether its β/γ decomposition, process-factorial adjustments and emergent constraints add information about TCRE beyond the raw CMIP ensemble.
- **Workstream B, build a convention-free alternative line.** Implement a simpler estimate, TCRE = TCR × R, where R is the per-model carbon factor (TCRE/TCR) diagnosed directly. Constrain R with target-matched historical observables (airborne fraction and/or sinks), and validate the historical-to-doubling transfer in a perfect-model setting. This parallels unpublished work by Séférian, Ribes and Qasmi (EGU25-12981: AF at 2×CO2 41–59 %, TCRE 1–2.1 K/EgC, very likely); a collaboration has been proposed. Build our implementation independently so it can serve as validation or contribution either way.

The deliverable is a short paper's worth of results: diagnostics, figures, and a clean, tested codebase.

---

## 2. Scientific background (what we established before this repo existed)

### 2.1 The identity

Jones & Friedlingstein (2020, ERL 15 074019; "JF2020") write, for concentration-driven 1pctCO2 at time t (normally year 70, 2×CO2):

```
E      = k·ΔCa + ΔC_L + ΔC_O                      (compatible cumulative emissions, PgC)
ΔC_x   = β_x·ΔCa + γ_x·ΔT                         (x = land, ocean)
α      = ΔT / ΔCa                                  (K ppm⁻¹; ≈ TCR / 280 ppm at 2×CO2)
AF     = k / (k + β + α·γ)                         (dimensionless airborne fraction)
TCRE   = α / (k + β + α·γ)                         (K PgC⁻¹)
k      = 2.12 PgC ppm⁻¹
```

The C4MIP convention (Arora et al. 2013, 2020) takes β from the biogeochemically coupled run (1pctCO2-bgc) and γ from the difference between the fully coupled run (COU, 1pctCO2) and BGC, solving a 2×2 system that accounts for the small warming in BGC:

```
ΔC_COU = β·ΔCa + γ·ΔT_COU
ΔC_BGC = β·ΔCa + γ·ΔT_BGC
```

**Key consequence: under this convention the reconstruction of AF and TCRE is exact by construction for each model** at the diagnosis time. It is bookkeeping, not a test of linearity. A good fit of reconstructed against actual TCRE across models (JF2020 Fig. 2; CMIP6 AF 0.532 actual vs 0.534 reconstructed, TCRE 1.78 vs 1.75) carries no evidential weight about the framework.

### 2.2 The convention problem

With the radiatively coupled run (1pctCO2-rad):

```
ΔC_COU = β·ΔCa + γ_RAD·ΔT + ε
γ_COU−BGC = γ_RAD + ε/ΔT
```

ε is the non-additive interaction term. It is mostly temperature degradation of the CO2-driven sink, ε ≈ (∂β/∂T)·ΔCa·ΔT. So the C4MIP γ contains part of the β response evaluated in a warmer state. Literature:

- Gregory et al. 2009 (J. Clim. 22, 5232): feedbacks do not combine linearly; RAD needed.
- Zickfeld et al. 2011 (J. Clim.): nonlinearity in UVic.
- Schwinger et al. 2014 (J. Clim. 27, 3869): ocean RAD+BGC exceeds COU uptake by 19–58 PgC in CMIP5, the same order as γ_O itself; cause is circulation and mixing weakening reducing transport of anthropogenic carbon to depth.
- Arora et al. 2013, 2020 (Biogeosciences 17, 4173): report both conventions; ocean γ is always larger in magnitude under COU−BGC.
- Anecdotally (P. Friedlingstein, pers. comm.), RAD-based γ "does not work" for reconstruction. One aim is to demonstrate this quantitatively.

Implications:
- β/γ attribution of spread is convention-dependent.
- γ is not a transferable system constant: ε scales with ΔCa·ΔT, so γ_COU−BGC drifts with forcing level. At fixed ΔCa, swapping in an external α is first-order consistent; transfer across forcing levels or pathways is not.
- State dependence is visible directly: in the upstream data, γ_L(4×) is more negative than γ_L(2×) in 10 of 13 CMIP6 models (CNRM −52 → −85, UKESM −27 → −37); CanESM5 changes sign.

### 2.3 Covariance

JF2020 eqs. 7–9 propagate dα, dβ, dγ as independent. Back-of-envelope from their tables suggests that overstates CMIP6 AF variance by about a factor of four (implied corr(β, αγ) ≈ −0.8; noisy, multi-generation). In the upstream CMIP6 data (11 common models, 2×CO2) the correlations are weaker but material:

| pair | corr |
|---|---|
| β_L, γ_L | −0.38 |
| β, α·γ | −0.29 |
| α, α·γ | −0.53 |
| β_O, γ_O | −0.56 |
| α, β_O | +0.52 |

Recombining independently sampled components (as the upstream notebook does for β_L, γ_L, α) inflates the carbon-term spread and **truncates the upper tail of TCRE**, because high-TCR models also lose more carbon per unit CO2:

| sampling | sd of AF term | TCRE mean | TCRE sd | TCRE 5–95 % (K/EgC) |
|---|---|---|---|---|
| model-paired (joint) | 0.036 | 1.84 | 0.475 | 1.31–2.85 |
| upstream pattern (β_L, γ_L, α independent; ocean paired) | 0.044 | 1.84 | 0.443 | 1.26–2.59 |
| fully independent | 0.044 | 1.84 | 0.442 | 1.26–2.59 |

Use these as regression-test targets (script logic reproduced in `scratch/cov_check.py` if present; otherwise rebuild from §5.1).

### 2.4 Preliminary result from the esm-deck C4MIP files (crude; to be redone properly in A1/A2)

Quick computation on the 11-model CMIP6 C4MIP global-mean files (§6.1), year 70 (20-year mean), temperature baseline = first 10 years of each run (not piControl), CO2 baseline = year 1, flux signs harmonised:

- COU−BGC reconstruction closes exactly (mean AF error 0.0000), as expected.
- RAD-based γ reconstruction: AF error mean +0.008, sd 0.025; TCRE error mean +0.03, sd 0.08 K/EgC. Most models are within ±0.05 K/EgC; the outlier is ACCESS-ESM1.5 (+0.23 K/EgC, land ε ≈ +116 PgC).
- ε_L mostly positive (coupled land uptake exceeds BGC+RAD), ε_O mostly negative (consistent with Schwinger et al. 2014), CNRM ocean the exception (+49 PgC). Land and ocean ε partly cancel in the total.

Implication: at 2×CO2 the "RAD does not work" statement is model-specific rather than general, while the land/ocean and β/γ *partition* is strongly convention-dependent. Check whether this holds with proper piControl baselines and at 4×CO2. Do not quote these numbers until A1/A2 reproduce them.

TCRE across these 11 models (year 70, COU): mean 1.70, sd 0.36 K/EgC.

### 2.5 Scale of the problem

TCR dominates TCRE variance. Upstream TCR prior sd/mean ≈ 21 %; CMIP6 carbon-factor spread ≈ 8 %. So roughly 85–90 % of TCRE variance comes from α. The carbon-cycle machinery mostly moves the central estimate and the upper tail, and drives the *narrative* about which processes matter. Keep this proportion in view when presenting results.

---

## 3. The community framework (upstream)

- Repo: https://github.com/netzeroasap/TCRE_Assessment, inspected at commit `98c74cb8d4ccea8327cf8d325c9fe6266a2ff122` (2026-08-12). **Pin this commit** (git submodule under `external/TCRE_Assessment`). Never edit upstream code; wrap or re-implement in our package.
- Companion methods paper: Marvel, Jones, Cooper, Steinert 2026, "Using climate models in climate assessments: a Bayesian perspective", Frontiers in Climate 8, doi:10.3389/fclim.2026.1839130 (general interpretive models; no TCRE results, no validation tests). Related package: github.com/netzeroasap/CMIP_BAYES.
- WCRP context: TCRE/ZEC assessment report (Vienna 2025) lists a Bayesian TCRE framework as high priority; interpretive model "not yet fully agreed".

### 3.1 Structure

`TCRE = α·AF_k`, with `AF_k = 1/(k + β + αγ)` (note: upstream calls this `AF`, but it is AF/k with units ppm/PgC; TCRE is still correct).

| Component | Upstream model | File |
|---|---|---|
| α | TCR ~ N(1.7, 0.36) K, α = TCR/280 | `config/priors.py`, notebook |
| β_L | multiplicative process model: β_L,m = base·Π η_p^{L_pm}, η_N pseudo-obs 0.69 (ACCESS CN/C), 0.77 (UKESM ctrl/nonlim), σ = 0.1 | `bayes/model_builders.py::build_beta_land_model` |
| γ_L | additive process model, split tropics/extratropics: γ = base + Σ(δ_p + ε_pm)·L_pm; permafrost forced ≈ 0 in tropics; δ_permafrost −LogNormal(log 30, 1); δ_N LogNormal(log 10, 1) | `build_gamma_land_model` |
| γ_L tropics EC | linear EC on IAV of CO2 growth rate (Zechlau et al. 2022 values), obs −4.3 ± 0.67 PgC yr⁻¹ K⁻¹ (or "OBS" row of sheet) | `add_emergent_constraint` |
| β_O, γ_O | joint multivariate EC: [log β_O, standardised γ_O] → (AMOC, SSS, CUC), LKJ covariance; obs AMOC 17.02 ± 0.35, SSS 34.04 ± 0.04, CUC 0.546 ± 0.003 (Terhaar et al.) | `build_joint_ocean_model` |
| Combination | land posteriors (separate PyMC models) + ocean joint posterior + independent α samples, added sample-wise | `notebooks/TCRE_total.ipynb` cells 31–35 |

Upstream headline (notebook output): TCRE posterior mean 1.766 K/EgC, 95 % HDI [0.927, 2.737].

Core process-model code: `bayes/CMIP_process.py::build_vectorized_process_model`. Model expectation `mu = X_unscaled·mult_eff + delta_eff + bias_common + bias_m`; reported truth `X_true = X_unscaled·Π η + Σ δ` (i.e. "a model with every process on at ensemble-average effect").

### 3.2 Data sources upstream

- `DATA/TCREsource_betagamma.csv`: land β_L, γ_L at 2× and 4×CO2 for 13 CMIP6/CMIP6+ models plus process-ensemble members (ACCESS-ESM C/CN/CNP, UKESM1-1 variants, IPSL_Perm_LandN, CNRM). Convention not documented; assume COU−BGC, verify.
- `DATA/beta_gamma_terhaar.nc`: ocean β, γ, annual, 140 years, 13 models; global attribute says both come from "1pctCO2-BGC", which cannot be right for γ. **Verify convention before use.**
- `DATA/Terhaar.txt`: AMOC, SSS, CUC per model.
- `DATA/lookup_table.nc`: process flags (nitrogen, permafrost, fire, veg) × 13 models.
- Notebook loads main inputs from a Google Sheet (key `17_GeL1M1g0ilXK8w2j8803ewMeJxapWcf3mCFx4LGBE`, sheets `beta_gamma`, `ocean predictors`, `gamma_iav`). The CMIP reference distribution uses local files instead. Snapshot the sheet into `data/upstream_snapshot/` with a retrieval date; do not depend on it live.

### 3.3 Known issues found in upstream (to test, not to "fix" upstream)

1. Land and ocean β/γ from different sources; convention and time-point consistency unverified (notebook itself notes "we need a spreadsheet with all β and γ calculated in the same way").
2. No closure test: reconstructed TCRE never compared with each model's actual 1pctCO2 TCRE.
3. β_L and γ_L fitted in separate models; α independent; sample-wise addition discards covariance (§2.3).
4. Process effects on β (η) and γ (δ) estimated separately with independent priors, although under COU−BGC part of a process's γ effect arises mechanically via β (the ε channel).
5. Process effects identified only by presence-flag regression on ~13 heterogeneous models; δ_permafrost prior-dominated (the notebook says so). Paired-experiment γ evidence (`data_utils.load_process_evidence()` computes ν_nitrogen, ν_fire, ν_vegetation) is never used. Paired experiments come from models already in the ensemble.
6. IAV observable (fast response of existing pools, ~present CO2) is closer to γ_RAD than to γ_COU−BGC at 2×CO2; relationship has shifted between CMIP generations (confirm against Zechlau et al. 2022).
7. Terhaar constraint developed for anthropogenic uptake (β-like); extension to γ_O is an extra assumption; 3×2 regression + LKJ on 11–13 models has few degrees of freedom; posterior tightness may come partly from regression priors.
8. Hard-coded α array (11 values) assumed to match `np.intersect1d` alphabetical order; verify.
9. No out-of-sample or perfect-model validation anywhere.
10. Double-counting risk at the assessment level: if permafrost and other missing feedbacks enter TCRE here, they must not also enter the remaining-carbon-budget "unrepresented Earth system feedback" adjustment.

---

## 4. Workstream A: diagnostics of the community framework

Priority order. A1–A3 are quick and should come first; A4 is the decisive test.

### A1. Data harmonisation and closure
- Recompute β, γ for land and ocean from raw CMIP6 output with one documented convention (Arora 2020 2×2 solve), same 20-year window centred on year 70 (and year 140 for 4×).
- Compare with upstream CSV and Terhaar netCDF; report discrepancies per model.
- Closure test: per model, TCRE_reconstructed (from α, β, γ) vs TCRE_actual (ΔT/E at year 70). Expected ≈ exact under consistent COU−BGC; any mismatch diagnoses inconsistency in upstream inputs.

### A2. Convention sensitivity
- For models with 1pctCO2-rad: compute γ_RAD and ε (land, ocean) at 2× and 4×.
- Reconstruct AF from (β_BGC, γ_RAD) and show the bias against actual AF (expect AF too low, ocean especially). This is the quantitative version of "RAD-based γ does not work".
- Redo a JF2020-style variance attribution under three allocations: COU−BGC (C4MIP), COU−RAD (puts ε in β), and symmetric factor separation (Stein & Alpert 1993 / Shapley split of ε). Report whether the "β dominates, γ minor" ranking survives.

### A3. Recombination and covariance
- Reproduce §2.3 numbers from upstream data.
- Rerun the upstream pipeline with (i) a joint bivariate land model for (β_L, γ_L) and (ii) α–carbon covariance retained (α conditioned on carbon terms across models). Report change in TCRE tails.
- Report both with the assessed TCR prior swapped in, and state explicitly which covariances survive that swap.

### A4. Perfect-model (leave-one-out) test: decisive
For each CMIP6 model m with the required observables:
1. Remove m from the training ensemble (and from process-ensemble evidence if m's family supplies paired experiments; also run a variant leaving out the whole model family).
2. Replace real-world observations with m's own simulated observables: IAV sensitivity (γ_IAV), AMOC, SSS, CUC (and TCR if testing α).
3. Refit the upstream framework unchanged (call upstream builders via the submodule).
4. Score the posterior for m's TCRE (and AF, β, γ components) against m's actual 1pctCO2 value.

Metrics: coverage of 66 % and 90 % intervals; CRPS; log score; skill relative to the raw ensemble (leave-one-out multi-model mean ± sd) and to an unweighted random-effects baseline. Pass criterion: calibrated coverage and CRPS skill > 0 against the raw ensemble.

Notes:
- With ~11–13 models, report per-model results, not just aggregates.
- Sampling cost: upstream notebook samples each component in ~10–65 s (4 chains, 1000 draws). The full LOO is roughly 13 × 3 component fits; feasible locally. Parallelise across models with joblib/multiprocessing; fix seeds.
- Process-flag component in LOO: the held-out model's flags are known; its "truth" is its own value, so this tests whether process regression predicts out of sample.

### A5. Prior versus data information
- Prior-predictive vs posterior for each component and for TCRE (`pm.sample_prior_predictive`).
- Sensitivity to the key priors: δ_permafrost, δ_nitrogen, η priors, σ_struct, regression priors, LKJ η.
- Quantify the share of posterior shift from priors (e.g. KL(prior‖posterior) per component; or rerun with flat-ish priors).

### A6. Process identifiability (synthetic)
- Synthetic truth experiment: generate pseudo-ensembles with known δ_p, η_p on the actual lookup table (13 models × 4 flags), plus realistic structural noise and confounders (e.g. nitrogen models also differ in base β). Refit and check recovery of δ_p, η_p. This shows directly whether the design can identify process effects.
- Test the ε channel: generate data where nitrogen scales β and γ contains κ·β·ΔCa/ΔT; show the separate η/δ model misattributes it.
- Incorporate the unused ν evidence and report its effect.

### A7. State dependence
- Rerun A1–A4 at 4×CO2 (upstream has `kind="4xCO2"` paths for land).

---

## 5. Workstream B: convention-free line, TCRE = TCR × R

### 5.1 Definitions
For each model m from 1pctCO2 (COU) at year 70 (20-year mean):
- TCR_m = ΔT_m
- E_m = compatible cumulative emissions = k·ΔCa + ΔC_L + ΔC_O, from integrated `nbp` and `fgco2` (1pctCO2 has no land-use change; check sign conventions)
- TCRE_m = TCR_m / E_m · 1000 (K/EgC)
- R_m = TCRE_m / TCR_m = 1/E_m (equivalently AF_m/(k·ΔCa))

Also compute from esm-flat10 (flat10MIP, CMIP7 fast track where available): TCRE at 1000 PgC, and the pathway difference TCRE_1pct − TCRE_flat10 per model (about 0.1 K/EgC lower in flat10 on average; Sanderson et al. 2025 GMD 18, 5699).

### 5.2 Estimator
- TCRE_real = TCR_assessed × R_real.
- R constrained by a historical observable H with cross-model relationship R_m ~ f(H_m), conditioned on observed H. Candidates for H:
  1. cumulative airborne fraction over 1959–2023 (NOAA growth rate + GCB fossil and LUC emissions);
  2. the time-varying AF–cumulative-emission relationship (the Séférian et al. approach);
  3. cumulative ocean + land sink (GCB), or ocean sink alone (well-constrained, β_O-like).
- Model-side historical quantities from `esm-hist` (emissions-driven, preferred) and from `historical` with diagnosed compatible emissions; document which models have which.
- Retain TCR–R covariance: fit (TCR_m, R_m, H_m) jointly (bivariate or trivariate hierarchical model); when substituting the assessed TCR, condition R on TCR via the cross-model relationship rather than sampling independently.
- Use PyMC for consistency with upstream, but keep the model small and transparent; compare with a simple Gaussian-conditioning (kriging, Ribes et al. 2021 style) closed form.

### 5.3 Validation
- **Historical-to-doubling transfer perfect-model test:** for each model m, use m's own historical H_m as pseudo-observation, train on the others, predict R_m (and TCRE_m given TCR_m). Same metrics as A4. This is the step reviewers will probe hardest (historical AF reflects near-exponential emissions growth; 1pctCO2 does not).
- Repeat with flat10 as the target pathway.
- Report the pathway-transfer term explicitly when comparing with historical TCRE estimates (e.g. DAMIP attribution ≈ 1.39 ± 0.09 K/EgC in upstream `TCRE_from_DAMIP.ipynb`).

### 5.4 Comparison
Final figure: TCRE distributions from (i) raw CMIP6, (ii) upstream framework as published, (iii) upstream with joint covariance, (iv) TCR × R unconstrained, (v) TCR × R historically constrained, (vi) historical attribution line, with perfect-model skill annotated for (ii), (iii), (v).

---

## 6. Data

### 6.0 Primary source: esm-deck (use first)

Repo: https://github.com/benmsanderson/esm-deck (Ben's; pin as a submodule at `external/esm-deck`, commit `eef04dc63ee128d3fa0b84ceabd101d8849c0ea4`). It already holds global-mean data for most of what A1, A2 and B need, so ESGF/Pangeo downloads become a fallback for gaps (historical/esm-hist runs, extra models).

Contents and caveats found on inspection:

- `C4MIP_data/c4mip_cmip6_ascii_files_11_models/<model>/<model>_1pc_{BGC,COU,RAD}.txt` and `_PI_CONTROL.txt`: annual global means, columns `YEAR ST CO2 LA_CO2_FLUX OA_CO2_FLUX GPP NPP HETRES CLIVE CDEAD`, for 11 CMIP6 models (ACCESS-ESM1.5, BCC, CESM2, CNRM-ESM2.1, CanESM5, GFDL, IPSL, JAMSTEC [= MIROC-ES2L], MPIESM, NorESM2LM, UKESM). This is the Arora et al. 2020 set, **including RAD runs**.
  - Formats differ: whitespace vs tab separated; some headers quoted (CNRM); MPIESM spells `HETERES`; UKESM uses `Year` and mid-year stamps; run lengths 140–152 rows. Normalise column names on read.
  - **Flux sign conventions differ by model.** Positive = uptake for CanESM5, CESM2, CNRM, GFDL, IPSL, JAMSTEC; negative = uptake for ACCESS, BCC, MPIESM, NorESM2LM, UKESM. Harmonise explicitly with a per-model table and assert that cumulative ocean uptake at year 70 is positive.
  - Fluxes are PgC/yr; `ST` in K; `CO2` in ppm (initial values 284–287 ppm, differ by model; use each model's own year-1 value as baseline, or piControl).
  - Use `_PI_CONTROL.txt` for baselines and drift correction rather than the first years of the forced run.
  - Combined tables `la_co2_flux_data_integer_years.csv`, `oa_co2_flux_data_integer_years.csv`, `co2_data_models.csv` (semicolon separated) carry the same mixed signs plus UVic; they are not harmonised.
  - Model-name mapping to upstream TCRE_Assessment names is needed (JAMSTEC ↔ MIROC-ES2L, GFDL ↔ GFDL-ESM4, IPSL ↔ IPSL-CM6A-LR, MPIESM ↔ MPI-ESM1-2-LR, NorESM2LM ↔ NorESM2-LM, UKESM ↔ UKESM1-0-LL, BCC ↔ BCC-CSM2-MR, ACCESS-ESM1.5 ↔ ACCESS-ESM1-5, CNRM-ESM2.1 ↔ CNRM-ESM2-1). Upstream has EC-Earth3-CC and CMCC-ESM2, which esm-deck lacks; esm-deck has JAMSTEC/MIROC and others.
- `output/{temperatures,land_carbon_sink,ocean_carbon_sink,atmospheric_co2}_flat10.csv`: harmonised annual global means for esm-flat10 for 10 models (ACCESS-ESM1-5, CESM2, GFDL-ESM4, GISS, NorESM2-LM, MPI-ESM1-2-LR, CNRM-ESM2-1, UKESM, UVic-ESCM-2-10, HadCM3LC-Bris), 318 rows. Land and ocean sinks appear to be **cumulative PgC**; check whether temperatures are anomalies and over what baseline, and what the 318 rows span (flat10 plus continuation). `atmospheric_co2_flat10.csv` has at least one column stored as bracketed strings (`[284.9072]`, MPI); parse defensively.
- `esm_metrics_flat10.csv`: per-model flat10 metrics (tcre, zec50, zec100, tpw, tnz, tr1000, tr0). Use as a regression target for our own flat10 TCRE computation.
- `ESM_data/<model>/{flat10,flat10-zec,flat10-cdr}_<var>_<model>.nc`: raw per-model files, heterogeneous (monthly vs annual, variable names such as `TAS_GLOB`, `temp_mm_1_5m`, `nep` stored in an `nbp` file; units K vs C). Prefer the harmonised `output/` CSVs; go to raw only if needed, and reuse the processing logic in `flat10_esm.ipynb`.
- `ZECMIP_data/`: ZECMIP collated co2, fgco2, nbp, tas time series; `zecmip_emis.csv`, `zecmip_tas.csv`.
- `c_data/1pctCO2_gm.nc`, `abrupt-4xCO2_gm.nc`: global means (tas, rsut, rlut, rsdt) for 29 ensemble members, 165 years; useful for TCR across a wider model set.
- `SCM_data/`: FaIR, MAGICC, CICERO-SCM flat10 runs and metrics (for the emulator line of evidence later).
- `flat10_protocol.pdf`, `flat10*.ipynb`: protocol and processing reference.

### 6.1 CMIP6 (and CMIP7 where available), fallback for gaps
Experiments: `1pctCO2`, `1pctCO2-bgc`, `1pctCO2-rad`, `historical`, `esm-hist`, `piControl`, `esm-piControl`; flat10MIP: `esm-flat10` (and `esm-flat10-zec`, `esm-flat10-cdr` optionally).

Variables (monthly → annual): `tas` (Amon), `nbp` (Lmon), `fgco2` (Omon), `co2`/`co2mass` (Amon, for esm runs), optionally `cLand`/`cVeg`/`cSoil`, `dissicos` not needed. For the ocean EC predictors use upstream Terhaar values unless recomputing.

Access, in order of preference:
1. Pangeo CMIP6 on Google Cloud (zarr, `gs://cmip6`, catalogue `https://storage.googleapis.com/cmip6/pangeo-cmip6.json`) via `intake-esm`.
2. ESGF (`intake-esgf` or `esgf-pyclient`).
3. flat10MIP: global means are in esm-deck (§6.0); ask Ben if gridded fields are needed.

Still needed from ESGF/Pangeo: `historical` and `esm-hist` global-mean tas, nbp, fgco2 (and co2 for esm-hist) for the esm-deck models plus EC-Earth3-CC and CMCC-ESM2, for Workstream B's historical observables.

Baselines: subtract piControl (or esm-piControl) over the corresponding period; drift-correct carbon fluxes using a linear fit to the control.

Target model list: the 13 CMIP6/CMIP6+ models in upstream (ACCESS-ESM1-5, BCC-CSM2-MR, CanESM5, CESM2, CMCC-ESM2, CNRM-ESM2-1, EC-Earth3-CC, GFDL-ESM4, IPSL-CM6A-LR, MIROC-ES2L, MPI-ESM1-2-LR, NorESM2-LM, UKESM1-0-LL), plus MRI-ESM2-0 and GISS-E2-1-G where ocean data exist, plus any other model with BGC and COU runs.

### 6.2 Observations
- Global Carbon Budget (latest release): fossil, LUC emissions, atmospheric growth, ocean and land sinks.
- NOAA GML global CO2 growth rate.
- HadCRUT5 (for TCR or attribution only).
- Upstream Terhaar observational values (§3.1) and Zechlau et al. 2022 IAV values.

### 6.3 Reference values for checks
- Arora et al. 2020 CMIP6 means: β_L 0.97 ± 0.40, β_O 0.79 ± 0.07 PgC ppm⁻¹; γ_L −45.1 ± 50.6, γ_O −17.2 ± 5.0 PgC K⁻¹; TCRE 1.77 ± 0.37 K/EgC (11 models).
- JF2020 Table 2 CMIP6: AF 0.532 actual vs 0.534 reconstructed; TCRE 1.78 vs 1.75.

---

## 7. Proposed repo layout

```
tcre-evidence/
  HANDOVER.md               this file
  CLAUDE.md                 short operating notes for Claude Code (see §9)
  README.md
  pyproject.toml            package `tcre_evidence`, Python ≥ 3.11
  environment.yml           pymc, arviz, xarray, intake-esm, gcsfs, zarr, netCDF4, pandas, scipy, matplotlib, joblib, pytest
  external/TCRE_Assessment  git submodule pinned to 98c74cb
  external/esm-deck         git submodule pinned to eef04dc (global-mean C4MIP 1pct BGC/COU/RAD, flat10)
  data/
    raw/                    gitignored; download scripts write here
    derived/                small per-model tables (csv/nc), committed
    upstream_snapshot/      Google Sheet + upstream DATA copies, with retrieval date
    obs/                    GCB, NOAA, etc. (or download scripts)
  src/tcre_evidence/
    cmip/                   fetch, annualise, baseline, compatible emissions
    feedbacks.py            β, γ (BGC/COU 2×2, RAD), ε, at year 70 / 140
    identity.py             AF, TCRE reconstruction, closure checks
    upstream.py             thin wrappers calling external/TCRE_Assessment builders
    perfect_model.py        LOO harness shared by A4 and B
    attribution.py          variance attribution under 3 conventions
    synthetic.py            A6 synthetic process experiments
    tcr_r.py                Workstream B model(s)
    scoring.py              coverage, CRPS, log score
    plotting.py
  scripts/                  CLI entry points, one per analysis, writing to results/
  notebooks/                thin; call package functions only
  results/                  figures + tables (gitignore large traces; keep summaries)
  tests/                    pytest: identity closure, regression values in §2.3, scoring
```

---

## 8. Milestones

| # | Deliverable | Rough effort |
|---|---|---|
| M0 | Repo scaffold, env, submodule pinned, upstream notebook reproduces 1.766 [0.927, 2.737] (within sampling noise) | 0.5 day |
| M1 | Upstream data snapshot; regression test of §2.3 table | 0.5 day |
| M2 | Reader for esm-deck C4MIP files (sign/format harmonisation, piControl baselines); β/γ/ε/TCRE/R per model at 2× and 4×; closure test (A1); reproduce §2.4 | 1–2 days |
| M3 | A2 convention sensitivity + attribution figure; A3 covariance rerun | 2 days |
| M4 | A4 perfect-model harness + results | 2–3 days |
| M5 | B estimator + historical observables + transfer perfect-model test | 3–4 days |
| M6 | A5, A6, A7; comparison figure; short write-up of results (results/SUMMARY.md) | 3 days |

Do M0–M1 before anything else. A1, A2 and the flat10 part of B can run on esm-deck alone. Only the historical observables for B need ESGF/Pangeo; if that is slow, proceed with A3–A6 on upstream inputs meanwhile.

---

## 9. Operating notes (also put a condensed version in CLAUDE.md)

- Do not modify `external/`. Wrap upstream functions; if behaviour must change for a test, re-implement in `src/` and document the difference.
- Every number that appears in a figure or summary must be produced by a script in `scripts/` and reproducible from `data/derived/` with a fixed seed.
- Record PyMC/ArviZ versions and sampler settings with every result; check r̂ < 1.01 and ESS; report divergences.
- Keep β/γ convention explicit in every variable name or attribute (`gamma_coubgc`, `gamma_rad`, `eps`).
- Units: β PgC ppm⁻¹; γ PgC K⁻¹; α K ppm⁻¹; TCRE K EgC⁻¹ (equivalently K per 1000 PgC). Watch upstream's AF/k naming.
- Small-ensemble caution: report per-model results and leave-one-family-out variants (shared code bases: e.g. NorESM/CESM, UKESM/ACCESS land components, IPSL/CNRM ocean heritage).
- Ask Ben before: contacting anyone, pushing to a public remote, or framing results as conclusions about the upstream authors' work. Results summaries should state findings about methods, with the numbers.

---

## 10. Key references

- Arora, V. K. et al. 2013, J. Clim. 26, 5289 (CMIP5 feedbacks, both conventions).
- Arora, V. K. et al. 2020, Biogeosciences 17, 4173 (CMIP6 feedbacks).
- Gregory, J. M. et al. 2009, J. Clim. 22, 5232 (quantifying carbon-cycle feedbacks).
- Hall, A. et al. 2019, Nat. Clim. Change 9, 269 (emergent-constraint credibility).
- Jones, C. D. & Friedlingstein, P. 2020, ERL 15, 074019 (TCRE decomposition).
- Marvel, K. et al. 2026, Front. Clim. 8, doi:10.3389/fclim.2026.1839130.
- Ribes, A. et al. 2021, Sci. Adv. 7, eabc0671 (Bayesian/kriging constraints).
- Sanderson, B. M. et al. 2025, GMD 18, 5699 (flat10MIP).
- Schwinger, J. et al. 2014, J. Clim. 27, 3869 (ocean feedback nonlinearity).
- Séférian, R., Ribes, A., Qasmi, S. 2025, EGU25-12981 (historical AF constraint; abstract only).
- Stein, U. & Alpert, P. 1993, J. Atmos. Sci. 50, 2107 (factor separation).
- Terhaar, J. et al. 2022 (ocean anthropogenic carbon uptake constraint; AMOC/SSS/CUC). Confirm exact citation from upstream.
- Zechlau, S. et al. 2022 (IAV constraint on tropical land γ, CMIP5/6). Confirm exact citation from upstream.
- Zickfeld, K. et al. 2011, J. Clim. (nonlinearity of carbon-cycle feedbacks).
- WCRP TCRE and ZEC Assessment report, Vienna 2025: https://www.wcrp-climate.org/WCRP-publications/2025/TCRE-Report-Vienna-2025-final.pdf

Verify all references against Crossref before citing in any manuscript.
