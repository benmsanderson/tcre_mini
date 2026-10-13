# Results summary: TCRE lines-of-evidence testbed

Status as of 2026-10-13. Every number below is produced by the script named in its section, from files in `data/` with fixed seeds; the JSON/CSV files in `results/` hold the full values. Findings describe the behaviour of methods and inputs; they are not conclusions about the authors of the upstream framework, and they should be checked with Ben before being communicated.

Conventions: β in PgC ppm⁻¹, γ in PgC K⁻¹, α in K ppm⁻¹, TCRE in K EgC⁻¹, R = TCRE / TCR in EgC⁻¹. "COU−BGC" is the C4MIP convention (Arora et al. 2020 2×2 solve); "RAD" uses 1pctCO2-rad; ε = ΔC_COU − ΔC_BGC − ΔC_RAD. Default diagnosis: esm-deck C4MIP global means, 20-year mean over simulation years 61–80 (2×CO₂), temperature and carbon fluxes relative to a linear fit of piControl.

## Overview

<!-- OVERVIEW -->

## M0. Upstream reproduction (`scripts/m0_reproduce_upstream.py`)

- At the pinned commit (98c74cb) the land builders do not construct (TypeError in `build_beta_land_model`; duplicate variable names in `build_gamma_land_model`; `add_emergent_constraint` reads a variable that is never created). The saved notebook output lists variable names that only exist in a version with `var_name`-prefixed internals, so the published headline was produced with code that differs from the commit.
- A minimal patch (names only; priors, likelihoods and structure unchanged; exact diff in `src/tcre_evidence/upstream_patched/upstream.diff`) reproduces the headline: TCRE mean 1.771 (sd across 3 seeds 0.002), 95 % HDI [0.903, 2.714], against the notebook's 1.766 [0.927, 2.737]. All fits r̂ < 1.004, bulk ESS > 1600, no divergences (PyMC 5.28.5, 4 × 1000 draws, target_accept 0.99).

## M1. Upstream inputs and the §2.3 covariance table (`scripts/m1_covariance_check.py`)

- The HANDOVER §2.3 numbers are reproduced exactly by bootstrap recombination of the 11 land∩ocean models (regression tests in `tests/test_covariance.py`).
- The hard-coded α array in the notebook's CMIP reference does not follow the `np.intersect1d` model order it is paired with: it follows the order of the 11 models that have α in the sheet (which include BCC-CSM2-MR and CNRM-ESM2-1, not in the common set). Six of the nine checkable models receive another model's α; CMCC-ESM2 and EC-Earth3-CC have no α.
- With α matched by name (9 models): corr(α, α·γ) changes from −0.53 to +0.28 and corr(α, β_O) from +0.52 to −0.02; paired recombination gives a narrower upper tail (95th percentile 2.35) than upstream-style recombination (2.58). The §2.3 statements that high-TCR models lose more carbon and that independent recombination truncates the upper tail depend on the misalignment.
- The notebook's posteriors use the sheet values, its CMIP reference uses the CSV/Terhaar files; these differ by up to about 4 PgC K⁻¹ in γ_L and 1.5 PgC K⁻¹ in γ_O.

## M2 / A1. Harmonised feedbacks and closure (`scripts/m2_feedbacks_closure.py`)

- Reader for the esm-deck C4MIP files with an explicit per-model sign, unit and format table and assertions (`src/tcre_evidence/c4mip.py`). Data issues handled: mixed flux signs; UKESM forced runs in K but piControl in °C with a different header; CNRM BGC/RAD runs 10 years shorter; CNRM piControl CO₂ missing; MIROC (JAMSTEC) piControl CO₂ column holds the 1 %/yr ramp.
- 11 models, 2×CO₂: TCRE 1.74 ± 0.41 K EgC⁻¹; TCR 1.96 (CV 23 %); R = 0.884 EgC⁻¹ (CV 6 %); AF 0.533.
- HANDOVER §2.4 is reproduced with its first-10-year baseline (TCRE 1.70 ± 0.36; RAD reconstruction bias +0.009 ± 0.023 in AF and +0.03 ± 0.07 K EgC⁻¹ in TCRE; ACCESS ε_L 116 PgC; CNRM ε_O +49 PgC). With piControl drift correction the CNRM ocean exception disappears (ε_O −10 PgC; its control ocean flux is −0.8 PgC yr⁻¹) and BCC ε_L changes sign (+25 to −33 PgC). ε is sensitive to drift correction.
- COU−BGC reconstruction closes to machine precision (bookkeeping, as expected).
- Convention of upstream inputs: our COU−BGC values correlate with the upstream land CSV at r = 0.96 (β_L) and 0.97 (γ_L), and with the Terhaar ocean file at 0.97 (β_O) and 0.99 (γ_O), against 0.31 for our RAD γ_O. The Terhaar file is COU−BGC despite its "1pctCO2-BGC" attribute.
- Closure of upstream inputs against each model's actual 1pctCO2 TCRE: sheet β, γ with our α reproduce TCRE to 0.012 ± 0.026 K EgC⁻¹; with the sheet's own α the error is 0.03 ± 0.14 (up to 0.23 for IPSL and ACCESS). The notebook's CMIP reference (CSV/Terhaar with hard-coded α) misses by up to 0.89 K EgC⁻¹ (GFDL), sd 0.47.

## M3 / A2. Convention sensitivity and attribution (`scripts/m3_convention_attribution.py`)

- TCR dominates: 89 % of the cross-model variance of ln TCRE comes from α (94 % with the first-10-year baseline), consistent with HANDOVER §2.5.
- Carbon-term variance (S = β + αγ, exact per-model decomposition) under three allocations of ε, 2×CO₂: the β terms carry 71–84 % by independent shares in all three, and land carries 84–99 %. Covariance shares exceed one for β because β and αγ covary negatively; the C4MIP allocation gives the largest β share (1.41) and COU−RAD the smallest (1.06). The "β dominates" ranking survives the change of convention; the size of the β/γ split does not.
- RAD-based reconstruction at 2×CO₂: AF +0.004 ± 0.019, TCRE +0.018 ± 0.062 K EgC⁻¹; 8 of 11 models within 0.05; ACCESS-ESM1-5 +0.16. The bias is model-specific rather than systematic, and AF is not systematically too low.
- ε_O is negative in all 11 models (−8.8 ± 3.5 PgC at 2×, −50 ± 12 PgC at 4×), consistent with Schwinger et al. (2014). ε_L is +14 ± 38 PgC at 2× (7 of 11 positive) and +48 ± 142 PgC at 4×.

## M3 / A3. Recombination and covariance (`scripts/m3_covariance_rerun.py`)

- Deck models with α, β, γ from the same runs: corr(β_L, γ_L) = −0.82, corr(β, αγ) = −0.84, corr(β_O, γ_O) = −0.71, corr(α, αγ) = +0.42, corr(α, β_O) = −0.28. Independent recombination doubles the sd of AF (0.060 against 0.030 paired, a factor of four in variance, as estimated in HANDOVER §2.3) and widens both TCRE tails: 5–95 % [1.21, 2.65] against [1.15, 2.52] paired. With consistent inputs, discarding covariance inflates rather than truncates the upper tail.

<!-- A3 -->

## M4 / A4. Perfect-model test of the upstream framework (`scripts/m4_perfect_model.py`)

<!-- A4 -->

## M5 / Workstream B. TCRE = TCR × R (`scripts/m5_tcr_r.py`)

- Data: CMIP6 historical global means fetched from Pangeo/ESGF with QC (`data/derived/cmip_globalmean/`, availability table there); GCB 2025, NOAA GML, HadCRUT 5.2 in `data/obs/`. Historical observables 1960–2014, defined so that model and observed quantities match (model nbp includes land-use change, so the model airborne fraction is that of fossil emissions): cumulative ocean sink (obs 113 ± 22 PgC; models 77–112), cumulative net land flux (obs 36 ± 27 PgC as the budget residual; models −15 to 110), fossil airborne fraction (obs 0.539 ± 0.027; models 0.47–0.69).
- R varies little across models (0.884 EgC⁻¹, CV 6 %) and is uncorrelated with TCR (r = 0.07, n = 11). Correlations of R with the observables are weak: ocean sink +0.48 (90 % CI −0.09 to +0.82, n = 10), net land −0.10, airborne fraction 0.00.
- Perfect-model test (leave one model out, the held-out model's own H as pseudo-observation, TCR known): no observable improves on the plain ensemble. CRPS skill against the ensemble is +0.00 (ocean), −0.11 (land), −0.17 (airborne fraction); with perfect observations at most +0.03. Conditioning on TCR as well costs about 0.3 in skill (estimating a 3 × 3 covariance from 10 models adds noise when the TCR–R correlation is near zero). Leave-one-family-out gives the same picture.
- Estimates with the upstream TCR prior N(1.7, 0.36): TCR × R unconstrained 1.50 [0.96, 2.07] (mean, 5–95 %); with the ocean-sink constraint 1.54 [1.02, 2.07]. With an AR6-like TCR N(1.8, 0.365): 1.59 [1.03, 2.18] and 1.63 [1.09, 2.17].
- Pathway: TCRE(1pctCO2) − TCRE(flat10) = +0.05 ± 0.09 K EgC⁻¹ (7 models). The historical attribution in the upstream DAMIP notebook (quoted, not rerun; its attribution module is not in the pinned commit) is 1.386 ± 0.087.
- With the same TCR prior, the upstream posterior mean (1.77) implies an effective carbon factor near 1.04 EgC⁻¹, about 2.5 ensemble standard deviations above the CMIP6 R (0.884 ± 0.054). A5 identifies where this comes from.

## M6 / A5. Prior versus data (`scripts/m6_prior_sensitivity.py`)

<!-- A5 -->

## M6 / A6. Process identifiability (`scripts/m6_identifiability.py`)

- In the lookup table the permafrost flag is set for exactly three models, CESM2, NorESM2-LM and CMCC-ESM2, which all use the CLM land model and also carry the nitrogen and fire flags. The permafrost effect is therefore not separable from a CLM-family offset in this design.

<!-- A6 -->

## M6 / A7. State dependence, 2× to 4×CO₂ (`scripts/m6_state_dependence.py`)

- β falls with forcing level (β_L 1.22 → 1.00, β_O 0.92 → 0.80; 10/11 and 11/11 models), so TCRE at 4× reconstructed from 2× parameters is too low by 0.17 ± 0.08 K EgC⁻¹ (−11 %) under either convention.
- γ_O(COU−BGC) roughly doubles from 2× to 4× (−8.8 → −16.7, 11/11 models) while γ_O(RAD) changes little (−3.7 → −5.0): the drift of the C4MIP ocean γ is almost entirely the interaction term. ε_O / (ΔCa·ΔT) is nearly level-independent (−0.016 and −0.014 PgC ppm⁻¹ K⁻¹; r = 0.91 across models between levels).
- γ_L is more negative at 4× in 9 of 11 models (mean change −10 PgC K⁻¹; upstream CSV −13).

## Comparison (`scripts/m6_comparison.py`)

<!-- COMPARISON -->

## Caveats and open items

<!-- CAVEATS -->
