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

- Upstream pipeline variants (seed 0; `results/m3_covariance.json`, `fig_m3_covariance.png`):

  | variant | TCRE mean ± sd | 5–95 % | 95 % HDI |
  |---|---|---|---|
  | U, as published | 1.77 ± 0.47 | 1.07–2.58 | 0.90–2.71 |
  | J, joint land model (correlated per-model biases) | 1.77 ± 0.46 | 1.07–2.58 | 0.95–2.78 |
  | U + TCR–carbon rank correlation from the deck | 1.74 ± 0.38 | 1.10–2.35 | 1.02–2.51 |
  | J + TCR–carbon rank correlation from the deck | 1.73 ± 0.38 | 1.10–2.34 | 1.01–2.48 |
  | U + TCR–carbon rank correlation from the hard-coded α order | 1.79 ± 0.51 | 1.05–2.68 | 0.84–2.76 |

- Fitting β_L and γ_L jointly does not change the result. The per-model biases are correlated (posterior mean −0.46 between β_L and tropical γ_L, −0.39 with extratropical γ_L), but the framework's reported truth (baseline × Πη + Σδ) excludes the biases, so the posterior correlation of the land truths is 0.01 to 0.04 either way. In this structure the cross-model β–γ covariance does not reach the TCRE estimate. (Joint fit: r̂ 1.001, bulk ESS > 2500, 8 divergences in 4000 draws.)
- Retaining the TCR–carbon covariance does matter: imposing the deck's cross-model rank correlations (TCR–β −0.15, TCR–γ +0.56 [90 % CI 0.06, 0.84], β–γ −0.44; n = 11) on the TCR prior and carbon posteriors lowers the 95th percentile from 2.58 to 2.35 and the HDI upper bound from 2.71 to 2.51. With the misaligned hard-coded α the same procedure moves the upper tail the other way.
- Which covariances survive substituting the assessed TCR: only TCR–γ is distinguishable from zero across 11 models, and only at the 90 % level. Imposing it assumes the cross-model relationship reflects a physical link that also holds for the real world; it is reported as a sensitivity, not a preferred estimate.

## M4 / A4. Perfect-model test of the upstream framework (`scripts/m4_perfect_model.py`)

- Setup: 9 target models with an actual 1pctCO2 TCRE and ocean predictors (ACCESS-ESM1-5, CESM2, CanESM5, GFDL-ESM4, IPSL-CM6A-LR, MIROC-ES2L, MPI-ESM1-2-LR, NorESM2-LM, UKESM1-0-LL). For each, the model (or its family) is removed from the land, IAV and ocean training data and from the paired nitrogen experiments; its own γ_IAV, AMOC, SSS and CUC replace the observations (with the real observational σ; IPSL has no simulated IAV, so that constraint is made uninformative); the three upstream component models are refitted unchanged; α is the held-out model's own ΔT/ΔCa, for the framework and the baselines alike, so the test isolates the carbon-cycle part. Baselines: Gaussian of the training models' R times the held-out TCR ("raw ensemble"), and a Student-t random-effects predictive. Because R varies by only 6 % across models, these baselines are narrow and hard to beat. All 36 component fits: r̂ ≤ 1.009, no divergences.
- Two framework outputs are scored. (a) The reported truth (baseline × Πη + Σδ, "all processes on"), which is what the framework reports for the real world. (b) A flag-conditional predictive for a new model with the held-out model's own process flags, which tests the process regression out of sample.

  | leave one model out, TCRE | CRPS | CRPSS vs raw | better than raw | 66 % / 90 % coverage | mean bias | mean PIT |
  |---|---|---|---|---|---|---|
  | raw ensemble | 0.057 | 0 | | 0.78 / 1.00 | −0.03 | 0.61 |
  | random effects | 0.057 | | | 0.89 / 1.00 | | |
  | framework, reported truth | 0.278 | −3.9 | 0 / 9 | 0.22 / 0.89 | +0.43 | 0.11 |
  | framework, own flags | 0.111 | −0.96 | 0 / 9 | 0.89 / 1.00 | +0.20 | 0.46 |

  Leave one family out gives the same picture (CRPSS −3.6 and −0.84).
- (a) The reported truth lies above the held-out model's TCRE for all 9 models, by 0.43 K EgC⁻¹ on average (about 25 %; CanESM5 +1.31, UKESM1-0-LL +0.70, IPSL-CM6A-LR +0.54). This is largely by construction: the reported truth adds every process effect, including a permafrost effect whose prior mean is about −50 PgC K⁻¹, to every model. If a CMIP model were the real world, the framework would not recover its TCRE.
- (b) With the model's own flags the predictive is centred (mean PIT 0.46) but much wider than the ensemble (median sd 0.43 against 0.11 K EgC⁻¹) and has no skill. The widest case is CESM2 (5–95 % 1.30–4.76), one of the three models with the permafrost flag.
- Components (leave one model out, CRPSS against the raw training ensemble of the same input values): β_L −0.61 (reported truth) and −0.16 (own flags); γ_L −0.98 and −1.18; β_O −0.05; γ_O −0.10. The ocean emergent constraint is roughly neutral out of sample; the land process models lose skill.

## M5 / Workstream B. TCRE = TCR × R (`scripts/m5_tcr_r.py`)

- Data: CMIP6 historical global means fetched from Pangeo/ESGF with QC (`data/derived/cmip_globalmean/`, availability table there); GCB 2025, NOAA GML, HadCRUT 5.2 in `data/obs/`. Historical observables 1960–2014, defined so that model and observed quantities match (model nbp includes land-use change, so the model airborne fraction is that of fossil emissions): cumulative ocean sink (obs 113 ± 22 PgC; models 77–112), cumulative net land flux (obs 36 ± 27 PgC as the budget residual; models −15 to 110), fossil airborne fraction (obs 0.539 ± 0.027; models 0.47–0.69).
- R varies little across models (0.884 EgC⁻¹, CV 6 %) and is uncorrelated with TCR (r = 0.07, n = 11). Correlations of R with the observables are weak: ocean sink +0.48 (90 % CI −0.09 to +0.82, n = 10), net land −0.10, airborne fraction 0.00.
- Perfect-model test (leave one model out, the held-out model's own H as pseudo-observation, TCR known): no observable improves on the plain ensemble. CRPS skill against the ensemble is +0.00 (ocean), −0.11 (land), −0.17 (airborne fraction); with perfect observations at most +0.03. Conditioning on TCR as well costs about 0.3 in skill (estimating a 3 × 3 covariance from 10 models adds noise when the TCR–R correlation is near zero). Leave-one-family-out gives the same picture.
- Estimates with the upstream TCR prior N(1.7, 0.36): TCR × R unconstrained 1.50 [0.96, 2.07] (mean, 5–95 %); with the ocean-sink constraint 1.54 [1.02, 2.07]. With an AR6-like TCR N(1.8, 0.365): 1.59 [1.03, 2.18] and 1.63 [1.09, 2.17].
- Pathway: TCRE(1pctCO2) − TCRE(flat10) = +0.05 ± 0.09 K EgC⁻¹ (7 models). The historical attribution in the upstream DAMIP notebook (quoted, not rerun; its attribution module is not in the pinned commit) is 1.386 ± 0.087.
- With the same TCR prior, the upstream posterior mean (1.77) implies an effective carbon factor near 1.04 EgC⁻¹, about 2.5 ensemble standard deviations above the CMIP6 R (0.884 ± 0.054). A5 identifies where this comes from.

## M6 / A5. Prior versus data (`scripts/m6_prior_sensitivity.py`)

- Prior predictive against posterior for the reported truths (`results/m6_prior_info.csv`):

  | quantity | prior median (IQR) | posterior mean ± sd | CMIP6 inputs mean ± sd | KL(post ‖ prior) |
  |---|---|---|---|---|
  | β_L | 0.82 (0.81) | 0.74 ± 0.27 | 1.12 ± 0.43 | 0.30 |
  | γ_L | −14 (211) | −52 ± 41 | −26 ± 26 | 0.94 |
  | β_O | 1.01 (0.68) | 1.00 ± 0.13 | 0.90 ± 0.07 | 1.06 |
  | γ_O | −8.9 (3.6) | −9.7 ± 2.4 | −9.0 ± 2.7 | 0.05 |
  | TCRE | 1.46 (0.80) | 1.77 ± 0.48 | | 0.22 |

- γ_O is prior-dominated: its prior (standardised N(0, 1)) is the CMIP distribution itself, and the posterior barely departs from it (IQR ratio 0.89). β_O narrows strongly but is centred on the prior median of exp(θ), 1.0, about 1.4 CMIP standard deviations above the CMIP mean; with tighter regression priors (σ 1 instead of 10) it moves to 1.16.
- Relative to the TCRE implied by the CMIP6 input means (1.55 K EgC⁻¹ at TCR 1.7 K), the posterior means give 1.74 (`results/m6_shift_decomposition.json`). Of the change in the denominator k + β + αγ, 86 % comes from the lower β_L (the reported truth applies all multiplicative process factors, chiefly η_nitrogen ≈ 0.7, to every model), 36 % from the more negative γ_L (process deltas, chiefly permafrost), and −24 % from the higher β_O.
- Sensitivity (one prior changed, other components at default; `results/m6_prior_sensitivity.csv`; all fits r̂ ≤ 1.007, no divergences): TCRE mean 1.77 by default, 1.67 with no permafrost effect, 1.71 and 1.80 with permafrost prior medians of −10 and −60 PgC K⁻¹, 1.89 with wide η priors (σ 1), 1.80 without the nitrogen pseudo-observations, 1.72 with tight ocean regression priors. The 95th percentile ranges from 2.39 (no permafrost) to 2.85 (wide η). Structural-variance and LKJ priors have small effects (≤ 0.02 in the mean).

## M6 / A6. Process identifiability (`scripts/m6_identifiability.py`)

- In the lookup table the permafrost flag is set for exactly three models, CESM2, NorESM2-LM and CMCC-ESM2, which all use the CLM land model and also carry the nitrogen and fire flags. The permafrost effect is therefore not separable from a CLM-family offset in this design.

- Synthetic truth experiments (12 replicates per scenario; data generated on the real 13-model lookup table with model-level process scatter at a quarter of the upstream σ_process_struct, structural biases of 8 PgC K⁻¹ and noise of 5 PgC K⁻¹; refitted with the upstream additive process model and the upstream γ_L priors; all fits r̂ ≤ 1.01, 1 divergence in 36 fits; `results/m6_identifiability.csv`, `fig_m6_identifiability.png`):

  | scenario | δ_N (truth → mean estimate, 90 % coverage) | δ_permafrost | δ_fire | δ_veg |
  |---|---|---|---|---|
  | S1 data from the model's own structure | 10 → 12.4 (100 %) | −30 → −36.4 (100 %) | −5 → −0.5 (100 %) | 5 → 0.1 (100 %) |
  | S2 CLM-family offset ±15, no permafrost effect | 10 → 13.2 (100 %) | 0 → −29.4 (0 %) | −5 → 0.4 | 5 → 0.3 |
  | S3 ε channel, no direct nitrogen effect on γ | 0 → 8.5 (0 %) | −30 → −35.0 | −5 → −0.4 | 5 → −0.1 |

- δ_permafrost returns about −30 PgC K⁻¹ (the prior median) whether the truth is −30 or 0, and absorbs a CLM-family offset; its prior (negative-definite) excludes zero, so a zero effect can never be covered. With three flagged models that share one land model, the design cannot identify it.
- δ_fire and δ_veg shrink to zero with posterior sd equal to the prior sd (about 9 against 10): not identified, though nominally covered because the intervals are wide.
- δ_N is recovered when the data follow the model's structure, but when the nitrogen effect acts on β and reaches γ only through the COU−BGC interaction term (S3), it is attributed to a direct δ_N of +8.5 ± 5.4 with zero coverage. The separate η (β) and δ (γ) parameterisation double-counts the ε channel.
- Paired-experiment ν evidence (computed in upstream `data_utils.load_process_evidence` but not used in the notebook): ν_N = +7.2 (ACCESS-ESM CN−C) and +0.2 (UKESM1-1 ctrl−nonlim), ν_fire = −13.0, ν_veg = +1.2 PgC K⁻¹. Added as pseudo-observations (σ 5) on the summed tropical and extratropical δ_p, they change the real-data posterior of δ_N from 18.4 ± 10.2 to 7.0 ± 2.4 and of δ_fire from 0.6 ± 12.9 to −11.2 ± 4.6, and γ_L from −50.0 to −62.7 PgC K⁻¹ (about +0.04 K EgC⁻¹ in TCRE, first order at TCR 1.7 K). δ_permafrost is unchanged (−29 ± 22): no paired experiment informs it. These experiments come from models already in the ensemble (ACCESS, UKESM).

## M6 / A7. State dependence, 2× to 4×CO₂ (`scripts/m6_state_dependence.py`)

- β falls with forcing level (β_L 1.22 → 1.00, β_O 0.92 → 0.80; 10/11 and 11/11 models), so TCRE at 4× reconstructed from 2× parameters is too low by 0.17 ± 0.08 K EgC⁻¹ (−11 %) under either convention.
- γ_O(COU−BGC) roughly doubles from 2× to 4× (−8.8 → −16.7, 11/11 models) while γ_O(RAD) changes little (−3.7 → −5.0): the drift of the C4MIP ocean γ is almost entirely the interaction term. ε_O / (ΔCa·ΔT) is nearly level-independent (−0.016 and −0.014 PgC ppm⁻¹ K⁻¹; r = 0.91 across models between levels).
- γ_L is more negative at 4× in 9 of 11 models (mean change −10 PgC K⁻¹; upstream CSV −13).

## Comparison (`scripts/m6_comparison.py`)

<!-- COMPARISON -->

## Caveats and open items

<!-- CAVEATS -->
