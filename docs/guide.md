# Guide: the TCRE process framework and how this repository tests it

This guide explains what the framework does, where it lives in the code, and how to read each test,
with the statistical terms spelled out. The notebooks in `notebooks/` carry out the analysis, and the
technical note in `note/` has the wider results.

## 1. TCRE from feedback parameters

In a concentration-driven 1pctCO2 run, the carbon taken up by land and ocean is written as a response
to CO₂ (β) plus a response to warming (γ):

    ΔC_land  = β_L ΔCa + γ_L ΔT        ΔC_ocean = β_O ΔCa + γ_O ΔT

With α = ΔT/ΔCa (warming per ppm) and k = 2.12 PgC/ppm, this gives

    TCRE = α / (k + β + αγ),     β = β_L + β_O,  γ = γ_L + γ_O.

β is in PgC/ppm and γ in PgC/K. A larger β (stronger sink) lowers TCRE; a more negative γ raises it.

**Conventions matter.** γ can be diagnosed in two ways:
- **COU−BGC:** from the coupled run minus the biogeochemically coupled run. This is what the
  development code uses.
- **RAD:** from the radiatively coupled run.

They differ by a non-linear term ε = ΔC_COU − ΔC_BGC − ΔC_RAD, which is large for some models. Variable
names in this repository always say which convention they use (`gamma_land_coubgc`, `gamma_rad`).

## 2. Interpretive models (Marvel et al. 2026)

An interpretive model is an explicit statement of how the models relate to the real world. The paper
builds up in steps:

- **Random effects.** Each model's value is the real-world value plus a model-specific error. The real
  world is then near the ensemble mean, with uncertainty set by the spread between models.
- **Correlated random effects.** Part of the error is shared by all models (common heritage, common
  tuning targets). The fraction of the error that is shared is called ρ. A shared error cannot be seen
  in the spread, so ρ comes almost entirely from its prior.
- **Process model.** Each model has a flag for each process (1 if it represents nitrogen limitation,
  fire, permafrost carbon, dynamic vegetation, and so on). A process changes the value
  multiplicatively by a factor η (used for β_L), or additively by δ (used for γ_L). The difference
  between flagged and unflagged models is what informs η and δ.
- **"All processes on."** The development code's estimate for the real world is the baseline with every
  listed process switched on, at its ensemble-average effect: U × Πη for β_L, U + Σδ for γ_L.
- **Paired experiments.** Runs of one model with a process switched on and off give direct evidence on
  η or δ, entered as extra data.
- **Emergent constraints.** A relationship across models between the target and an observable quantity
  is combined with an observation of that quantity.

In code, the process model is `build_vectorized_process_model` in
`src/tcre_evidence/upstream_patched/CMIP_process.py`, and `process_tests.build()` wraps it. In the model
each model's value is

    additive:        y_m = U + Σ_p L_mp (δ_p + e_mp) + b_common + b_m
    multiplicative:  y_m = U · exp(Σ_p L_mp (log η_p + e_mp)) + b_common + b_m

where:
- U is the value with no processes;
- e_mp lets a process act differently in each model;
- b_common is the shared bias;
- b_m is the bias specific to model m.

## 3. The development configuration

`src/tcre_evidence/upstream.py` runs it end to end (`upstream.run()`), with the inputs frozen in
`data/upstream_snapshot/`. It has three parts, fitted separately and combined draw by draw:

| Component | Model | Extra evidence | Builder |
|---|---|---|---|
| β_L | multiplicative; nitrogen, fire, veg | paired nitrogen runs (η = 0.69, 0.77) | `build_beta_land_model` |
| γ_L | additive; tropics and extratropics separately; nitrogen, permafrost, fire, veg | interannual-variability (IAV) constraint on tropical γ_L | `build_gamma_land_model` |
| β_O, γ_O | emergent constraint only (no process model) | AMOC, sea-surface salinity, chemical uptake capacity | `build_joint_ocean_model` |

TCRE is then computed from each posterior draw, with TCR drawn independently from N(1.7, 0.36) K. The
pinned commit does not construct as found; `upstream_patched/` changes only internal variable names
(`upstream.diff` shows the exact patch).

## 4. The tests and how to read them

`src/tcre_evidence/process_tests.py` runs these on any ensemble; `scripts/process_tests.py` is the
command line. The CMIP6 versions behind the note are the `m*` scripts listed in the README.

1. **Flag budget** (`flag_budget`). A least-squares regression of the quantity on the flags.
   - `r2` is the share of the spread between models that the flags explain.
   - `loo_q2` is the same share when each model is predicted from the others. If it is negative, the
     flags predict worse than the ensemble mean.
   - For CMIP6 γ_L: r2 = 0.20, loo_q2 = −0.70.
2. **Detection power** (`detection_power`). Given the scatter the flags leave unexplained, this is the
   standard error of each effect, and the number of independent models needed to detect an effect of a
   given size at two standard errors. It is cheap, and it is the quickest way to see whether a flag can
   be informative at all.
3. **Prior against posterior** (`prior_vs_posterior`). One fit, comparing each term with its prior.
   - `width_ratio` (posterior over prior 90 % interval width) near 1 and `shift` near 0 mean the
     ensemble has not informed the term; its value in the result then comes from the prior.
   - A narrower posterior is not proof of identification: the data may only rule out the prior's tail.
     For the CMIP6 permafrost δ the width halves, yet the synthetic test shows the estimate stays at
     the prior median whatever the true effect.
4. **Synthetic recovery** (`synthetic_recovery`). Pseudo-ensembles are made on the real flag table with
   known effects and the same scatter as the real data, then refitted.
   - A term is identified if its estimate follows the truth.
   - If the estimate stays at the prior mean whatever the truth, the term is not identified.
   - `cov90` is the share of replicates whose 90 % interval contains the truth.
5. **Leave one model out** (`leave_one_out`). Each model in turn plays the real world and is removed from
   training (with `groups`, its whole family is removed). Three predictions of its value are scored:
   - `all_processes_on`: the real-world estimate;
   - `own_flags`: the prediction for a model with its flags;
   - `raw_ensemble`: the training models' mean and spread.

   `crps_skill` below 0 means worse than the raw ensemble. The real-world estimate is not expected to
   match a model that lacks some processes; the test shows how far the adjustment moves the answer, and
   that nothing in the ensemble checks it.

## 5. Practicalities

- **Sampler checks.** Each fit reports `rhat_max` (should be below 1.01), effective sample size and
  divergences. A few divergences in many fits are usually harmless; many in one fit are not.
- **Runtimes.** A single fit of one quantity takes a few minutes; the full γ_L component (tropics,
  extratropics and the IAV constraint) about ten. `--quick` (two chains of 500 draws) is for trying things out, not for reported numbers.
- **Parallel fits.** Set `OMP_NUM_THREADS=1` before starting Python when using `--jobs` or
  `TCRE_JOBS` (the scripts do this). Forked workers with multithreaded BLAS can stall.
- **Seeds.** Reported numbers use fixed seeds; change `--seed` to check that a result is not a
  sampling accident.

## 6. Glossary

| Term | Meaning |
|---|---|
| α, β, γ, k | warming per ppm; carbon uptake per ppm; carbon uptake per K; atmospheric conversion (2.12 PgC/ppm) |
| R | the carbon factor TCRE/TCR, used in the TCR × R line |
| η, δ | multiplicative and additive process effects |
| U | value with no processes ("unscaled" in the code) |
| ρ | share of structural error common to all models |
| σ_struct | total structural spread between models |
| IAV | interannual variability; here the sensitivity of the CO₂ growth rate to tropical temperature |
| CRPS | continuous ranked probability score; lower is better; skill is relative to the raw ensemble |
| r̂, ESS | convergence and effective sample size diagnostics of the sampler |
| HDI | highest-density interval |
