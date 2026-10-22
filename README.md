# tcre_mini: testing process-flag interpretive models for TCRE

This repository tests the Bayesian process framework that the TCRE assessment is developing
([Marvel et al. 2026](https://doi.org/10.3389/fclim.2026.1839130); development code in
`netzeroasap/TCRE_Assessment`). It asks whether flag-based process models, applied to an ensemble of
about 13 Earth system models, learn process effects from the ensemble or mainly express their priors.
It also builds a simpler line of evidence, TCRE = TCR × R.

## The analysis, in notebooks

The notebooks in `notebooks/` are the main presentation of the process-model analysis and the place to
explore it. The data, the regression, the synthetic experiments and every figure are written out in the
cells, so tables, flags, priors and effect sizes can be changed and rerun directly. Only the development
code's model builder and the sampler are imported. Sampling results are cached in `results/precomputed/`,
so the notebooks run in seconds; `REFIT = True` samples live instead.

1. `01_process_model.ipynb`: the CMIP6 land feedbacks by process flag, the process model and its priors,
   and how far the 13 models move each term away from its prior.
2. `02_what_the_ensemble_can_measure.ipynb`: the flag regression, detection power, recovery of known
   effects from pseudo-data, and whether the real-world estimate predicts held-out models.
3. `03_other_ensembles.ipynb`: the same tests on any other ensemble (for example WIEMIP), using the
   packaged functions.

The technical note in `note/` (a submodule, `benmsanderson/tcre_bayes_paper`) covers the wider set of
analyses (`m*` scripts below); `results/SUMMARY.md` lists every number and the script it comes from.
`docs/guide.md` gives the background and `docs/new_ensemble.md` the input format for other ensembles.

## Quick start

```bash
git clone --recurse-submodules https://github.com/benmsanderson/tcre_mini
cd tcre_mini
conda env create -f environment.yml && conda activate tcre-evidence   # or: pip install -e ".[test,notebooks]"
python -m pytest                           # fast tests, under a minute
jupyter lab notebooks/                     # or open the notebooks in VS Code
```

To redo the cached sampling, or to run it for another ensemble:

```bash
python scripts/precompute.py --jobs 4                     # CMIP6 example, about an hour
python scripts/precompute.py --name wiemip --values my_values.csv --flags my_flags.csv \
    --quantities gamma_land_coubgc:additive --jobs 4     # your ensemble; --quick for a first look
```

`scripts/process_tests.py` runs the same tests from the command line without the notebooks.

## Repository map

| Path | What it holds |
|---|---|
| `src/tcre_evidence/` | The package. All logic lives here; scripts are thin. |
| `notebooks/` | The analysis of the process model, as editable notebooks (see above). |
| `src/tcre_evidence/process_tests.py` | The process-model tests, for any ensemble. |
| `src/tcre_evidence/process_plots.py` | Figures for `03_other_ensembles.ipynb`. |
| `src/tcre_evidence/upstream_patched/` | The development code's model builders, patched only so that they construct (diff and reasons inside). |
| `src/tcre_evidence/upstream.py`, `upstream_data.py` | The full development configuration (land, ocean, TCRE) and its frozen inputs. |
| `scripts/` | One script per analysis; each writes to `results/` and uses a fixed seed. |
| `results/` | Outputs and `SUMMARY.md`; `results/precomputed/` holds the notebooks' cache. |
| `data/examples/` | CMIP6 input tables in the format `process_tests.py` reads. |
| `data/upstream_snapshot/` | Frozen copy of the development code's inputs (Google Sheet tabs and `DATA/`). |
| `data/c4mip/`, `data/flat10/` | C4MIP 1pctCO2 BGC/COU/RAD and flat10 global means for 11 CMIP6 models (copied from esm-deck at eef04dc; read through `tcre_evidence.c4mip`). |
| `data/derived/`, `data/obs/` | Feedbacks diagnosed from C4MIP runs; observational series (GCB, NOAA, HadCRUT5). |
| `external/TCRE_Assessment` | Development code, pinned at 98c74cb. Never edited. |
| `note/` | LaTeX technical note (submodule); figures from `scripts/note_figures.py`. |
| `docs/` | `guide.md` (concepts), `new_ensemble.md` (applying the tests elsewhere). |
| `tests/` | `python -m pytest`; add `-m slow` for the MCMC tests. |
| `STATUS.md` | Current state, setup on a new machine, next steps and open items. |
| `HANDOVER.md`, `CLAUDE.md` | Original project brief and working rules. |

## Analyses

Script names keep the milestone numbers of the original work plan (M0 to M11). Runtimes are for 4 cores.
Scripts that sample with PyMC pin BLAS to one thread; `TCRE_JOBS` sets the number of parallel fits.

| Script | Question | Note section | Runtime |
|---|---|---|---|
| `precompute.py` | Sampling for the notebooks: prior against posterior, synthetic recovery, leave one out (any ensemble) | 1, 3, 6, 7 | 1 h |
| `process_tests.py` | The same tests from the command line | 1, 3, 6, 7 | seconds to 1 h |
| `m0_reproduce_upstream.py` | Does the development configuration reproduce its own headline (1.77 K/EgC)? | Context | 20 min |
| `m1_covariance_check.py` | Are its inputs aligned and is its covariance argument right? | Methods | seconds |
| `m2_feedbacks_closure.py` | Which β/γ convention are the inputs in, and do they close to each model's TCRE? | Methods | seconds |
| `m3_convention_attribution.py` | How much does the convention matter, and which term drives TCRE spread? | Methods | seconds |
| `m3_covariance_rerun.py` | What changes if covariances between components are kept? | Implications | 40 min |
| `m4_perfect_model.py` | Does the configuration predict held-out models (full TCRE pipeline)? | 3 | 1.5 h |
| `m5_tcr_r.py` | TCRE = TCR × R, and do historical observables constrain R? | Implications | 2 min |
| `m6_prior_sensitivity.py` | Which priors set the result? | 4 | 1 h |
| `m6_shift_decomposition.py` | Which component moves TCRE away from the CMIP inputs? | 4 | seconds |
| `m6_identifiability.py` | Are known process effects recovered? (development priors) | 1 | 2 h |
| `m6_state_dependence.py` | Do 2×CO₂ feedbacks predict 4×CO₂ TCRE? | Implications | seconds |
| `m6_comparison.py` | All lines of evidence side by side | Methods | seconds |
| `m7_nitrogen_historical.py` | Does the nitrogen flag show up in β_L, R and the historical sink? | 2 | seconds |
| `m8_power.py` | How many models would detect each effect? | 6 | seconds |
| `m9_tuning_selection.py` | Can calibration practices hide or distort a flag effect? | 6 | 30 min |
| `m10_nitrogen_decomposition.py` | Does the configuration's nitrogen decomposition reproduce the R difference? | 2, 4 | 20 min |
| `m11_completeness.py` | How much do the flags miss, and what happens to an unlisted process? | 7 | 1.5 h |
| `note_figures.py` | Figures for the note, from saved results | all | seconds |
| `fetch_cmip_globalmeans.py` | Downloads CMIP6 historical global means (network; only needed for M5, M7) | | hours |

## Conventions

- **Identity:** TCRE = α / (k + β + αγ), with k = 2.12 PgC/ppm, β = β_L + β_O, γ = γ_L + γ_O.
- **Units:** β in PgC/ppm, γ in PgC/K, α in K/ppm, TCRE in K/EgC.
- **Naming:** variable names state the β/γ convention: `gamma_land_coubgc` (COU−BGC), `gamma_rad`
  (RAD), `eps` (the non-linearity ε = ΔC_COU − ΔC_BGC − ΔC_RAD). The development code's inputs are COU−BGC.
- **Development code variable `AF`** is the airborne fraction divided by k.
- **Writing:** UK English and a neutral technical register; describe methods, not people.

## Reproducibility

- Every reported number comes from a script in `scripts/` with a fixed seed.
- Every PyMC fit records r̂ (required < 1.01), effective sample size and divergences.
- PyMC is pinned to 5.x and ArviZ to 0.x, matching the API the development notebook uses.
- The development code is a snapshot at one commit. Its priors, inputs and structure may have changed
  since, and the results here describe that snapshot, not a finished assessment.

## Adding an analysis

1. Put reusable logic in `src/tcre_evidence/` (or extend `process_tests.py`), with a docstring that
   says what the function computes and in which units.
2. Add a thin script in `scripts/` with a module docstring stating the question, method and outputs. It
   should use a fixed seed and write to `results/`.
3. Record r̂, effective sample size and divergences for every fit.
4. Add the key numbers to `results/SUMMARY.md` under a heading naming the script, and a fast test in `tests/`.
