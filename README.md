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
4. `04_wiemip.ipynb`: what the WIEMIP 1pctCO2 land runs already show (β_L, γ_L, climate pattern,
   paired factorial process effects, detection power); built by `scripts/make_notebook_04.py` from
   `results/wiemip/`.

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
| `docs/` | `guide.md` (concepts), `new_ensemble.md` (applying the tests elsewhere), `off_hub.md` (WIEMIP data and runs without the JupyterHub). |
| `tests/` | `python -m pytest`; add `-m slow` for the MCMC tests. |
| `STATUS.md` | Current state, setup on a new machine, next steps and open items. |
| `HANDOVER.md`, `CLAUDE.md` | Original project brief and working rules. |

