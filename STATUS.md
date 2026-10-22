# Status and handover (2026-10-23)

Where the work stands, for picking it up on the machine that hosts the WIEMIP output. Read with
`README.md` (repository map and analysis catalogue) and `CLAUDE.md` (working rules).
`HANDOVER.md` is the original brief and is kept for background.

## What is done

**CMIP6 analysis of the development configuration (M0 to M11).** All scripts have run, with fixed seeds.
Results are in `results/`, summarised in `results/SUMMARY.md`. The main findings:
- Only nitrogen (on β_L) is close to identifiable from presence flags.
- The "all processes on" estimate does not predict held-out models.
- About half the configuration's shift above CMIP comes from priors.
- Calibration practices can hide flag effects.
- The four flags explain about 20 % of the spread in γ_L.

**Technical note.** LaTeX in the `note/` submodule (`benmsanderson/tcre_bayes_paper`, commit ed2b4b2),
with figures from `scripts/note_figures.py`. It has not been circulated.

**Toolkit for any ensemble.** `src/tcre_evidence/process_tests.py` runs:
- flag budget;
- detection power;
- prior against posterior;
- synthetic recovery;
- leave one model or family out;
- an evidence table that combines them.

It wraps the development code's own builder (`upstream_patched/`). `scripts/precompute.py` caches the
sampling; `scripts/process_tests.py` is a command-line front end.

**Notebooks** (the main presentation of the process-model analysis):
- `notebooks/01_process_model.ipynb`
- `notebooks/02_what_the_ensemble_can_measure.ipynb`
- `notebooks/03_other_ensembles.ipynb` (the template for WIEMIP)

They are executed against the committed CMIP6 cache in `results/precomputed/cmip6/` (full-length runs,
all fits r̂ ≤ 1.005, no divergences).

**Data.** The C4MIP and flat10 global means previously read from the esm-deck submodule are in
`data/c4mip/` and `data/flat10/`, copied unchanged from esm-deck commit eef04dc. `data/examples/` holds the
CMIP6 β/γ and flag tables in the format the toolkit reads.

## Setting up on the new machine

```bash
git clone --recurse-submodules https://github.com/benmsanderson/tcre_mini   # private; needs access
cd tcre_mini
conda env create -f environment.yml && conda activate tcre-evidence          # or pip install -e ".[test,notebooks]"
python -m pytest                                                             # 99 fast tests, under a minute
```

- **Versions used here:** Python 3.13, PyMC 5.28.5, ArviZ 0.23.4, NumPy 2.4, pandas 3.0. PyMC must stay
  below 6 and ArviZ below 1.
- **The `note/` submodule** is private and not needed to run anything; the clone works without it
  (`git submodule update --init external/TCRE_Assessment` fetches only the public development code).
- **Parallel fits.** Scripts set `OMP_NUM_THREADS=1` themselves. Run parallel fits through `--jobs`, not
  PyMC's own multi-core sampling: chains run one after another (`cores=1`), which avoids stalls and a
  joblib pickling failure.

## Next steps for WIEMIP

1. **Build the values table** from WIEMIP global means: one row per model, columns such as
   `gamma_land_coubgc`, `beta_land_coubgc` (and, if possible, R = TCRE/TCR).
   - The formulas are in `docs/new_ensemble.md`; `tcre_evidence.feedbacks.diagnose()` implements them for
     the C4MIP ascii format.
   - A reader for WIEMIP's own file format does not exist yet. Write it in `src/tcre_evidence/` (for
     example `wiemip.py`), with an explicit per-model table of sign conventions and units and an
     assertion on signs, as `c4mip.py` does.
   - Check which experiments WIEMIP provides: COU, BGC and RAD give both γ conventions; with RAD only, use
     `gamma_rad` and say so.
2. **Build the flags table:** one 0/1 column per process, with `demography` separate from dynamic
   vegetation (`veg`), plus `family` for shared land-model code. Becky Varney is compiling this from the
   WIEMIP modellers and the WIEMIP protocol table.
3. **Explore** with `notebooks/03_other_ensembles.ipynb`: set the paths at the top and run the cheap
   tests.
4. **Precompute:**
   `python scripts/precompute.py --name wiemip --values <values.csv> --flags <flags.csv> --quantities gamma_land_coubgc:additive beta_land_coubgc:multiplicative --jobs <n>`.
   - Use `--quick` first.
   - A full run took about 2.5 hours on 4 cores for two quantities and 13 models.
   - Commit `results/precomputed/wiemip/` if the input data may be shared; otherwise commit only the
     summaries.
5. **The demography question** (raised by Rosie Fisher and Becky):
   - Are demographic models' γ_L larger in magnitude?
   - Separate demography from family (CLM-FATES against CLM big-leaf with the same host model is the
     cleanest comparison).
   - State the γ convention.
   - Show β, γ and an end-to-end quantity (TCRE/TCR or the airborne fraction), since the β and γ effects
     partly cancel.

## Open items

- **Note Section 1 qualifier** (not yet applied; it is a statement about the upstream method, so it needs
  Ben's sign-off). With synthetic scatter equal to that left unexplained in CMIP6 (28.7 PgC/K), the
  nitrogen effect on γ_L is not recovered (tracking 0.21; `results/precomputed/cmip6/gamma_land_coubgc/`).
  The M6 test behind "only nitrogen is recovered" used scatter of about 10 PgC/K. Add a
  `results/SUMMARY.md` entry for the precompute results before citing them.
- **SUMMARY.md** has no section yet for `scripts/precompute.py`. The numbers are in the cache CSVs and the
  executed notebooks.
- **Upstream contact.** A draft email to the Frontiers authors and a reply to the WIEMIP thread were
  written in conversation and are not in the repository. Ben is handling communication; do not contact
  anyone without asking him.
- **Access.** Becky needs read access to this repository to run the notebooks; she does not need the
  note.
- **Not yet tested:**
  - the IAV emergent constraint on tropical γ_L in isolation (with and without, plus leave one out);
    Liu et al. (2024, doi:10.1038/s43247-024-01504-6) find it weak in CMIP6 (R² 0.40, the same value as
    in the development inputs);
  - the prior on the regression scatter in `add_emergent_constraint` (HalfNormal(1)) is tight relative to
    the scatter in the data (about 1.8).
- **Housekeeping.**
  - Rerunning `scripts/note_figures.py` rewrites the figure PDFs even when nothing changed (embedded
    timestamps); restore them with `git checkout` if the content is unchanged.
  - The fast test `test_detection_power_matches_m8` relies on the order of effects, because the
    balanced-design search draws random designs.
