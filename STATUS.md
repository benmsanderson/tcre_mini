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
- `notebooks/04_wiemip.ipynb` (WIEMIP summary without a flag table; text in `scripts/make_notebook_04.py`)

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

**Done on the hub (2026-10-22):**
- **Environment.** The `tcre-evidence` env is at `~/.conda/envs/tcre-evidence` (Python 3.14, PyMC 5.28.5, ArviZ 0.23.4) and is registered as a Jupyter kernel; pytest passes.
- **Extraction.** `scripts/wiemip_extract.py` runs under the hub Python, which has the maintained `wiemip_registry`; the copy in `/opt/wiemip-data-processing` is older. It writes the global land carbon of all 208 1pctCO2 runs to `data/wiemip/`.
- **Feedbacks.** `src/tcre_evidence/wiemip.py` and `scripts/wiemip_values.py` write `results/wiemip/`: values tables (UKESM pattern and pattern mean), a long table and paired factorial effects.
  - ΔT is the FaIR median GMST, the same for every pattern; ΔT_BGC = 0.
  - Both directories are committed (about 1.5 MB). The repository is private; share it only with WIEMIP authors.
- **Power before the flag table.** `scripts/wiemip_power.py` writes `results/wiemip/power_*.csv` and `flags_from_factorials.csv`.
  - It compares the cross-model standard error on hypothetical designs (18 models, k flagged) with the precision of the paired factorial effects.
  - The flag budget and the regression against real flags wait for the full flag table.
- **Known data issues.**
  - The JULES noFire_noNitrogen run has non-finite values (dropped).
  - The BEPS ctrl drifts by 62 PgC per century.
  - The open `data-quality` issues on `WIEMIP/wiemip-data-registry` have not yet been read (private repository, no `gh` on the hub).

**Done on the hub (2026-11-15):**
- **Paired factorials, hierarchical.** `process_tests.build(paired_spread=...)`: `None` keeps the development
  code's treatment (each factorial observes the mean effect with a fixed sigma, so n values narrow the
  posterior to sigma/sqrt(n) however much they disagree); `"estimate"` treats each factorial as one model's
  own effect and estimates the between-model scatter tau_p. `"struct"` uses the prior process scatter and is
  not usable for beta (`sigma_mult_struct` = 0.01 against an observed 0.1 to 0.6 in log eta).
  `wiemip.paired_by_model` gives one value per model (the three patterns are not independent).
- **Flag table.** `data/wiemip/process_inclusion_sheet.csv` is the 'WIEMIP process inclusion' tab of the
  WIEMIP spreadsheet as downloaded on 2026-11-15, unchanged. `wiemip.read_flag_sheet` harmonises names, sets
  '?' to unknown and fills it with 1 where the model ran noX; no cell contradicts the factorials. Resolved
  flags: `results/wiemip/flags_resolved.csv`. No demography or family column yet.
- **Re-extraction (2026-11-15).** The registry (code hash 27460c3b29d530f9) no longer exposes BEPS,
  DVM_DOS_TEM, JULES, VISIT_UT or ORCHIDEE_MICT_CALIPSO, although their files are under `/srv/wiemip`
  (JULES has been reorganised as `JULESwiemipV2_*` with several permafrost variants). 182 runs, against 208
  on 2026-10-22; the earlier extraction is in git history (commit 303c67a). Models still exposed are
  unchanged except the re-uploaded LPJmL6 (beta) and LPX_Bern (gamma). New: CARDAMOM_JPL, noDynVeg for
  LPJ_EOSIM and LPJ_GUESS (flag `veg`).
- **Re-extraction (2026-11-17).** JULES is exposed again (baseline, noDynVeg, noFire_noNitrogen; beta 1.20,
  gamma -34.0, as on 2026-10-22); BEPS, DVM_DOS_TEM, VISIT_UT and ORCHIDEE_MICT_CALIPSO are still missing.
  206 runs. The LPJ_EOSIM noDynVeg ctrl drifts by 133 PgC/century (anomalies are taken from its own ctrl,
  but one of three veg factorials rests on it).
- **Comparison.** `scripts/wiemip_paired_spread.py --ensemble {cmip6,wiemip}` writes
  `results/wiemip/paired_spread[_wiemip].csv`. On WIEMIP (14 models with complete flags and values), with tau
  estimated, only fire stays constrained; nitrogen on gamma_L and veg return close to prior width. The
  permafrost factorials give about zero effect on gamma_L (-0.4 PgC/K from two models; +13 before the
  LPX_Bern re-upload), a value the NegLogNormal prior (median -30) excludes, so the narrow upstream
  posterior just below zero reflects the prior's bound, not the data. In WIEMIP the same models give the regression values and
  the factorials, so the two are not independent. The upstream every-pattern-row gamma fit has 2
  divergences on 2026-11-15; on 2026-11-17 three WIEMIP gamma fits have 1 or 2 divergences (upstream rows,
  upstream models, estimate; r-hat <= 1.004).
- **Spreadsheet differences to raise with Becky** (not yet raised): JSBACH beta 1.30 in the sheet against
  1.10 from the registry; CARDAMOM_JPL beta 0.56 against 0.95; JULES pairs mix backgrounds (fire is
  baseline vs noFire, nitrogen and veg use noFire as "on"); the sheet's JULES baseline (beta 0.94, gamma
  -59.7) is not the registry's baseline (1.20, -34.0).
- **Power tables** regenerated on the new extraction (16 baseline models on 2026-11-17). `wiemip_power.py` now skips
  designs with a flag on every model (they never reach full rank, and `design_se_factor` looped forever).
- **Stale:** `notebooks/04_wiemip.ipynb` and the text in `scripts/make_notebook_04.py` quote numbers from the
  2026-10-22 extraction (18 models, fire +5.1 +/- 2.6 PgC/K); regenerate and re-execute before sharing.
- **Recovery and leave one out with factorials.** `synthetic_recovery(paired_design=...)` draws pseudo
  factorials with known between-model scatter (`synthetic_paired`); `leave_one_out(paired=...)` removes the
  held-out model's (or family's) factorials; `own_flags_predictive` uses the estimated tau where present.
  `scripts/wiemip_paired_tests.py {synthetic,loo}` writes `results/wiemip/paired_tests_*.csv`. Full runs 2026-11-18 (16 replicates; 4 x 1000 draws; r-hat <= 1.008;
  gamma fits only: 5 of 160 synthetic and 7 of 84 leave-one-out fits have 1 to 3 divergences).
  - Synthetic, factorials scattered as observed: 90 % coverage for upstream is 0.19 (nitrogen) and 0.56
    (fire) on gamma_L, 0.25 and 0.75 on beta_L; for estimate 0.69, 0.94, 0.75, 0.94. With tau = 0 both
    are calibrated and estimate is sharper (posterior sd 0.8 against 1.6 PgC/K for fire on gamma_L).
  - Permafrost on gamma_L: coverage 0 in every treatment. The truth (-0.5 PgC/K) lies where the
    NegLogNormal prior has almost no mass, so no treatment of the factorials recovers it.
  - Leave one out: no treatment predicts held-out WIEMIP models clearly better than the raw ensemble (best
    crps_skill 0.06 for gamma_L, 0.14 for beta_L); with upstream, the all-processes estimate is biased by
    +22 PgC/K on gamma_L and covers 57 % at 90 %.
- **Next:** re-extract when the registry restores the missing models; a family column for the flag table.

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
   vegetation (`veg`), plus `family` for shared land-model code. As of 2026-10-22 no such table is known
   to exist. `data/wiemip/flags_template.csv` lists the registry models with the flags the factorial runs
   prove (1 where a model ran noX) and blanks to be filled by the modelling groups.
3. **Explore** with `notebooks/03_other_ensembles.ipynb`: set the paths at the top and run the cheap
   tests.
4. **Precompute:**
   `python scripts/precompute.py --name wiemip --values <values.csv> --flags <flags.csv> --quantities gamma_land_coubgc:additive beta_land_coubgc:multiplicative --jobs <n>`.
   - Use `--quick` first.
   - A full run took about 2.5 hours on 4 cores for two quantities and 13 models.
   - Commit `results/precomputed/wiemip/` (lightweight WIEMIP data may go in this private repository); otherwise commit only the
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
