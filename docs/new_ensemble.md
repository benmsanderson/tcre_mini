# Applying the process tests to a new ensemble (for example WIEMIP)

This walks through running `scripts/process_tests.py` on an ensemble other than CMIP6. You need two
CSV files and a few minutes for the cheap tests; the Bayesian tests take longer. Read `docs/guide.md`
first if the process model is new to you.

## 1. The values table

One row per model, a `model` column, and one column per quantity. `data/examples/cmip6_land_feedbacks.csv`
is the CMIP6 version:

```
model,beta_land_coubgc,gamma_land_coubgc,gamma_land_tropics_coubgc,gamma_land_extratropics_coubgc
ACCESS-ESM1-5,0.7042,-10.1079,-14.3683,4.2604
...
```

- **Units.** β in PgC/ppm, γ in PgC/K.
- **Convention in the name.** Say which convention each column uses. With global means of the 1pctCO2
  runs (COU fully coupled, BGC biogeochemically coupled, RAD radiatively coupled), averaged over a window
  (years 61 to 80 for 2×CO₂) and taken as anomalies from piControl:

      gamma_coubgc = (ΔC_COU − ΔC_BGC) / (ΔT_COU − ΔT_BGC)
      beta_coubgc  = (ΔC_BGC − gamma_coubgc · ΔT_BGC) / ΔCa
      gamma_rad    = ΔC_RAD / ΔT_RAD

  `tcre_evidence.feedbacks.diagnose()` does this for the C4MIP ascii files in `data/c4mip/`, including the
  sign checks; for other formats, compute the same quantities from your own global means.
  COU−BGC matches the development code. RAD differs from it by the non-linear term ε, so do not mix
  the two in one column.
- **End-to-end quantities.** If you can, also include TCRE/TCR (the carbon factor R) or the airborne
  fraction. Process effects on β and γ partly cancel, and the end-to-end quantity is what matters for
  TCRE.

## 2. The flags table

One row per model, a `model` column, one 0/1 column per process, and optionally a `family` column.
`data/examples/cmip6_process_flags.csv` is the CMIP6 version:

```
model,nitrogen,permafrost,fire,veg,family
CESM2,1,1,1,0,NCAR
...
```

- **Add the processes you want to test.** For example, separate `demography` (cohort or individual
  based vegetation) from `veg` (dynamic biogeography). A process the development code does not have
  gets a wide prior centred on no effect.
- **`family`** groups models that share code (land model, ocean model). It is used by the
  leave-one-family-out test (`--groups`). It also lets you check whether a flag only appears within one
  family. In CMIP6 the permafrost flag was set only for the three CLM-based models, so it cannot be told
  apart from a CLM offset.
- **Avoid flags no model differs on.** A flag set for every model, or for none, cannot be estimated, and
  `load_ensemble` warns about it.
- **Calibration information.** Worth recording in separate columns if groups will share it (for example
  `face_calibrated`, `recalibrated_after_historical`). These are not processes, but they can be tested
  in the same way, and they matter for interpreting a process effect.

## 3. Run the cheap tests

```bash
python scripts/process_tests.py --values my_values.csv --flags my_flags.csv \
    --quantity gamma_land_coubgc --kind additive --tests budget power \
    --effects demography=-20 nitrogen=10
```

- `--kind additive` for γ (effects δ in PgC/K); `--kind multiplicative` for β (effects are factors η,
  e.g. `nitrogen=0.73`).
- `budget` reports how much of the spread the flags explain. A negative leave-one-out share means the
  flags predict new models worse than the ensemble mean.
- `power` reports the standard error of each effect and how many independent models would be needed to
  detect the size given in `--effects`. Without `--effects`, it uses the regression's own estimate.

If `power` says a few hundred models are needed, the Bayesian fit cannot learn that effect from this
ensemble. Its posterior will then follow its prior, and the next step will show that.

## 4. Run the Bayesian tests

```bash
python scripts/process_tests.py --values my_values.csv --flags my_flags.csv \
    --quantity gamma_land_coubgc --kind additive --tests prior synthetic loo \
    --effects demography=-20 --groups --jobs 4
```

- `prior` (one fit): compares each term's prior and posterior. `width_ratio` near 1 means not informed;
  a narrower posterior can still mean only that the data rule out the prior's tail.
- `synthetic` (`--nrep` fits, default 12): refits pseudo-ensembles with the effects in `--effects` and
  the same scatter as your data. Are the effects recovered?
- `loo` (one fit per model): does the real-world estimate, or the model's own-flag prediction, beat the
  raw ensemble for held-out models?
- `--quick` uses short chains for a first look. Report numbers from full runs, and check that `rhat_max`
  is below 1.01 in the output files.

Outputs go to `results/process_tests/<quantity>/` (or `--out`), with `run.json` recording the settings,
priors and package versions.

## 5. From Python

```python
from tcre_evidence import process_tests as pt

y, flags, families = pt.load_ensemble("my_values.csv", "my_flags.csv", "gamma_land_coubgc")
print(pt.flag_budget(y, flags))

priors = pt.default_priors("additive", list(flags.columns))
priors["delta"]["demography"] = dict(dist="Normal", mu=-10.0, sigma=10.0)   # your own prior
table, diagnostics = pt.prior_vs_posterior(y, flags, priors)

# paired experiments (same model, process on minus off) as direct evidence on delta
model = pt.build(y, flags, priors, paired={"demography": [-15.0]}, paired_sigma=5.0)
idata, diagnostics = pt.fit(model)
```

Prior specifications follow the development code (`external/TCRE_Assessment/config/priors.py`):
- **delta:** `Normal`, `LogNormal` or `NegLogNormal`, with `mu` and `sigma`;
- **eta:** `mu` and `sigma` of log η;
- **process scatter:** `sigma_process_struct` (additive) or `sigma_mult_struct` (multiplicative), which
  sets how much a process's effect varies between models.

## 6. Reading the results together

A process effect is supported by the ensemble when:
1. `power` says the ensemble is large enough for an effect of that size;
2. `synthetic` recovers a known effect of that size;
3. `prior` shows the posterior moved away from the prior;
4. the flagged and unflagged models do not simply sit in different families;
5. where possible, a paired experiment agrees.

When these fail, the effect in the result comes from the prior and should be reported as a judgement.
Keep in mind that calibration can hide a real effect (see `scripts/m9_tuning_selection.py`), so a
failed test does not show the process is unimportant, only that this ensemble cannot measure it.
