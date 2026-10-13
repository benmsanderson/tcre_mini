# tcre_mini: TCRE lines-of-evidence testbed

Diagnostics of the Bayesian process-based TCRE framework (`netzeroasap/TCRE_Assessment`) and a convention-free TCRE = TCR × R line of evidence. Scope, background and milestones are in `HANDOVER.md`; operating rules in `CLAUDE.md`.

## Setup

```
git clone --recurse-submodules https://github.com/benmsanderson/tcre_mini
conda env create -f environment.yml    # or: pip install -e ".[test]"
pytest                                  # fast tests; `pytest -m slow` runs MCMC tests
```

PyMC is pinned to 5.x and ArviZ to 0.x, matching the API used by the upstream notebook (`az.hdi(..., hdi_prob=)`).

## Layout

- `external/TCRE_Assessment`: upstream, pinned at 98c74cb (never edited).
- `external/esm-deck`: C4MIP 1pctCO2 BGC/COU/RAD and flat10 global means, pinned at eef04dc.
- `data/upstream_snapshot/`: frozen copy of the upstream Google Sheet and `DATA/` (see its README).
- `src/tcre_evidence/`: package. `upstream_patched/` is a minimally patched copy of the upstream builders (the pinned commit does not construct; the patch and reasons are documented there with an exact diff).
- `scripts/`: one script per analysis, writing to `results/`.
- `tests/`: pytest.

## Scripts

Results and their interpretation are collected in `results/SUMMARY.md`. Scripts that sample with PyMC pin BLAS to one thread; `TCRE_JOBS` sets the number of parallel fits where relevant.

| Script | Milestone | Main outputs (in `results/`) | Runtime (4 cores) |
|---|---|---|---|
| `scripts/fetch_cmip_globalmeans.py`, `python -m tcre_evidence.obs` | data | `data/derived/cmip_globalmean/`, `data/obs/` | hours (network) |
| `scripts/m0_reproduce_upstream.py` | M0 | `m0_upstream_reproduction.json` | ~20 min |
| `scripts/m1_covariance_check.py` | M1 | `m1_covariance.json`, `m1_alpha_alignment.csv`, `m1_input_sources.csv` | seconds |
| `scripts/m2_feedbacks_closure.py` | M2 / A1 | `m2_*`, `data/derived/c4mip_feedbacks_*.csv` | seconds |
| `scripts/m3_convention_attribution.py` | A2 | `m3_attribution.*`, `fig_m3_attribution.png`, `fig_m3_eps.png` | seconds |
| `scripts/m3_covariance_rerun.py`, `scripts/plot_m3_covariance.py` | A3 | `m3_covariance.json`, `m3_samples.npz`, `fig_m3_covariance.png` | ~40 min |
| `scripts/m4_perfect_model.py` | A4 | `a4/`, `m4_scores.csv`, `m4_summary.json`, `fig_m4_perfect_model.png` | ~1.5 h |
| `scripts/m5_tcr_r.py` | B | `m5_*`, `fig_m5_relationships.png` | ~2 min |
| `scripts/m6_prior_sensitivity.py`, `scripts/m6_shift_decomposition.py` | A5 | `m6_prior_*`, `m6_shift_decomposition.json`, `fig_m6_prior_posterior.png` | ~1 h |
| `scripts/m6_identifiability.py` | A6 | `m6_identifiability.csv`, `m6_nu_evidence.json`, `fig_m6_identifiability.png` | ~2 h (1 job) |
| `scripts/m6_state_dependence.py` | A7 | `m6_state_dependence.*`, `fig_m6_state_dependence.png` | seconds |
| `scripts/m6_comparison.py` | comparison | `fig_m6_comparison.png`, `m6_comparison.csv` | seconds |
