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

| Script | Milestone | Output |
|---|---|---|
| `scripts/m0_reproduce_upstream.py` | M0 | `results/m0_upstream_reproduction.json` |
| `scripts/m1_covariance_check.py` | M1 | `results/m1_covariance.json`, `m1_alpha_alignment.csv`, `m1_input_sources.csv` |
