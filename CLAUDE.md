# CLAUDE.md

Read `STATUS.md` first for the current state and next steps, then `HANDOVER.md` for the original scope and background.

Operating rules:
- UK English; never use em-dashes; neutral technical register. Describe methods, not people.
- `external/TCRE_Assessment` is a pinned submodule (commit 98c74cb). Never edit it; wrap or re-implement in `src/tcre_evidence/`.
- Reusable logic lives in `src/tcre_evidence/`; notebooks may show analysis code in their cells but read sampled results from caches written by scripts; every reported number comes from a script in `scripts/` with a fixed seed.
- `data/c4mip/` and `data/flat10/` (copied unchanged from esm-deck at commit eef04dc) are the primary data source for 1pctCO2 BGC/COU/RAD and flat10 global means. The C4MIP ascii files use mixed flux sign conventions and formats by model: harmonise with an explicit per-model table and assert signs (HANDOVER §6.0).
- Name the β/γ convention in every variable (`gamma_coubgc`, `gamma_rad`, `eps`). Units: β PgC/ppm, γ PgC/K, α K/ppm, TCRE K/EgC. Upstream's `AF` is AF/k.
- Check r̂ < 1.01, ESS and divergences for every PyMC fit; record versions and sampler settings.
- Report per-model results and leave-one-family-out variants; the ensemble is ~11–13 models.
- Milestones M0 to M11 are complete for CMIP6; the next step is the WIEMIP application (STATUS.md).
- Ask Ben before contacting anyone, pushing to a public remote, or writing conclusions about the upstream authors' work.
