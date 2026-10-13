# CLAUDE.md

Read `HANDOVER.md` first; it is the source of truth for scope, background and priorities.

Operating rules:
- UK English; never use em-dashes; neutral technical register. Describe methods, not people.
- `external/TCRE_Assessment` is a pinned submodule (commit 98c74cb). Never edit it; wrap or re-implement in `src/tcre_evidence/`.
- Logic lives in `src/tcre_evidence/`; notebooks are thin; every reported number comes from a script in `scripts/` with a fixed seed.
- `external/esm-deck` (commit eef04dc) is the primary data source for 1pctCO2 BGC/COU/RAD and flat10 global means. Its C4MIP ascii files use mixed flux sign conventions and formats by model: harmonise with an explicit per-model table and assert signs (HANDOVER §6.0).
- Name the β/γ convention in every variable (`gamma_coubgc`, `gamma_rad`, `eps`). Units: β PgC/ppm, γ PgC/K, α K/ppm, TCRE K/EgC. Upstream's `AF` is AF/k.
- Check r̂ < 1.01, ESS and divergences for every PyMC fit; record versions and sampler settings.
- Report per-model results and leave-one-family-out variants; the ensemble is ~11–13 models.
- Start with milestones M0–M1 in HANDOVER.md §8. If CMIP data access is slow, proceed with analyses that run on upstream inputs (A3–A6).
- Ask Ben before contacting anyone, pushing to a public remote, or writing conclusions about the upstream authors' work.
