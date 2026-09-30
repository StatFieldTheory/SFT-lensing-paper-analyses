# Changelog

## 2026-09-30: reproduction package cleanup

- Rewrote the README, reproduction guide, figure map and calculation overview for
  readers of the revised paper. The documentation identifies its AI-assisted
  authorship and distinguishes saved-product regeneration from numerical reruns.
- Preserved the active R1 products, their recorded hashes, numerical code,
  configurations, runtime pins and symbolic checks. Historical inputs that remain
  dependencies of current entry points are retained.
- Archived old reorganization records, independent August audit and revision
  drivers, the superseded PSD investigation, the unused shear-3PCF backup, and
  generated logs. These files remain locally under
  `_archive/cleanup_2026-09-30/`, with a hash-verified manifest. The 548 previously
  tracked files can also be recovered from commit `866ea12`.
- Stopped tracking the root `docs/` tree and agent plans, handoffs and run-state
  notes (74 files). They remain at their original local paths and are ignored by
  Git. Numerical product records and user-facing documentation remain tracked.
- Preserved the vertex-rebuild provenance note byte-for-byte at
  `reproduce/provenance/vertex_rebuild_20260906.md` and updated the replay driver
  to use that tracked copy. Historical run manifests retain their original paths.
- Made the number-check entry point locate its own package and honor
  `CANOES_ROOT` for the external power-spectrum table.

The cleanup does not change scientific inputs or numerical formulas. Reproduction
validation and its scope are recorded in
`reproduce/provenance/cleanup_20260930_validation.json`.

## 2026-09-29: referee revision products

Commit `866ea12` records the R1 numerical products and reproducibility evidence:
the actual source-plane endpoints, SFT-Wick 0.6.1, corrected pair-phase inputs,
the permutation-aware cutoff ladder, and the aligned Order-0/PyCCL comparison.
`reproduce/active_products.json` selects the accepted figure inputs.

## Earlier history

The detailed September reorganization changelog is retained in
[CHANGELOG_REORG_2026-09.md](CHANGELOG_REORG_2026-09.md). Its historical paths and
validation statements describe that earlier state; some supporting records are
now in the local cleanup archive or Git history.
