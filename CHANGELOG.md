# Changelog

## 2026-10-05: accepted sign correction and public saved-product selection

- Select the accepted C23 main FF/FK products, corrected FK at five redshifts,
  native vertex inputs, and explicit tree/auto and BiHalofit auto/given cutoff
  series. The selected C callable now binds its table inside the package.
- Compact non-scientific session metadata in new public derivatives. Every
  scientific NPY member retains its original bytes; a separate lineage record
  distinguishes accepted-original and public-file hashes.
- Preserve the author-selected approximate diagonal panels and four historical
  reference PDFs while using corrected main products for the cross panel.
  Current calculated values and retained approximate quotes remain distinct.
- These products have finite-amplitude and sampled-sign acceptance. Global
  accuracy, new multiz FF folds, and a full numerical replay are not certified.

## 2026-10-03: vertex band builds take the P(k) table as an option

- `kappa3_vertex/rebuild/build_band.py --pk-table` chooses between
  `PCAMBz0.txt` (the default, used by every existing piece) and
  `PCAMB_pyccl_stf_fid_z0.txt` (the two-point side's table). The choice feeds
  the bispectrum model, `pk_matter_today` and the BiHalofit sigma8
  (`b_model.cosmology.fiducial_for_table`), and the build metadata records the
  table with its sha256. `merge_bands.py` and `assemble.py` refuse to combine
  bands built on different tables; an assembled table records `pk_table_low`
  and `pk_table_high`.
- The P(k) default is unchanged. BiHalofit leg ordering now defaults to `auto`;
  use `--leg-order given` for the historical ordering. The T-001 vertex rebuild
  uses the two-point power table, and mixed P(k) or leg-order bands are rejected.

## 2026-10-02: superseded products archived

- Moved 53 superseded files (17 MB) out of the working tree into the local,
  untracked `_archive/cleanup_2026-10-02/`: the old Order-0, FF and FK sweeps,
  the old cutoff-ladder and multi-redshift arrays, the August n_phi = 64 vertex
  tables, the 2026-09-06 folds at the historical source plane, superseded
  cross-check arrays, and the earlier Order-0 comparison with its two scripts.
  `reproduce/archived_inputs.json` lists each file with its SHA256; all 53 can
  be recovered from commit `ac7ff24`.
- The figure and number readers now require `reproduce/active_products.json`.
  Their historical fallback paths are gone. `check_numbers.py` reports 32
  summaries; the dropped one described the archived June FK sweep.
- The referee-revision drivers (`replay_suite.py`, `companion_suite.py`,
  `assemble_figure_products.py`) read archived inputs through the new
  `reproduce/archived_inputs.py`, which checks each archived copy against its
  recorded hash and keeps the recorded path, so existing run records still
  verify. A control that runs on an archived vertex table needs that table
  restored in place first.
- Kept in place: the run configurations (templates for replay) and the cut1000
  vertex table with its callable, which the accepted FF configurations import.
- 37 older analysis and probe scripts name an archived array and need it
  restored before they run. `_archive/cleanup_2026-10-02/dependents.json`
  lists them.

No scientific input or formula changed. The validation is recorded in
`reproduce/provenance/archive_20261002_validation.json`.

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
