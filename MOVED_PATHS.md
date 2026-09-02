# MOVED_PATHS.md

Old path to new path for the 2026-09-02 reproducibility reorganisation of
`~/Documents/MyDrafts/STF_lensing`. Written for the sessions of `SFT-WL-B` and `rezeta`,
which vendor files from this tree and pin them by SHA-256, and for anything else holding
a path into it.

All paths are relative to `~/Documents/MyDrafts/STF_lensing/`. The manuscript itself did
not move: `main.tex`, `sections/`, `figures/` and `biblio.bib` are where they were.

**Two things changed shape.** The `driver_field_emulators/` folder no longer exists: it
was the 2026-08 working tree and its contents are now inside the analysis package, which
is a git repository (`SFT-lensing-paper-analyses`, GitHub `StatFieldTheory/SFT-lensing-paper-analyses`).
And `SFT-lensing-paper-analyses/_archive/` no longer exists: the 2026-05 canoes pipeline
and the 2026-06 cleanup archive were retired to the Trash, their contents listed with md5
in `reorg_2026-09/trash_manifest_2026-09-02.txt`. Everything they still supplied is
vendored in the live tree (see the last section).

## 1. Files vendored by SFT-WL-B

Checked 2026-09-02 at the new paths. "unchanged" means the bytes are the same as the
copy vendored on 2026-08-31, so the pinned SHA-256 still matches and only the path in
`vendor/*/PROVENANCE.md` needs editing.

| vendored as | old path | new path | bytes |
|---|---|---|---|
| `stf_lensing_mirror/.../callable_fixed/perm_aware_kappa3_callable.py` | `driver_field_emulators/code/callable_fixed/perm_aware_kappa3_callable.py` | `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable.py` | unchanged (`4d7fefaa129e1f1e`) |
| `.../callable_fixed/README.md` | `driver_field_emulators/code/callable_fixed/README.md` | same folder as above, `README.md` | unchanged (`a0d0001a84f39c3c`) |
| `.../callable_fixed/table_permclosed.npz` | `driver_field_emulators/products/table_permclosed_cut15360.npz` | same folder as above, `table_permclosed.npz` | unchanged (`b400dc5b11e1033b`); it now carries the name the callable resolves, so the mirror's rename step is no longer needed |
| `.../callable_fixed/compare_callables.py` | `driver_field_emulators/code/callable_fixed/compare_callables.py` | same folder as above | **changed** (`73ac4f587532924d` to `b9f7989e0c0b7754`): the hardcoded Order-0 path became a relative one |
| `.../callable_fixed/tests/test_perm_aware.py` | `driver_field_emulators/code/callable_fixed/tests/test_perm_aware.py` | same folder as above, `tests/` | **changed** (`1eaf330b5712081b` to `d0cbaaa51ee05a89`): the `parents[...]` index that locates the deployed table |
| `.../equal_time_limber_kappa3_callable.py` | unchanged path under `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber/` | unchanged | unchanged (`4b9c47cc8771ab31`) |
| `.../equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz` | unchanged path | unchanged | unchanged (`76342b35c5637465`) |
| `analysis3_cl/.../xi_C_corr_op_O0.npz` | unchanged path | unchanged | unchanged (`517d19aa55c36aab`) |
| `analysis3_cl/.../xi_C_corr_op_K_limber_FF.npz` | unchanged path | unchanged | unchanged (`cc3df79f8ea12427`) |
| `analysis3_cl/.../xi_C_corr_op_K_limber_FK.npz` | unchanged path | unchanged | unchanged (`c3eae8cbe8515668`) |
| `analysis3_cl/.../plot_cl_EB_polarization.py` | unchanged path | unchanged | unchanged (`7d92bc6b449ad9f1`) |
| `analysis3_cl/.../_plot_style.py` | unchanged path | unchanged | unchanged (`b9e7b0cc6bbd9d5c`) |
| `analysis3_cl/.../plot_analysis3_cl_decomposition.py` | unchanged path | unchanged | **changed** (`455dd479c0278147` to `f28f9f5d21424bc9`): see the warning below |
| `cutoff_study/cutoff_vertex.npz` | `driver_field_emulators/products/cutoff_vertex.npz` | `SFT-lensing-paper-analyses/sachs_sft/analyses/fk_audit_2026-08/products/cutoff_vertex.npz` | unchanged (`c429cb9b88807e49`); untracked by git, local only |
| `consumer_pk/PCAMBz0.txt` | `~/projects/angular_statistics/canoes/examples/data/PCAMBz0.txt` | unchanged (a canoes file, not ours) | unchanged |

**Warning about `plot_analysis3_cl_decomposition.py`.** Its `FK_NPZ` default now points at
the manuscript's converged sweep instead of the superseded June one:

```
old: sftwick_outputs/2PCF/C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK.npz
new: sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz
```

A vendored copy that reads the module-level constant will change which FK it loads if it
is re-vendored. The vendored *data* files are unaffected, and the B/E numbers a gate
re-derives from the cut15360 sweep are the manuscript's (`0.95` at `ell = 60` falling to
`0.42` at `ell = 1500`). The June sweep is still present, with a `SUPERSEDED` banner in
its README, so a copy that names it explicitly keeps working.

## 2. Paths cited in `SFT-WL-B/docs/decision_log_bmode.md`

| `[P]` reference | status |
|---|---|
| `sections/conclusion.tex`, `sections/insights.tex`, `sections/path_int.tex`, `sections/appendix.tex` (line numbers) | unchanged; the manuscript was not edited |
| `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber/equal_time_limber_kappa3_callable.py:155-157` (the `np.sort` then `np.round(..., 8)` lookup) | unchanged file and line numbers. The decision log's note that it is "still unfixed in production" is now out of date: that callable is no longer production. The production vertex is `equal_time_limber_cut15360_permaware/`, and `equal_time_limber/DEPRECATED.md` says so beside the old file. |
| `driver_field_emulators/notes/finding_grid_artifact_measured.md` | `SFT-lensing-paper-analyses/docs/fk_audit_2026-08/notes/finding_grid_artifact_measured.md` |
| `driver_field_emulators/notes/emulator_candidates.md` | `.../docs/fk_audit_2026-08/notes/emulator_candidates.md` |
| `driver_field_emulators/notes/input_ready_b_to_zeta_spec.md` | `.../docs/fk_audit_2026-08/notes/input_ready_b_to_zeta_spec.md` |
| `driver_field_emulators/note` (the LaTeX audit note) | `.../docs/fk_audit_2026-08/note/` |
| `driver_field_emulators/code/rebuild/assemble.py` and the rest of `code/rebuild/` | `.../sachs_sft/callables/kappa3_vertex/rebuild/` |
| `driver_field_emulators/code/figures_corrected/regenerate.py` | `.../sachs_sft/analyses/revision_2026-08/regenerate.py` |
| `driver_field_emulators/products/` | split by role: paper products into the run folders under `sachs_sft/sftwick_outputs/2PCF/`, rebuild inputs into `.../rebuild/products/`, audit data into `.../analyses/fk_audit_2026-08/products/` |
| `SFT-lensing-paper-analyses/figure_manifest.md` | unchanged path, now tracked by git |
| `SFT-lensing-paper-analyses/REPRODUCE.md` (referenced prospectively) | now exists |
| `SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt/fk_analytic.py` | unchanged |
| `talk/assets/figures/_data/multiz_components.npz` | unchanged in the talk repository; a byte-identical copy with provenance is at `.../sachs_sft/analyses/analysis3/inputs/multiz_components_talk_2026-06-10.npz` |
| `figures/archived_2026-08-*`, `figures/_archived/` | removed from the working tree; recoverable from the paper repository's history at tag `pre-reorg-2026-09` |

Every note under `driver_field_emulators/notes/` moved to
`SFT-lensing-paper-analyses/docs/fk_audit_2026-08/notes/` keeping its file name. Thirteen
of them gained a dated correction block saying that their absolute FK amplitudes are
`ell_max = 1000` values; the key is `docs/fk_audit_2026-08/FK_BASELINE_NUMBERS.md`.

## 3. Paths cited by rezeta

`rezeta` references this tree only by its root (`~/Documents/MyDrafts/STF_lensing/`, read
only) and by two files:

| old | new |
|---|---|
| `sections/appendix.tex` | unchanged |
| `driver_field_emulators/notes/input_ready_b_to_zeta_spec.md` | `SFT-lensing-paper-analyses/docs/fk_audit_2026-08/notes/input_ready_b_to_zeta_spec.md` |

## 4. Everything else that moved

| old | new |
|---|---|
| `driver_field_emulators/code/rebuild/` | `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/rebuild/` |
| `driver_field_emulators/code/b_model/` | `.../sachs_sft/callables/kappa3_vertex/b_model/` |
| `driver_field_emulators/code/callable_fixed/` | `.../sachs_sft/callables/kappa3_vertex/equal_time_limber_cut15360_permaware/` |
| `driver_field_emulators/code/mc_fk_complete/` | `.../sachs_sft/analyses/mc_fk_complete/` |
| `driver_field_emulators/code/figures_corrected/` | `.../sachs_sft/analyses/revision_2026-08/` |
| `driver_field_emulators/code/{probe_*,redshift_*,sweep_ell_high_max,plot_fk_variants*,compare_fk_variants,mc_crosscheck_dense,verify_*.wl}` | `.../sachs_sft/analyses/fk_audit_2026-08/` |
| `driver_field_emulators/code/{normalisation,reduced_shear}/` | `.../sachs_sft/analyses/fk_audit_2026-08/{normalisation,reduced_shear}/` |
| `driver_field_emulators/code/{run_fk_variant,make_dense_triples,make_collapsed_triples}.py` | `.../sachs_sft/callables/kappa3_vertex/rebuild/` |
| `driver_field_emulators/code/README.md` | `.../sachs_sft/analyses/fk_audit_2026-08/README_code_2026-08.md` |
| `driver_field_emulators/products/table_permclosed_cut15360_permfix_xi.npz` | `.../sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz` |
| `driver_field_emulators/products/table_{tree,bihalofit}_cut*_r4_xi.npz` | `.../sachs_sft/sftwick_outputs/2PCF/cutoff_ladder/` |
| `driver_field_emulators/products/table_{tree,bihalofit}_cut*_r4.npz` | `.../sachs_sft/callables/kappa3_vertex/rebuild/products/` |
| `driver_field_emulators/products/{multiz_ff_real*,order0_multiz_xi,source_distances}.npz` | `.../sachs_sft/sftwick_outputs/2PCF/multiz/` |
| `driver_field_emulators/products/zeta_bands_cut15360*.npz` | `.../sachs_sft/callables/kappa3_vertex/equal_time_limber/outputs/` |
| `driver_field_emulators/products/{triples_*,fk_kernel_crosscheck,mean_kappa_plateaus}.npz`, `logs/` | `.../sachs_sft/callables/kappa3_vertex/rebuild/{products,logs}/` |
| `driver_field_emulators/products/pieces/` | `.../sachs_sft/callables/kappa3_vertex/rebuild/products/pieces/` (untracked) |
| `driver_field_emulators/{README,SUMMARY,HANDOVER,SESSION_STATE_2026-08-26,CHANGES_OUTSIDE_THIS_FOLDER}.md` | `.../docs/fk_audit_2026-08/` (the README as `README_driver_field_emulators_2026-08.md`) |
| `driver_field_emulators/PROMPT_*.md`, `sftwick_outputs/PROMPT_*.md`, `shear3pcf_fastnc/PROMPT*.md` | `.../docs/agent_briefs/` |
| `driver_field_emulators/papers/` | `.../docs/papers/arxiv_2026-08/` (untracked) |
| `driver_field_emulators/revision_plan/` | `.../docs/private/revision_plan_2026-08/` (untracked) |
| `driver_field_emulators/figures/fk_variant*.{pdf,png}` | `.../sachs_sft/analyses/fk_audit_2026-08/figures/` |
| `driver_field_emulators/figures/corrected/` | removed; the generators now write to their own `outputs/` folders |
| `talk/scripts/run_multiz_components.py` | code at `.../sachs_sft/analyses/analysis3/multiz_sweep.py`; the talk file is now a wrapper that re-exports it |
| `SFT-lensing-paper-analyses/_archive/canoes_pipeline/scripts/plot_coord_maps.py` | `.../figures/make_cosmology_coord_maps.py` (canoes dependency vendored away) |
| `SFT-lensing-paper-analyses/_archive/canoes_pipeline/scripts/config_FK_LIMBER_z5.yaml` | `.../sachs_sft/scripts/inputs/config_FK_LIMBER_z5_archived_2026-05.yaml` |
| the L2 lambda grid inside `_archive/canoes_pipeline/inputs/l2_callables_2026_05_26/` | `.../sachs_sft/callables/kappa3_vertex/equal_time_limber/inputs/l2_lambda_grid_z5covgrid16.npz` (the one array the builds read) |
| the Feynman diagrams' generator (`sft-wick/examples/nonlocal_vertex_2pt.ipynb`, cells 29 and 38) | extracted to `.../figures/make_feyndiag.py`; the notebook stays in sft-wick |
| root `analysis_setup.md`, `findings.md`, `progress.md`, `task_plan.md`, `outputs/` | retired to the Trash folder, listed in `reorg_2026-09/trash_manifest_2026-09-02.txt` |
| `SFT-lensing-paper-analyses/_archive/` (canoes pipeline, 2026-06-17 cleanup, 2026-06-09 archive) | retired to the Trash folder, same manifest |

## 5. What was retired rather than moved

`~/.Trash/STF_lensing_reorg_2026-09-02/` holds, until the Trash is emptied: the 2026-05
canoes pipeline archive, the 2026-06-17 cleanup archive, the 2026-06-09 h^4 archive, the
May 2026 figure archive, the obsolete root planning documents, the LaTeX build residue,
and the untracked remainder of the imported `driver_field_emulators` copy. Every file is
listed with its md5 and original path in `reorg_2026-09/trash_manifest_2026-09-02.txt`
and `reorg_2026-09/dfe_manifest_2026-09-02.txt`. Anything tracked by git before the move
is also recoverable from the tag `pre-reorg-2026-09` in the repository that held it.
