# REORG_PLAN.md: STF_lensing reproducibility reorganisation

**Phase 0 deliverable (inventory and freeze). Date 2026-09-02. Status: awaiting operator sign-off.**

Nothing has been moved, committed, deleted or edited in the manuscript. Every
claim below that says "verified" was verified by execution during this session;
the machine-readable evidence is in
`SFT-lensing-paper-analyses/reorg_2026-09/phase0/`
(`figure_checks.md|json`, `numbers.md|json`, `harness/`, `logs/`). The harness
scripts there (`p0_regen.py`, `pdfcmp.py`, `p0_numbers.py`) are the seed of the
Phase-2/5 reproduction checker.

Goal restated: a reader with the latest manuscript plus this repository must be
able to reproduce every figure and every quoted number without oral tradition.
Constraints honoured: manuscript read-only; archive rather than delete; Overleaf
layout untouched; nothing moved before a baseline commit; vendored files kept
byte-identical; no canoes build fan-out; no re-computation longer than an hour.

## 0. Corrections to the 2026-09-02 external survey (the brief)

The brief was checked item by item. These points differ from it and change the plan.

1. **The two Feynman diagrams do have a generator.** `figures/FeynDiag_2pt_order2.pdf`
   and `figures/FeynDiag_3pt_order1.pdf` are byte-identical (md5 `c26347c4...`,
   `66a713ae...`) to `/Users/zzhang/projects/SFT/sft-wick/examples/FeynDiag_*.pdf`,
   which are written by cells 29 and 38 of
   `sft-wick/examples/nonlocal_vertex_2pt.ipynb`
   (`result.draw_diagrams(order=2, ...)` and `result_3pt.draw_diagrams(order=1, ...)`),
   sft-wick commit `dd03748` (2026-06-03). They are not hand-made vector art. They
   belong to the sft-wick repository, so the plan pins them rather than moves them.
2. **`driver_field_emulators/` feeds 7 of the 17 figures, not 5**: rows 2, 3, 4, 6, 7
   through the rebinding drivers, row 5 through `regenerate_multiz.py` plus the real FF
   sweeps, row 17 directly. The parent `CLAUDE.md` is stale on this point.
3. **The multi-z compute step does not read the talk npz**; it imports code
   (`build_lambda_of_z`, `run_vertex_sweep`, `CONFIGS`) from
   `talk/scripts/run_multiz_components.py`. But the deployed figure's Order-0 planes
   ARE PCHIP-in-z interpolants of the talk cache `talk/assets/figures/_data/multiz_components.npz`
   (verified: maximum relative difference 0.0 against the June npz that the August
   regeneration reused with `--reuse`). So the talk data is a provenance dependency
   after all, through a different route than the manifest states.
4. **The FK 2PCF sweep with the permutation-aware `cut15360` table already exists**
   (`driver_field_emulators/products/table_permclosed_cut15360_permfix_xi.npz`,
   2026-08-26 09:48, fold time about 2.5 min) and IS the file the deployed FK figures
   read. Phase 3 therefore becomes "promote, then verify by one re-run", not "compute".
5. **Repository sync state**: the paper repo is 1 commit behind Overleaf
   (`1f0ec58 Update on Overleaf.`, which flips `\revhltrue` to `\revhlfalse` in
   `main.tex`); the analyses repo is 1 commit ahead of GitHub (`e27df7c`, unpushed);
   the talk repo is in sync. Phase 1 must fast-forward the paper repo first.
6. **Two external runtime dependencies are in dirty states**: sft-wick
   (`/Users/zzhang/projects/SFT/sft-wick`, HEAD `94cc5af` 2026-08-05) has uncommitted
   changes in `src/sft_wick/evaluate.py`, `workflow/config.py`, `workflow/cli.py`;
   canoes (`/Users/zzhang/projects/angular_statistics/canoes`, branch
   `feat/nonlinear-bdelta-vertex`, HEAD `d7b20ce`) has an uncommitted
   `src/canoes/sachs/kappa3.py` (the `b_delta_fn` patch). Any Phase-3 re-sweep runs on
   these working trees, so their state must be pinned as patch files before any run.
7. **The 18 uncommitted items** are 16 modified files plus 2 untracked files
   (`mathematica/beam_width_regulator.wl`, `_smoke_equal_time_limber.npz`).
8. **`figures/` holds two tracked orphans** not referenced by any live `.tex`:
   `appendix_3cumulant_fastnc.pdf` (referenced only by `sections/_archive/input_validation_fastnc.tex`)
   and `cumulant_hierarchy.pdf`. Both are on Overleaf.
9. **A local `block-no-verify` git hook false-positives on shell commands whose text
   contains "commit" and "--"**; Phase 1 commit messages must be single-line and plain,
   and long documents must be written with the file tool, not shell heredocs.

## 1. Figure provenance (17 figures, all exercised)

Environments: `pyccl` = `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`
(Python 3.11.15, matplotlib 3.10.8, pyccl 3.3.1, camb 1.6.5, scipy 1.17.1);
`sftw` = `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`
(Python 3.14.3, matplotlib 3.10.8, sft_wick from `~/projects/SFT/sft-wick/src`);
`+canoes` = `PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src`.
Verdicts: **IDENTICAL** = regenerated PDF equals the deployed one byte for byte
except the `/CreationDate` string (matplotlib writes nothing else non-deterministic);
**EQUIVALENT** = curves and text overlay, raster identical at 60 dpi, but the tight
bounding box differs by 0.2 to 1.0 pt (font-metric drift); **DIFFERS** = content differs.

| # | figure | generator (repo-relative) | key inputs | env that reproduces | reproducible today without manual rebinding? | result of the check |
|---|---|---|---|---|---|---|
| 1 | `analysis1_O0_vs_pyccl.pdf` | `SFT-lensing-paper-analyses/sachs_sft/analyses/analysis1/plot_analysis1_O0_vs_pyccl_fkem.py` | `sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz` (2026-06-01), `analysis1/pyccl_xi_shear_reference_fkem.npz` | `sftw` | yes | IDENTICAL in `sftw`; in `pyccl` the page height differs by 1 pt (EQUIVALENT). Deployed file is a raw copy of `analysis1/outputs/analysis1_O0_vs_pyccl_fkem.pdf`. |
| 2 | `analysis3_NLO_FFFK.pdf` | `analysis3/plot_analysis3_nlo_decomposition.py` | O0 and FF sweeps as above; FK = `driver_field_emulators/products/table_permclosed_cut15360_permfix_xi.npz` | `pyccl` | **no**: default `FK_NPZ` is the superseded `C_corr_op_K_limber_FK` sweep (cut1000, old callable); FK(0.5') comes out 4.54x low (4.288e-6 vs 1.948e-5) | standalone DIFFERS (2.3% of pixels); with `FK_NPZ` rebound (what `regenerate.py` does) IDENTICAL. Deployed = raw copy of `driver_field_emulators/figures/corrected/analysis3_nlo_O0_FF_FK_corrected.pdf`. The tracked in-tree `analysis3/outputs/analysis3_nlo_O0_FF_FK.pdf` is the STALE June render (16% of pixels differ). |
| 3 | `analysis3_cl_full_vs_O0.pdf` | `analysis3/plot_analysis3_cl_decomposition.py` (needs pyccl only for a console validation print) | same three sweeps | `pyccl` | **no** (same trap) | standalone DIFFERS (1.6%); rebound IDENTICAL. Tracked `analysis3/outputs/analysis3_cl_O0_FF_FK.pdf` is STALE (June). |
| 4 | `cl_EB_polarization.pdf` | `analysis3/plot_cl_EB_polarization.py` (imports the row-3 module for its transform and its `FK_NPZ`) | FF and FK sweeps | `pyccl` | **no** (same trap) | standalone DIFFERS (1.1%); rebound IDENTICAL. Tracked `analysis3/outputs/cl_EB_polarization.pdf` is STALE (June). |
| 5 | `multiz_kappa_xi_cl.pdf` | plot stage `analysis3/plot_multiz_kappa_2x5.py`; compute stage `analysis3/compute_multiz_kappa_2pcf.py` driven by `driver_field_emulators/code/figures_corrected/regenerate_multiz.py --reuse` | `analysis3/outputs/multiz_kappa_2pcf_5z.npz` (uncommitted, 2026-08-26): O0 planes = PCHIP-in-z of `talk/assets/figures/_data/multiz_components.npz`; FK planes = perm-aware cut15360 fold at 5 planes (1034 s); FF planes from `driver_field_emulators/products/multiz_ff_real_all5.npz` (five single-threaded sweeps, about 3.1 h each, 2026-08-27) | `pyccl` (plot) | plot stage: yes, from the cached intermediates; compute stage: not without the talk repo (imports `talk/scripts/run_multiz_components.py` by absolute path) and not in under an hour | plot stage IDENTICAL. Compute stage not re-run in Phase 0. Deployed = raw copy of `analysis3/outputs/multiz_kappa_xi_cl_2x5.pdf`. |
| 6 | `appendix_mc_workflow.pdf` | `mc_sachs_2pt/fig_xi_channels.py`, deployed through `driver_field_emulators/code/figures_corrected/make_val_figure.py --fk-xi <cut15360 permfix> --fk-markers code/mc_fk_complete/_markers_pooled.npz` | O0/FF sweeps; FK cut15360 fold; `mc_sachs_2pt/outputs/appendix_mc_curve.npz` (FF Monte-Carlo cache, seed 11); `_markers_pooled.npz` (4 blocks x 24 seeds, 2026-08-28) | `pyccl +canoes` (`driver_stats.py` imports the corr_op callable, which imports `canoes.sachs.sft_input.corr_op.table`) | only because the cache already holds the corrected FK line and the FK markers: `--from-cache` reproduces; a fresh `_compute` would read the cut1000 FK and emit no FK markers | IDENTICAL by both routes (`--from-cache`, and the `make_val_figure.py` logic). Deployed = raw copy of `mc_sachs_2pt/figures/xi_kappa_channels.pdf`. |
| 7 | `zeta_driving_field_slices_draft.pdf` | `callables/kappa3_vertex/equal_time_limber/plot_zeta_figures.py`, deployed through `driver_field_emulators/code/figures_corrected/regenerate_zeta_slices.py --plot-only --gamma-max 85` | `driver_field_emulators/products/zeta_bands_cut15360_gmax85.npz` (built by `ell_band_decomp.py` with eight HIGH windows to 15360, 2026-08-27) | `pyccl` | **no**: default `_NPZ` is the June cut1000 `ell_band_decomp_results.npz`, and `_FIGDIR` is the repo `figures/` by absolute path (a bare run overwrites the paper file) | standalone DIFFERS (panel geometry differs); rebound IDENTICAL. Deployed = raw copy of `driver_field_emulators/figures/corrected/zeta_driving_field_slices_draft.pdf`. |
| 8 | `corr_operator_slices_draft.pdf` | `SFT-lensing-paper-analyses/figures/make_response_corr_operator_figures.py` | `callables/C_propagator/corr_op/corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz` (5.5 MB), `sachs_sft/scripts/D_callable.py` | `pyccl +canoes` | yes (writes straight into `figures/`) | EQUIVALENT: page 661.24 x 601.60 pt vs 661.45 x 601.82 pt, 0.8% of pixels at glyph and line edges, curves overlay. Not byte-reproducible today. |
| 9 | `response_operator_draft.pdf` | same script as row 8 | same | `pyccl +canoes` | yes | EQUIVALENT: 407.76 x 256.55 vs 408.00 x 256.79 pt, 2.1% of pixels at edges. |
| 10 | `screen_basis_spin2_pattern.pdf` | `SFT-lensing-paper-analyses/figures/make_screen_basis_spin2_figure.py` (`render()`) | none (procedural) | either | yes | EQUIVALENT in both envs: 509.69 x 275.22 vs 509.97 x 275.67 pt, 1.8% of pixels; embedded fonts identical (Cmr10/Cmsy10/Cmmi10). The script requests "CMU Serif", absent in both envs today; the fallback is consistent. |
| 11 | `driving_field_portrait.pdf` | `SFT-lensing-paper-analyses/figures/plot_driving_field_portrait.py --seed 42` | CAMB linear P(k) computed at run time | `pyccl` (camb 1.6.5) | yes | IDENTICAL. |
| 12 | `selection_rule_schematic.pdf` | `SFT-lensing-paper-analyses/figures/make_selection_rule_schematic.py --seed 42` | none | `pyccl` | yes | IDENTICAL in `pyccl`; `sftw` gives 0.14% edge pixels. |
| 13 | `null_geodesic_congruence.pdf` | `SFT-lensing-paper-analyses/figures/make_null_congruence_schematic.py` | none | `pyccl` | yes | IDENTICAL in `pyccl`; `sftw` gives 0.09%. |
| 14 | `cosmology_coord_maps.pdf` | `SFT-lensing-paper-analyses/_archive/canoes_pipeline/scripts/plot_coord_maps.py` (archive only) | `_archive/canoes_pipeline/cosmo.py` (`FiducialCosmology`) which imports constants from `canoes.cosmo.background` | `pyccl +canoes` plus a package shim `scripts/canoes_pipeline -> _archive/canoes_pipeline` (the script still imports the pre-2026-06-17 package name) | **no** without the shim and canoes | IDENTICAL with the shim in `pyccl`; in the canoes `.venv` (matplotlib 3.10.9) raster-identical but the producer string differs. |
| 15 | `FeynDiag_2pt_order2.pdf` | `~/projects/SFT/sft-wick/examples/nonlocal_vertex_2pt.ipynb` cell 29 | sft-wick `Action`/`compute_moment` diagram set, drawn by `DiagramRenderer` | `sftw` (notebook) | **no**: the generator is in another repository and is a notebook | byte-IDENTICAL to `sft-wick/examples/FeynDiag_2pt_order2.pdf`; notebook not executed in Phase 0. |
| 16 | `FeynDiag_3pt_order1.pdf` | same notebook, cell 38 | same | `sftw` | **no** (same) | byte-IDENTICAL to `sft-wick/examples/FeynDiag_3pt_order1.pdf`. |
| 17 | `fk_cutoff_convergence.pdf` | `driver_field_emulators/code/rebuild/fig_cutoff_paper.py --products ../../products --out ...` | `products/table_{tree,bihalofit}_cut{960,1920,3840,7680,15360}_r4_xi.npz` (ten folds, 2026-08-26) and the O0 sweep | `pyccl` | yes (writes straight into `figures/` by the documented `--out`) | IDENTICAL. |

Cross-cutting findings from the table:

- The FK trap in the brief is confirmed quantitatively for rows 2, 3, 4 and 7; the
  standalone default inputs give an FK channel 4.54x too low at 0.5' (rows 2 to 4) or
  the wrong panel geometry (row 7).
- Three tracked in-tree twins are stale June renders (rows 2, 3, 4 under
  `analysis3/outputs/`), so the repository currently ships two contradictory versions
  of those figures. Phase 2 regenerates them.
- Which conda env was used matters at the 1 pt bounding-box level (rows 1, 2, 12, 13).
  REPRODUCE.md must name the env per figure; the table above records the one that
  reproduces byte-for-byte.
- Rows 8, 9, 10 are not byte-reproducible in either env today (font-metric drift since
  June, same matplotlib 3.10.8 producer string). Decision D2 below.
- Rows 6, 8, 9, 14 need canoes importable; the corr_op callable imports it at module
  level. This is a hard external dependency of the reproduction package.

## 2. Number provenance

All 34 quantitative statements traced (full table with source products:
`reorg_2026-09/phase0/numbers.md`). Every number reproduces from products on disk.
Summary of the ones the brief named, recomputed 2026-09-02:

| statement in the manuscript | recomputed | source products (repo-relative) |
|---|---|---|
| FK at 0.5' is 2.31% of Order-0 at ell_max 15360; ladder 0.48 / 0.97 / 1.53 / 2.02 / 2.31% for 960 to 15360; final doubling +14%; +8 to 9% across 2' to 12'; extrapolation near 2.7% | 2.31%; 0.48/0.97/1.53/2.02/2.31; +14.2%; +7.9/8.4/8.7/9.1% at 2/5/8/12'; 2.71% | `driver_field_emulators/products/table_tree_cut*_r4_xi.npz`, `sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz` |
| BiHalofit ladder does not turn over | 0.78 / 2.56 / 8.01 / 21.5 / 41.1% | `products/table_bihalofit_cut*_r4_xi.npz` |
| FK/Order-0 is 1.3 to 1.7% across 2' to 12' in xi_kappa and xi_+ | xi_kappa 1.74/1.58/1.43/1.31% at 2/5/8/12'; xi_+ 1.79/1.68/1.56/1.47% | `products/table_permclosed_cut15360_permfix_xi.npz` |
| xi_-: about 9% at the 2' edge, near a percent above 3' (commit 44aade3: 8.60/3.91/2.19% at 2.06/2.61/3.31') | 8.60/3.91/2.19% at those grid points; 9.26% interpolated at 2.0' | same |
| xi_kappa_gamma_t suppressed to 1 to 2% | -2.06 to -1.05% (magnitude) over the band | same |
| FK several times FF (kk 6.9 to 3.5, ++ 11.9 to 10.9 across the band); FK above FF in all four | kk 7.0/5.4/4.3/3.3; ++ 12.0/11.8/11.2/10.9; xi_- 51/6.6/6.2/5.2; kg 12.5/7.7/7.3/6.6 (magnitudes) | same plus the FF sweep |
| FF about 0.2% of Order-0 at small gamma; flat plateau; overtakes Order-0 beyond about 200' | 0.24%; plateau spread 0.02% over the last 8 points; crossing at 180' (log interpolation between grid points 145' and 183') | FF and O0 sweeps |
| B-mode falls from 0.95 of the E-mode at ell=60 to 0.42 at ell=1500; about two thirds over 50 to 1500 | 0.951 at ell=60; 0.420 at ell=1500; median 0.629 | FK cut15360 fold through the paper's own Wigner-d transform |
| FF B-mode some thirty times below E | median EE/BB = 31.2 | FF sweep |
| C_EB: exactly zero for FK, FF residual eight orders down | max EB_FK = 0.0; max EB_FF / max EE_FF = 3.9e-9 | same |
| Total Order-2 correction to C_kk over 50 to 1500 "one to two percent, slowly rising" | 1.30% (ell=93) rising to 2.47% (ell=1500); EE+BB 1.49 to 2.45% | same |
| Multi-z: grows with source redshift; 0.9 to 1.1% at z_s=1; percent level at small separations | FK/O0 over 2' to 12': z=1 0.86 to 1.13%, z=1.7 1.07 to 1.41%, z=2.5 1.19 to 1.58%, z=3.2 1.23 to 1.64%, z=4 1.28 to 1.70%; at 0.5': 1.36/1.77/2.04/2.14/2.24% | `analysis3/outputs/multiz_kappa_2pcf_5z.npz`, `products/multiz_ff_real_all5.npz` |
| Monte-Carlo validation: FK to about a percent at 1', FF within 1.3 standard errors | markers/fold 1.011+/-0.005, 1.012+/-0.008, 1.008+/-0.014, 0.999+/-0.033 at 1/2.6/6.7/17.3'; pooled extrapolation 1.0022+/-0.0056 (`headline.py` re-run); FF max pull 1.27 sigma at 1' | `code/mc_fk_complete/_markers_pooled.npz`, `mc_sachs_2pt/outputs/appendix_mc_curve.npz` |
| Kernel quadrature cross-check: about a percent at arcminutes, within 7% out to 17' | ratios 0.992 to 1.015 for 0.5' to 5', 1.069 at 17.3' | `products/fk_kernel_crosscheck.npz` (script needs canoes; npz read, not re-run) |
| zeta slices change sign near 60' in the higher source shells | TTT and Bmod flip at 59' for z >= 4.1; TTP never; Dmod at 59' for every shell | `products/zeta_bands_cut15360_gmax85.npz` |
| n_eff crosses -2 at k = 0.21 h/Mpc (linear) and 8.4 h/Mpc (BiHalofit one-halo) | 0.204 h/Mpc from `PCAMBz0.txt`; the BiHalofit value was not re-derived (source: `driver_field_emulators/note/sec_uv.tex`) | canoes `examples/data/PCAMBz0.txt` |
| Fig. 12 truncation angle 83 deg | 83.33 deg | O0 sweep grid |

**Stop-and-report table.** No numeric contradiction was found. Four statements are
loose at the wording level; they are reported here for the operator, and this
reorganisation proposes no manuscript edit:

| statement | location | recomputed | assessment |
|---|---|---|---|
| FF overtakes Order-0 "beyond about 200'" | `sections/insights.tex`, FF paragraph | crossing at 180' (between grid points 145' and 183') | about 10% high; harmless |
| "a slowly rising one to two percent" total correction in C_kk over 50 to 1500 | Fig. 12 caption in `insights.tex` | 1.30% to 2.47% | upper end rounds to 2.5% |
| amplitudes are "lower bounds, low by of order ten percent" | cutoff subsection | the paper's own extrapolation (2.71%) is 17% above 2.31% | "of order ten percent" is loose, the 2.7% figure itself is exact |
| `t_final = 2313.029` labelled z_s = 5 | every production config | z = 4.70 on the closed-form background used by `D_callable.py` (lambda(5) = 2317.96); a 0.2% offset in lambda to which the fold is sensitive at about 1.1% per Mpc | all three runs and all figures share it, so every ratio is consistent; a reader computing lambda(z=5) themselves would not reproduce the absolute amplitudes to better than about 5%. Decision D3. |

Numbers that trace to one-off scripts rather than to a figure generator (to be
promoted into `reproduce/` in Phase 4): the 2' to 12' band
(`code/mc_fk_complete/numbers_2to12.py`), the Monte-Carlo headline
(`code/mc_fk_complete/headline.py`), the B/E values (originally
`driver_field_emulators/note/sec_eb.tex` machinery; recomputed here through the paper's
EB generator), the kernel cross-check (`code/rebuild/fk_kernel_crosscheck.py`), and
the n_eff crossings (`note/sec_uv.tex`).

## 3. Products: classification

Buckets: **A** paper-dependent (must ship, tracked); **B** intermediate and regenerable
(keep locally, gitignore, PROVENANCE line); **C** exploratory variant (archive under a
dated folder); **D** superseded production product (keep in place with a banner because
external repositories pin it, or archive if not pinned).

### 3.1 `driver_field_emulators/products/` (91 MB, 138 entries)

| item | size | bucket | role and provenance |
|---|---|---|---|
| `table_permclosed_cut15360.npz` | 1.70 MB | A | the production vertex table (perm-closed r4 grid plus extra rows, ell_max 15360; assembled 2026-08-26 09:46 from `pieces/low_permclosed_r4.npz` + `band_tree_permclosed_*` + `band_tree_extra_*`); vendored by SFT-WL-B (sha256 `b400dc5b11e1033b`) |
| `table_permclosed_cut15360_permfix_xi.npz` | 14 KB | A | THE FK 2PCF sweep behind figures 2, 3, 4, 6 and every FK number; `run_fk_variant.py --callable perm_aware --n-jobs 6`, 2026-08-26 09:48 |
| `_variant_table_permclosed_cut15360_permfix/` (config + materialised callable) | 92 KB | A | the run definition of the line above |
| `table_{tree,bihalofit}_cut{960,1920,3840,7680,15360}_r4_xi.npz` (10) | 10 x 14 KB | A | the cutoff ladder folds (figure 17, the 0.48 to 2.31% numbers) |
| `table_{tree,bihalofit}_cut*_r4.npz` (10) | 10 x 1.45 MB | A (regenerable from B `pieces/` in seconds with `assemble.py`) | ladder tables; recommended tracked (14.5 MB) |
| `_variant_table_{tree,bihalofit}_cut*_r4/` (10) | 10 x 72 KB | A | ladder run definitions (configs); caches inside are B |
| `zeta_bands_cut15360.npz`, `zeta_bands_cut15360_gmax85.npz` | 179 + 108 KB | A | figure 7 input (2026-08-27) |
| `multiz_ff_real_all5.npz` (+ per-z `multiz_ff_real_z*.npz`, `multiz_ff_real.npz`, `multiz_ff_real_z5p0ctl.npz`) | 2.7 KB + 6 x 1.4 KB | A | figure 5 FF planes; five single-threaded sweeps of about 3.1 h each plus the z=5 control (bit-identical to the June production FF) |
| `_variant_multiz_conv/` (config) | 88 KB | A | figure 5 FK run definition (perm-aware, cut15360, 5 planes, 1034 s) |
| `order0_multiz_xi.npz`, `source_distances.npz`, `_variant_order0_multiz/` | 77 KB | A | exact Order-0 at 7 source planes (lambda 1330.7 to 2318.0); used by the redshift-trend notes; candidate exact O0 for figure 5 (D4) |
| `multiz_kappa_2pcf_5z_preRevision.npz` | 6 KB | D | June multi-z npz (O0/FF/FK planes; O0 planes reused by the deployed figure) |
| `fk_kernel_crosscheck.npz`, `mean_kappa_plateaus.npz` | 3 KB | A | the appendix cross-check numbers; the FF plateau check in the footnote |
| `logs/` (76 logs) | 1.9 MB | A | build and fold provenance (timings, commands) |
| `pieces/` (LOW branch + octave bands for r4, permclosed, extra) | 50 MB | B | hours of canoes builds; needed to re-assemble any cutoff table; keep locally, gitignore, PROVENANCE |
| `cutoff_vertex.npz`, `cutoff_vertex_tree.npz`, `cutoff_study_partial.npz` | 3.7 + 1.9 + 0.9 MB | B | vertex-level cutoff study behind the audit note; SFT-WL-B vendored `cutoff_vertex.npz` and could not find its producer; the producer is to be identified in Phase 4 (candidates `code/rebuild/analyse_cutoff.py`, `fig_cutoff.py`, `fig_convergence.py`) |
| `redshift/{fold_trends,pershell}.npz`, `cl_ratio_kk*.npz`, `weight_identity.npz`, `solveQ_roundtrip.npz`, `mc_*.npz`, `reduced_shear_zs5.npz`, `grid_comparison_table.txt`, `uv_convergence_bands.txt` | < 100 KB total | B | data behind the audit notes (`notes/*.md`, `note/*.tex`); keep with the audit docs |
| `triples_permclosed_r4.npz`, `triples_permclosed_extra.npz`, `triples_full_r4.npz` | 150 KB | A | grids of the production and ladder tables |
| `triples_full_r{1,2}`, `triples_permclosed_r{1,2}`, `triples_collapsed_only`, `triples_production_only`, `triples_permutation_probe`, `collapsed_triples_dense` | < 200 KB | C | grids of exploratory variants |
| `_variant_collapsed_dense_*` (5), `collapsed_dense_*.npz` (8), `_low_cache_collapsed_dense.npz` | ~1 MB | C | 2026-08-25 dense-collapsed-family variants (cut1000/2000/8192, BiHalofit, low shells) |
| `table_tree_cut1000_r{1,2,4}*.npz` (+ `_banded`, `_extshells`, `_multiz_xi`), `_variant_table_tree_cut1000_*` (5) | ~6 MB | C | resolution-study and extended-shell variants at cut1000 |
| `table_permclosed_cut1000{,_r1,_r2}.npz` and their `_legacy`/`_permfix`/`_permfixmultiz` folds and `_variant_*` folders | ~5 MB | C | the callable-fix study at cut1000 (note: `table_permclosed_cut1000_permfixmultiz_xi.npz` is a cut1000 product, factor about 5 low; never plotted in the paper) |
| `table_reproduce_deployed*.npz`, `_variant_table_reproduce_deployed/` | 1.4 MB | C | the rebuild-reproduces-deployed check (median relative difference 3e-16) |
| `_ff_caches/`, `_ff_dt1/`, `xi_FF_dt1.npz`, `_l2_lambda_grid_*.npz`, `xi_kappa_channels_preRevision.pdf` | < 200 KB | B or C | caches and one archived figure |

### 3.2 `SFT-lensing-paper-analyses/sachs_sft/sftwick_outputs/2PCF/`

| run folder | bucket | note |
|---|---|---|
| `C_corr_op_O0/` (xi 2026-06-01, config) | A, live | Order-0 for every figure and number |
| `C_corr_op_K_limber_FF/` (xi 2026-06-01, config) | A, live | FF for every figure and number; the z=5 control re-run of 2026-08-27 reproduced it bit for bit |
| `C_corr_op_K_limber_FK/` (xi 2026-06-10, config, `_PROD_REF_2026-06-10.npz` byte-identical twin) | D | cut1000, old callable, corner-frozen vertex; vendored by SFT-WL-B (sha256 `c3eae8cb8515668`); the config is the template every variant run copies. Keep in place with a SUPERSEDED banner; archive the `_PROD_REF` twin. |
| `.cache_expand_*`, `.cache_propagators_*` | B | joblib caches, gitignored already |

### 3.3 Other data

| item | bucket | note |
|---|---|---|
| `callables/C_propagator/corr_op/corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz` (5.5 MB, tracked) | A | the C propagator table |
| `callables/kappa3_vertex/equal_time_limber/equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz` (1.3 MB, tracked, vendored) | D | cut1000 table, still the default of the production callable and of `build_equal_time_limber_table.py`; keep, banner |
| `callables/kappa3_vertex/equal_time_limber/ell_band_decomp_results.npz` (117 KB, tracked) | D | June cut1000 input of figure 7 |
| `callables/kappa3_vertex/equal_time_limber/_smoke_equal_time_limber.npz` (120 KB, untracked) | C | 2026-08-25 smoke build; archive |
| `analysis3/outputs/multiz_kappa_2pcf_5z.npz` (tracked, modified) | A | figure 5 intermediate (see row 5) |
| `mc_sachs_2pt/outputs/appendix_mc_curve.npz` (tracked, modified) | A | figure 6 cache (FF Monte-Carlo, corrected FK line, FK markers) |
| `code/mc_fk_complete/_markers_pooled.npz`, `_gate2_highstat*.npz`, `_papergrid.npz`, `_blk*.npz`, `_reseed*.npz`, `_gate*.npz`, `_exact_*.npz` | A (markers) / B (blocks) | the FK Monte-Carlo markers and the seed blocks they pool |
| `code/mc_fk_complete/psd_fix/fold_bias/tables/corr_op_gl24_plus6.npz` (11 MB) and the rest of `psd_fix/` (14 MB) | B/C | review side-study; gitignore the npz, keep the docs |
| `code/normalisation/data/` (68 KB) | B | inputs of the normalisation chain check |
| `talk/assets/figures/_data/multiz_components.npz` (21 KB, tracked in talk) | A (provenance) | 20-plane talk cache; the deployed figure 5 Order-0 planes are PCHIP-in-z of it; vendor a copy with PROVENANCE |
| `talk/assets/figures/_data/{fk_multiz,mc_full_xi,fk_mc}*.npz` | C (for the paper) | talk-only |
| root `outputs/review_3pt_vs_pyccl_FINAL.npz` (2026-06-09) | C | archive |
| `figures/_archived/` (13 items, May 2026, untracked) | D | archive into the analyses repo |
| `figures/archived_2026-08-26/`, `figures/archived_2026-08-28/` (tracked, on Overleaf) | D | keep (decision D9) |
| `figures/appendix_3cumulant_fastnc.pdf`, `figures/cumulant_hierarchy.pdf` (tracked orphans) | D | decision D14 |
| `driver_field_emulators/papers/` (9 arXiv PDFs, 21 MB) | reference | keep locally, gitignore |

## 4. Documents: classification and contradictions

Status: **current** (true today, keep and maintain), **historical** (true when written,
keep for provenance, mark with a dated header line), **obsolete** (wrong or superseded,
archive).

| document | status | note |
|---|---|---|
| `CLAUDE.md`, `AGENTS.md` (root; identical except the memory path line) | current, needs update | "5 of the 17 figures" is wrong (7); Repository layout section to be rewritten in Phase 4 |
| `SFT-lensing-paper-analyses/figure_manifest.md` | current, needs update | row 5 misstates the talk dependency (code import, not npz); rows 15, 16 now have a generator; gitignored, must be tracked |
| `SFT-lensing-paper-analyses/README.md` | current, needs update | says "16 paper figures", "figure_manifest git-ignored" |
| `SFT-lensing-paper-analyses/CLEANUP_REPORT.md`, `CLEANUP_HANDOFF.md` (2026-06-17) | historical | the precedent format for the changelog; "16 figures" was correct then |
| `SFT-lensing-paper-analyses/docs/{onepz4_localization,pz4_fix_applied}.md` (2026-06-09) | historical | design records of the (1+z)^4 fix; gitignored today (`/docs/`), should be tracked |
| `mathematica/docs/*.md` | current | proof records |
| `sachs_sft/README.md`, `callables/README.md`, `analyses/README.md`, `scripts/README.md` (2026-05-30) | historical | describe the four-callable skeleton of May; two of the four are archived |
| `sftwick_outputs/README.md` and the three run READMEs, `equal_time_limber/README.md`, `corr_op/README.md`, `slices/README.md`, `scripts/gen_sftwick_configs.py` docstring | current but carry the stale `+3.086e-5` baseline stamp (contradiction 1) | fix the stamp text in Phase 3/4 |
| `equal_time_limber/{ELL_BAND_DECOMPOSITION,MATH_FLOW}.md` | historical | numbers at cut1000 |
| `sftwick_outputs/PROMPT_run_2PCF_three_runs.md`, `shear3pcf_fastnc/PROMPT*.md`, `driver_field_emulators/PROMPT_*.md` | historical (agent briefs) | keep out of the shipped package (`_archive/` or `docs/agent_briefs/`) |
| `analysis3/FK_GAMMA_TRENDS_ACROSS_PANELS.md`, `analysis3/WHY_FK_BREAKS_KAPPA_GAMMA_EQUIVALENCE.md` (2026-06-01) | obsolete | both explain features that were the cosine-grid artifact (contradiction 6) |
| `mc_sachs_2pt/{DESIGN,FF_NOTES,IMPLEMENTATION}.md` | current | the FF estimator is unchanged |
| `mc_sachs_2pt/FK_NOTES.md` | obsolete | describes the `simulate_fk_vr` calibration that was withdrawn (placement share 0.21) |
| `shear3pcf_fastnc/*.md` | historical (backup analysis, not in paper) | already labelled `BACKUP_NOT_IN_PAPER.md` |
| `driver_field_emulators/{README,SUMMARY,HANDOVER,SESSION_STATE_2026-08-26,CHANGES_OUTSIDE_THIS_FOLDER}.md` | historical with high provenance value | `CHANGES_OUTSIDE_THIS_FOLDER.md` is the only audit trail of what was deployed from where; SESSION_STATE carries two now-reversed decisions (FK markers, B/E max ratio) |
| `driver_field_emulators/notes/*.md` (24) | historical, several with stale absolute FK numbers | contradiction 1; the stamp notes named in the brief are `assumptions_and_limitations.md:55` and `input_ready_b_to_zeta_spec.md:170-173`; `finding_ell_max_uv_sensitivity.md` is self-labelled SUPERSEDED |
| `driver_field_emulators/note/` (LaTeX audit note, 22 pp, `main.pdf` 2026-08-27) | historical, the referee-grade record | keep tex + pdf; drop aux/log/fdb/fls/out/toc (regenerable) |
| `driver_field_emulators/revision_plan/` (letter, email, plan) | historical, **private** (editor correspondence) | must not be pushed to the public GitHub package (decision D15) |
| `driver_field_emulators/code/{README,rebuild/README,callable_fixed/README}.md`, `code/mc_fk_complete/{NOTES,SPEC,PAPER_PROPOSAL,...}.md`, `review/`, `psd_fix/*.md` | current for their code | move with the code |
| root `analysis_setup.md` (2026-05-19) | obsolete | pre-dates both reorganisations; references `sections/discussion.tex`, a PRL companion, the dead `/Users/zzhang/projects/canoes`, a 48-page draft |
| root `findings.md`, `progress.md`, `task_plan.md` (2026-06-16) | obsolete | an agent's Figure-16 verification scaffold, superseded by the 2026-08-26 finding that the FK estimator it verified measures a placement share |
| `figures/_archived/driving_field_schematics.README.md` (2026-05-17) | historical | references `scripts/canoes_pipeline/scripts/plot_driving_field_*.py` (now in `_archive/`) |
| `talk/README.md` | current for the talk | no paper dependency stated |

**Contradictions between documents (numbers that fight each other):**

1. **The FK kappa-kappa(0.5') baseline appears as five different numbers.**
   `+3.086e-5` (all `sftwick_outputs` READMEs, `equal_time_limber/README.md`, the
   production callable docstring, `gen_sftwick_configs.py`, `corr_op/README.md`) is the
   pre-h^4-fix stamp; `+1.219e-4` (`build_equal_time_limber_table.py` docstring) is the
   pre-Jacobian stamp; `+4.29e-6` (`notes/assumptions_and_limitations.md:55`,
   `notes/input_ready_b_to_zeta_spec.md:170-173`, `code/README.md`, `b_model/tree.py`,
   `notes/implementation_roadmap.md`) is the deployed cut1000 value with the old
   callable and linear-in-lambda interpolation, presented as "the true tree baseline"
   without its ell_max; `+3.84e-6` (`HANDOVER.md`, `SUMMARY.md` section 4d,
   `theory_redshift_dependence.md`, `finding_grid_artifact_measured.md`) is cut1000 with
   log-PCHIP lambda interpolation; the manuscript's value is `+1.948e-5` (2.31% of
   Order-0) at ell_max 15360 with the permutation-aware callable. Recomputed today:
   deployed sweep `+4.2877e-6` (0.508%), converged fold `+1.9480e-5` (2.31%).
   Fix (Phase 3): a dated correction block in the two stamp notes stating the
   ell_max, callable and lambda-interpolation of each number, pointing at the current
   run folder; the `+3.086e-5` / `+1.219e-4` stamps in the run READMEs and docstrings
   get the same treatment (those files are not vendored except the production
   callable, whose docstring stays untouched; its folder gets a DEPRECATED note).
2. **Figure count 16 vs 17**: `CLEANUP_REPORT.md`, `HANDOVER.md` ("six of the sixteen"),
   `SFT-lensing-paper-analyses/README.md`, `_archive/cleanup_2026-06-17/README.md` say 16
   (correct in June); `figure_manifest.md` and the build say 17.
3. **FK Monte-Carlo markers**: `SESSION_STATE_2026-08-26.md` and
   `CHANGES_OUTSIDE_THIS_FOLDER.md` item 7 record their removal; item 9 (2026-08-28)
   records their restoration from `code/mc_fk_complete/`; the deployed figure has them.
   `mc_sachs_2pt/FK_NOTES.md` still documents the withdrawn estimator as if valid.
4. **B/E**: the memory note says "B/E ~ 0.36"; `SESSION_STATE` says "B/E = 0.47" (ratio
   of maxima); `note/sec_eb.tex` says the callable fix moved it "from 0.36 to 0.42";
   the manuscript says 0.95 at ell=60 falling to 0.42 at ell=1500, about two thirds on
   average. These are different statistics of the same curve; only the manuscript's
   per-ell statement is recomputed here (0.951, 0.420, median 0.629).
5. **Where |E-B|/|E+B| reaches unity**: 22' (memory note) vs 31' (`SESSION_STATE`,
   `note/sec_eb.tex`). Not quoted in the manuscript.
6. **June analysis3 notes vs the August finding**: "FK suppressed by 5 to 6 orders in
   xi_+/- and kappa-gamma" and "FK flat until very large angles" were the corner-frozen
   vertex; the manuscript now quotes FK in all four 2PCFs.
7. **`t_final = 2313.029` labelled z_s = 5** (all configs, all figures) vs z = 4.70 on
   the closed-form background (`HANDOVER.md`, memory); see the stop-and-report table.
8. **`figure_manifest.md` row 5** says the compute step reads the talk npz; it imports
   the talk script's code (section 0, item 3).
9. **`CLAUDE.md`** says `driver_field_emulators` holds "the generators for 5 of the 17
   figures"; it is 7.

## 5. Target layout

Principle: the nested repository `SFT-lensing-paper-analyses/` (GitHub
`StatFieldTheory/SFT-lensing-paper-analyses`) becomes the single reproduction package;
the paper repository keeps only the manuscript; `driver_field_emulators/` is absorbed
into the package with its history preserved from a verbatim import commit; `talk/`
becomes a consumer of the package instead of a supplier.

```
STF_lensing/                                  paper repo (Overleaf master): main.tex, sections/, figures/, biblio.bib only
  REORG_PLAN.md                               (this file; moves into the package in Phase 4)
  driver_field_emulators/                     after Phase 2: a one-file MOVED.md pointer for external sessions (folder itself gitignored)
  SFT-lensing-paper-analyses/                 the reproduction package (nested repo)
    REPRODUCE.md  figure_manifest.md  MOVED_PATHS.md  CHANGELOG_REORG_2026-09.md   (Phase 4)
    reproduce/                                  Phase 2: regen_figure.py (per-figure driver, env-aware), check_figures.py (md5 + raster),
                                                         numbers.py (the quoted-number tracer), deploy.py (copy into ../figures/ under paper names)
    reorg_2026-09/phase0/                       this session's evidence (already present)
    figures/                                    schematic generators (rows 8 to 13), unchanged; outputs redirected to figures/outputs/ + deploy step
    sachs_sft/
      callables/
        C_propagator/corr_op/                   unchanged
        kappa3_vertex/
          equal_time_limber/                    build script, ell_band_decomp.py, plot_zeta_figures.py, the cut1000 table (vendored, unchanged),
                                                the production callable (vendored, byte-identical) + DEPRECATED.md sibling
          equal_time_limber_cut15360_permaware/ <- driver_field_emulators/code/callable_fixed/ (perm_aware_kappa3_callable.py byte-identical,
                                                tests/, README.md, compare_callables.py) + table_permclosed.npz (= products/table_permclosed_cut15360.npz,
                                                the name the callable resolves) + PROVENANCE.md
          rebuild/                              <- driver_field_emulators/code/rebuild/ (build_low, build_band, assemble, run_queue, jobs_*.json,
                                                sweep_cutoff.sh, fig_cutoff_paper.py, probes) + make_dense_triples.py, make_collapsed_triples.py,
                                                run_fk_variant.py + products/{triples_*, table_*_cut*_r4.npz, logs/} ; pieces/ gitignored (PROVENANCE)
          b_model/                              <- driver_field_emulators/code/b_model/
      sftwick_outputs/2PCF/
        C_corr_op_O0/  C_corr_op_K_limber_FF/   unchanged (live); README banners updated
        C_corr_op_K_limber_FK/                  SUPERSEDED banner; xi npz unchanged (vendored); _PROD_REF twin archived
        C_corr_op_K_limber_FK_cut15360_permfix/ NEW live FK run: config_L2.yaml (from _variant_table_permclosed_cut15360_permfix, paths rewritten),
                                                xi_C_corr_op_K_limber_FK_cut15360_permfix.npz (= products/table_permclosed_cut15360_permfix_xi.npz),
                                                README.md with provenance and the 2026-08-26 fold log
        ladder_cut{960..15360}_{tree,bihalofit}/ the ten cutoff-ladder folds with their configs (figure 17)
        multiz/                                 order0 (7 planes), fk_cut15360_permfix (5 planes), ff_real (5 planes + z5 control) with configs and logs
      analyses/
        analysis1/                              unchanged
        analysis3/                              FK_NPZ default -> the new run folder; compute_multiz_kappa_2pcf.py imports multiz_sweep.py (moved copy of
                                                talk/scripts/run_multiz_components.py); inputs/multiz_components_talk_2026-06-10.npz vendored with PROVENANCE
        mc_sachs_2pt/                           FK default -> new run folder; make_val_figure logic folded into fig_xi_channels.py (--fk-markers)
        mc_fk_complete/                         <- driver_field_emulators/code/mc_fk_complete/ (markers producer); psd_fix large npz gitignored
        revision_2026-08/                       <- driver_field_emulators/code/figures_corrected/ (regenerate*.py, make_val_figure.py, run_ff_one_lambda.py,
                                                run_ff_multiz_serial.py, extend_ff_markers*.py); path arithmetic fixed
        fk_audit_2026-08/                       <- driver_field_emulators/code/{probe_*, redshift_*, sweep_ell_high_max, plot_fk_variants*, compare_fk_variants,
                                                mc_crosscheck_dense, verify_*.wl, normalisation/, reduced_shear/} + products/{redshift/, cl_ratio*, cutoff_vertex*,
                                                weight_identity, solveQ_roundtrip, mc_*, reduced_shear_zs5, *.txt}
        shear3pcf_fastnc/                       unchanged (backup)
      scripts/                                  unchanged; run_FF_single.py gains a refuse-to-unlink guard (Phase 2)
    docs/
      fk_audit_2026-08/                         <- driver_field_emulators/{notes/, note/ (tex+pdf), README, SUMMARY, HANDOVER, SESSION_STATE, CHANGES_OUTSIDE..., task_a...}
      agent_briefs/                             <- every PROMPT_*.md (historical)
      papers/                                   <- driver_field_emulators/papers/ (gitignored)
      private/                                  <- driver_field_emulators/revision_plan/ (gitignored: editor correspondence)
    mathematica/                                unchanged (+ beam_width_regulator.wl committed)
    _archive/2026-09_reorg/                     README.md catalogue + dfe_products_exploratory/ (bucket C), root_docs_2026-05/, figures_may2026/,
                                                latex_residue/, smoke_builds/, _PROD_REF twin
  talk/                                         own repo; scripts/run_multiz_components.py becomes a thin wrapper importing the package's multiz_sweep.py;
                                                tracked __pycache__ removed, .gitignore added (decision D12)
```

### 5.1 Migration map (draft of `MOVED_PATHS.md`)

Every path that SFT-WL-B's `vendor/*/PROVENANCE.md` or `docs/decision_log_bmode.md` `[P]`
entries reference is listed; vendored files stay byte-identical (only their path
changes). rezeta references the tree only by its root and by `sections/appendix.tex` and
`driver_field_emulators/notes/input_ready_b_to_zeta_spec.md`.

| old path (under `STF_lensing/`) | new path | content |
|---|---|---|
| `driver_field_emulators/code/callable_fixed/perm_aware_kappa3_callable.py` | `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable.py` | byte-identical (vendored, sha `4d7fefaa129e1f1e`) |
| `driver_field_emulators/code/callable_fixed/{compare_callables.py,README.md,tests/test_perm_aware.py}` | same folder | byte-identical (vendored) |
| `driver_field_emulators/products/table_permclosed_cut15360.npz` | `.../equal_time_limber_cut15360_permaware/table_permclosed.npz` | byte-identical (vendored, sha `b400dc5b11e1033b`); renamed to the name the callable resolves, old name recorded in PROVENANCE.md |
| `driver_field_emulators/products/table_permclosed_cut15360_permfix_xi.npz` | `SFT-lensing-paper-analyses/sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz` | byte-identical |
| `driver_field_emulators/products/_variant_table_permclosed_cut15360_permfix/config__variant_table_permclosed_cut15360_permfix.yaml` | `.../C_corr_op_K_limber_FK_cut15360_permfix/config_L2.yaml` | paths rewritten (callable folder, output, caches); original kept in `_archive/2026-09_reorg/` |
| `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber/equal_time_limber_kappa3_callable.py` | unchanged | byte-identical (vendored, sha `4b9c47cc8771ab31`); `DEPRECATED.md` added beside it |
| `.../equal_time_limber/equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz` | unchanged | byte-identical (vendored) |
| `SFT-lensing-paper-analyses/sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK.npz` | unchanged | byte-identical (vendored, sha `c3eae8cb8515668`); folder README gets a SUPERSEDED banner |
| `.../C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK_PROD_REF_2026-06-10.npz` | `SFT-lensing-paper-analyses/_archive/2026-09_reorg/sftwick_prod_ref/` | byte-identical twin of the line above |
| `SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3/{plot_analysis3_cl_decomposition.py,plot_cl_EB_polarization.py,_plot_style.py}` | unchanged paths | edited in Phase 2 (default `FK_NPZ`); SFT-WL-B vendored the 2026-08-27/29 versions, which the Phase-1 baseline commit preserves in history |
| `SFT-lensing-paper-analyses/sachs_sft/sftwick_outputs/2PCF/{C_corr_op_O0,C_corr_op_K_limber_FF}/xi_*.npz` | unchanged | byte-identical (vendored) |
| `driver_field_emulators/products/cutoff_vertex.npz` | `SFT-lensing-paper-analyses/sachs_sft/analyses/fk_audit_2026-08/products/cutoff_vertex.npz` | byte-identical (vendored by SFT-WL-B, sha `c429cb9b...`) |
| `driver_field_emulators/notes/finding_grid_artifact_measured.md`, `notes/emulator_candidates.md`, `notes/input_ready_b_to_zeta_spec.md`, other `notes/*.md` | `SFT-lensing-paper-analyses/docs/fk_audit_2026-08/notes/<same name>` | unchanged except dated correction blocks in the two stamp notes (Phase 3) |
| `driver_field_emulators/note/` | `SFT-lensing-paper-analyses/docs/fk_audit_2026-08/note/` | tex + pdf unchanged; build residue dropped |
| `driver_field_emulators/code/rebuild/*` (incl. `assemble.py`) | `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/rebuild/*` | `_REPO` path arithmetic fixed where present |
| `driver_field_emulators/code/figures_corrected/regenerate.py` and siblings | `SFT-lensing-paper-analyses/sachs_sft/analyses/revision_2026-08/` | path arithmetic fixed |
| `driver_field_emulators/code/mc_fk_complete/` | `SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/` | unchanged |
| `driver_field_emulators/code/b_model/` | `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/b_model/` | unchanged |
| `driver_field_emulators/products/table_{tree,bihalofit}_cut*_r4{,_xi}.npz` and `_variant_*` configs | `SFT-lensing-paper-analyses/sachs_sft/sftwick_outputs/2PCF/ladder_cut<N>_<model>/` | byte-identical |
| `driver_field_emulators/products/pieces/` | `.../kappa3_vertex/rebuild/pieces/` (gitignored, PROVENANCE) | byte-identical, local only |
| `driver_field_emulators/products/{zeta_bands_cut15360.npz,zeta_bands_cut15360_gmax85.npz}` | `.../kappa3_vertex/equal_time_limber/outputs/` | byte-identical |
| `driver_field_emulators/products/multiz_ff_real*.npz`, `order0_multiz_xi.npz`, `source_distances.npz`, `_variant_multiz_conv`, `_variant_order0_multiz`, `_ff_caches` | `.../sftwick_outputs/2PCF/multiz/` | byte-identical |
| `driver_field_emulators/products/<bucket C items>` | `SFT-lensing-paper-analyses/_archive/2026-09_reorg/dfe_products_exploratory/` | byte-identical |
| `driver_field_emulators/{SUMMARY,HANDOVER,README,SESSION_STATE_2026-08-26,CHANGES_OUTSIDE_THIS_FOLDER}.md` | `SFT-lensing-paper-analyses/docs/fk_audit_2026-08/` | unchanged |
| `driver_field_emulators/revision_plan/` | `SFT-lensing-paper-analyses/docs/private/revision_plan_2026-08/` (gitignored) | unchanged |
| `driver_field_emulators/papers/` | `SFT-lensing-paper-analyses/docs/papers/` (gitignored) | unchanged |
| `talk/assets/figures/_data/multiz_components.npz` | unchanged in talk; copy at `SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3/inputs/multiz_components_talk_2026-06-10.npz` | byte-identical copy |
| `talk/scripts/run_multiz_components.py` | copy at `SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3/multiz_sweep.py` (absolute paths made `__file__`-relative); talk keeps a wrapper | edited copy |
| `sections/*.tex`, `figures/*.pdf`, `main.tex`, `biblio.bib` | unchanged | read-only |
| root `analysis_setup.md`, `findings.md`, `progress.md`, `task_plan.md`, `outputs/` | `SFT-lensing-paper-analyses/_archive/2026-09_reorg/root_docs_2026-05/` | unchanged |
| `figures/_archived/*` | `SFT-lensing-paper-analyses/_archive/2026-09_reorg/figures_may2026/` | unchanged |
| `figures/archived_2026-08-26/`, `figures/archived_2026-08-28/` | unchanged (tracked, Overleaf) | decision D9 |
| `SFT-lensing-paper-analyses/figure_manifest.md` | unchanged path, tracked from Phase 1 | rewritten in Phase 4 |
| `SFT-lensing-paper-analyses/REPRODUCE.md` (referenced prospectively by SFT-WL-B) | new in Phase 4 | |

### 5.2 Default-path corrections (Phase 2)

1. `analysis3/plot_analysis3_nlo_decomposition.py` and `plot_analysis3_cl_decomposition.py`:
   `FK_NPZ` default to the new run folder (`plot_cl_EB_polarization.py` and
   `plot_multiz_kappa_2x5.py` inherit through `C`). The standalone-run trap disappears.
2. `mc_sachs_2pt/fig_xi_channels.py`: `_load_kk("C_corr_op_K_limber_FK")` reads the new
   run; `--fk-markers <npz>` option absorbs `make_val_figure.py`; `driver_stats.py`
   imports the permutation-aware callable (affects only a Monte-Carlo recompute, not the
   cached figure; documented).
3. `equal_time_limber/plot_zeta_figures.py`: `_NPZ` default to the cut15360 gmax85 table
   (moved next to it); `_FIGDIR` to its own `outputs/`; deployment through
   `reproduce/deploy.py`.
4. `analysis3/compute_multiz_kappa_2pcf.py`: import `multiz_sweep` from its own folder
   instead of `talk/scripts` by absolute path; `plot_multiz_kappa_2x5.py`: `_FF_REAL`
   path to the new location.
5. `rebuild/fig_cutoff_paper.py`: `--products` default to the ladder run folders.
6. The three live `config_L2.yaml` (and the promoted ones) carry absolute
   `/Users/zzhang/...` paths; callable paths become YAML-relative (sft-wick resolves
   them against the YAML directory), output and cache paths stay managed by the runner
   (which chdirs to the YAML directory). Without this a reader cannot run a sweep at all.
7. `scripts/run_FF_single.py`: refuse to unlink an output inside a run folder that
   contains a `PRODUCTION` marker unless `--force` is given (the 2026-08-25 data-loss
   trap).
8. Schematic generators under `SFT-lensing-paper-analyses/figures/` and
   `make_response_corr_operator_figures.py`: default output to `figures/outputs/`,
   deployment only through `reproduce/deploy.py` (one path into the paper's `figures/`,
   with an md5 check against `REPRODUCE.md`).
9. The `_archive/canoes_pipeline/scripts/plot_coord_maps.py` generator: promote a copy to
   `SFT-lensing-paper-analyses/figures/make_cosmology_coord_maps.py` with the two canoes
   constants and the archived `cosmo.py` functions vendored (they are 40 lines), so the
   figure no longer needs canoes or the old package name. The archive copy stays.
10. The two Feynman diagrams: export cells 29 and 38 of the sft-wick notebook to
    `SFT-lensing-paper-analyses/figures/make_feyndiag.py` (pinning sft-wick commit
    `dd03748` and the `Action` definition), and record the sft-wick commit in
    REPRODUCE.md. The notebook itself stays in sft-wick.

## 6. Phases, risks, rollback

### Phase 1: version-control baseline (no moves)

1. Paper repo: `git fetch` and fast-forward `master` to `1f0ec58` (Overleaf's
   `\revhlfalse`). Then tag `pre-reorg-2026-09` in all three repos.
2. Analyses repo: commit the 18 items in four groups (Appendix C): canoes-path
   resolver and build flags; revision-era generator changes; regenerated outputs and
   caches; the new `.wl`. `_smoke_equal_time_limber.npz` is not committed (archived in
   Phase 2). `git push` is the operator's call (public GitHub).
3. Import `driver_field_emulators/` **by copy** into
   `SFT-lensing-paper-analyses/driver_field_emulators/` (137 MB on disk, about 45 MB in
   git after the `.gitignore` for `pieces/`, `papers/`, `revision_plan/`, `psd_fix/**/*.npz`,
   caches, `.benchmarks`, `note/*.aux|log|fdb_latexmk|fls|out|toc`), verify with an md5
   manifest of every file, commit "import driver_field_emulators verbatim (2026-09-02
   state)", tag `dfe-imported-2026-09`. The original folder stays in place, untouched,
   until Phase 4 delivers `MOVED_PATHS.md`; then it is reduced to a pointer file.
   Alternative (D1): `git init` inside the original folder as its own repository.
4. Root `.gitignore`: stop ignoring `*.md` and `*.markdown`; add `main_clean.pdf`,
   `main *.*`, `meascol.*`; keep ignoring `driver_field_emulators/` until it becomes a
   pointer. Analyses `.gitignore`: stop ignoring `/figure_manifest.md`,
   `/CLEANUP_REPORT.md`, `/CLEANUP_HANDOFF.md`, `/docs/`; add the large-product and
   private-document rules above. Commit.
5. talk repo (if D12 approved): `git rm --cached` the 7 tracked `__pycache__/*.pyc`,
   add a `.gitignore`.

Rollback: `git reset --hard pre-reorg-2026-09` in any repo; the import is an added
directory (`git rm -r`); nothing outside git changes.

### Phase 2: moves and default paths

Per group, in this order, one commit each: (a) callables (new permaware folder, rebuild,
b_model); (b) run folders (new FK run, ladder, multiz, banners); (c) analyses
(mc_fk_complete, revision_2026-08, fk_audit_2026-08, multiz_sweep, vendored talk cache);
(d) docs (fk_audit notes and note, agent briefs, private, papers); (e) archive bucket C
and the root leftovers with a catalogue README; (f) default-path edits of section 5.2;
(g) regenerate all 17 figures with `reproduce/regen_figure.py` into a scratch dir and
run `check_figures.py`: rows 1 to 7, 11 to 14, 17 must be IDENTICAL, rows 8 to 10
EQUIVALENT, rows 15 to 16 md5-matched; (h) if D2 approved, deploy fresh renders of rows
8, 9, 10 (old versions to `figures/archived_2026-09-<dd>/`) and regenerate the three
stale `analysis3/outputs/*.pdf` twins; (i) update `CLAUDE.md` layout section.

Risks and rollback: every move is `git mv` (revertible with `git revert`); archive moves
are `mv` with a catalogue; deployed PDFs replaced only with the old copy archived and
tracked; `driver_field_emulators/` original untouched until Phase 4.

### Phase 3: stale products and stamps

1. Pin the external runtime state first: record sft-wick and canoes commit hashes and
   save `git diff` of both working trees as patch files under
   `SFT-lensing-paper-analyses/reorg_2026-09/env/`, plus `conda list --explicit` of both
   envs.
2. Re-run the FK 2PCF sweep once through `scripts/run_FF_single.py` with a materialised
   config (never the production config directly), `n_jobs 6`, one fold at a time, into a
   scratch output; compare with the 2026-08-26 product (expected bit-identical, as the
   FF control was on 2026-08-27). About 2.5 min. Promote the verified file.
3. Optional (D4): exact Order-0 at the five multi-z planes at the figure's own lambda
   values (about 2 min) against the PCHIP planes.
4. Stamp notes: dated correction blocks with the ell_max qualifier and a pointer to the
   new run folder; the same for the `+3.086e-5` and `+1.219e-4` stamps in the run READMEs
   and docstrings (not in the vendored callable).
5. Re-run `reproduce/numbers.py`; put new and old values side by side; anything that
   changes against the manuscript goes to a stop-and-report table for the operator.
   Nothing longer than an hour is run; the FF multi-z sweeps (3 h each) and the table
   builds are documented, not repeated.

Risk: the `run_FF_single.py` unlink; mitigated by materialised configs and the guard.
Risk: memory; one fold at a time, `n_jobs 6`. Rollback: products are new files; the
2026-08-26 product stays.

### Phase 4: documents

`REPRODUCE.md` (envs with absolute interpreters and the canoes/sft-wick pins, one
command per figure with the env that reproduces byte-for-byte, expected runtime, expected
md5), `figure_manifest.md` rewritten and tracked, `CHANGELOG_REORG_2026-09.md` in the
`CLEANUP_REPORT.md` format, `MOVED_PATHS.md` (section 5.1 finalised, including every
`[P]` reference in SFT-WL-B's decision log: `sections/conclusion.tex:49-52,63-71`,
`sections/insights.tex:225-240`, the production callable `:155-157`,
`notes/finding_grid_artifact_measured.md`, the superseded FK sweep), `CLAUDE.md` and
`AGENTS.md` layout sections, memory entries. `driver_field_emulators/` reduced to
`MOVED.md`.

### Phase 5: acceptance

`latexmk -g -pdf -halt-on-error main.tex` exit 0, 32 pages, 17 figures (baseline
verified today from a scratch copy of the sources: exit 0, 32 pages, 17
`\includegraphics`, 17 `fig` labels, 0 undefined references or citations); every
figure reproduced through `REPRODUCE.md`; clean `git status` in all three repos;
`MOVED_PATHS.md` delivered.

### Risk register

| id | risk | mitigation | rollback |
|---|---|---|---|
| R1 | Overleaf sync: every push of `master` changes the Overleaf project; the paper repo is behind by one commit | fast-forward first; touch only `.gitignore`, new docs and re-deployed figures on `master` | `git revert` |
| R2 | `run_FF_single.py` unlinks the output path before running; a config pointing into production destroys the deployed sweep (happened 2026-08-25) | materialised configs only; guard in Phase 2 | `_PROD_REF` twin and git history |
| R3 | SFT-WL-B and rezeta pin SHA-256 and paths of 15 files | byte-identical moves; `MOVED_PATHS.md` | n/a |
| R4 | canoes is a runtime import for figures 6, 8, 9, 14 and all sweeps; it lives on a feature branch with an uncommitted patch and moved twice already | pin commit + patch file; document `CANOES_ROOT`; vendor the two constants for figure 14 | n/a |
| R5 | sft-wick working tree is dirty; the deployed sweeps were produced by some state of it | pin commit + patch file before any re-run; treat a non-bit-identical re-run as a signal, not noise | keep the 2026-08-26 product |
| R6 | memory: canoes HIGH builds fan-out caused two OOM reboots | no builds in this reorganisation; folds one at a time | n/a |
| R7 | environment drift produces sub-point bounding-box differences; byte-identical md5s in REPRODUCE.md may not survive a conda update | record env per figure and `conda list --explicit`; the checker reports EQUIVALENT separately from DIFFERS | n/a |
| R8 | talk repo edits break the talk | wrapper keeps the old entry point; only the import line changes | `git revert` in talk |
| R9 | figure archive policy touches Overleaf-visible files | D9 default keeps them | n/a |
| R10 | two diverged `_plot_style.py` copies | leave as-is in Phase 2 (both reproduce their figures); unify later only with the checker green | n/a |
| R11 | absolute `/Users/zzhang/...` paths in configs and four scripts make the package unrunnable elsewhere | section 5.2 items 4 to 6 | `git revert` |
| R12 | the public GitHub package would receive editor correspondence and the 21 MB of arXiv PDFs | `docs/private/` and `docs/papers/` gitignored | n/a |
| R13 | importing dfe by copy doubles 137 MB on disk for the duration of Phase 2 to 4 | acceptable; original removed only after `MOVED_PATHS.md` is delivered | n/a |
| R14 | the local `block-no-verify` hook blocks shell commands whose text contains "commit" and "--" (it blocked writing this file through a heredoc) | single-line commit messages; documents written with the file tool | n/a |

## 7. Decisions requested from the operator

| id | decision | recommendation |
|---|---|---|
| D1 | `driver_field_emulators` into the analyses repo (by verbatim import commit) or its own repository | into the analyses repo |
| D2 | Re-deploy figures 8, 9, 10 from a pinned env so their md5s are exact (visually identical; old copies archived and tracked), and regenerate the three stale `analysis3/outputs/*.pdf` twins | yes to both |
| D3 | The `t_final = 2313.029` (z = 4.70) convention: document in REPRODUCE.md only, or quantify the effect with a cheap Order-0 re-run at 2317.96 in Phase 3 (no manuscript change either way) | document, and quantify if cheap |
| D4 | Figure 5 Order-0 planes: keep the PCHIP-of-talk-cache provenance (vendored copy) or compute exact Order-0 at the five planes (about 2 min) and, if invisible, re-deploy | compute and compare; re-deploy only if the operator agrees after seeing the difference |
| D5 | Re-run the FK 2PCF sweep once for verification (2.5 min) before promoting the 2026-08-26 product | yes |
| D6 | Superseded cut1000 FK run folder: keep in place with a SUPERSEDED banner (vendored) vs archive | keep with banner |
| D7 | `_variant_*` folders (25): promote 3 (cut15360 permfix, multiz_conv, order0_multiz) plus the 10 ladder configs; archive the other 12 | as stated |
| D8 | `pieces/` (50 MB) gitignored with PROVENANCE, ladder tables (14.5 MB) tracked | as stated |
| D9 | Figure archives: keep the two tracked `figures/archived_2026-08-*` on Overleaf; move the untracked `figures/_archived/` (May) into the analyses archive; retire the `**/_archived/**` ignore rule | as stated |
| D10 | Root leftovers (`analysis_setup.md`, `findings.md`, `progress.md`, `task_plan.md`, `outputs/`, non-empty `main N.*` residue, `meascol.*`, `*Notes.bib`) to `_archive/2026-09_reorg/`; delete only the 20 zero-byte `.synctex(busy)` files and caches | as stated |
| D11 | Stamp notes: dated correction blocks (timeless phrasing) rather than rewrites | as stated |
| D12 | Touch the talk repo: wrapper import, remove tracked `.pyc`, add `.gitignore` | yes (three small changes) |
| D13 | Production callable: keep byte-identical with a `DEPRECATED.md` sibling; the permaware folder becomes production | as stated |
| D14 | The two tracked figure orphans (`appendix_3cumulant_fastnc.pdf`, `cumulant_hierarchy.pdf`): `git mv` to `figures/archived_2026-09_orphans/` so `figures/` holds exactly the 17 live PDFs | yes |
| D15 | `revision_plan/` (editor letter, email) and `papers/`: gitignored `docs/private/` and `docs/papers/` | yes |
| D16 | Push the analyses repo to GitHub after Phase 1 (it is 1 commit ahead now and would gain about 45 MB) | operator's call |

## Appendix A: what Phase 0 ran, and where

- Harness: `reorg_2026-09/phase0/harness/p0_regen.py --figure <key> --mode <standalone|rebound|cache|valfig>`
  imports a generator, rebinds its output constants into scratch (never the repo),
  optionally rebinds its FK input, runs `main()`; one figure per process.
  `pdfcmp.py <deployed> <candidate>`: raw md5, `/CreationDate`-stripped md5, raster
  diff at 60 dpi. `p0_numbers.py`: the 34-row table.
- Scratch outputs (not in the repo): `/private/tmp/claude-501/.../scratchpad/p0/`.
- Not run in Phase 0 (by design): the FK fold, the FF multi-z sweeps, any table build,
  the Feynman notebook, `fk_kernel_crosscheck.py`, the BiHalofit n_eff derivation, any
  LaTeX build into the repo root (the baseline build used a scratch copy of the sources).
- Both repos were left exactly as found (verified with `git status` after every run;
  the harness created only ignored `__pycache__` directories).

## Appendix B: environment facts for REPRODUCE.md

| figure rows | interpreter | extra |
|---|---|---|
| 1 | sft-wick env | none |
| 2, 3, 4, 5 (plot), 7, 11, 12, 13, 17 | PyCCL env | pyccl only for a validation print in row 3; camb for row 11 |
| 6, 8, 9, 14 | PyCCL env | `PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src`; row 14 also needs the `scripts.canoes_pipeline` shim until section 5.2 item 9 lands |
| 15, 16 | sft-wick env, Jupyter | sft-wick commit `dd03748` or later |
| sweeps (Phase 3) | sft-wick env | canoes importable through the callables' resolver; `n_jobs <= 6` |

Tools present: pdftoppm/pdfinfo (poppler), qpdf, gs, ImageMagick, latexmk, wolframscript
(licence expired per the brief; SymPy stands in).

## Appendix C: the 18 uncommitted items, grouped for Phase 1

| group | files | proposed message |
|---|---|---|
| A: environment repair and build flags | `equal_time_limber/{equal_time_limber_kappa3_callable.py,_local_cosmo_pk.py,build_equal_time_limber_table.py}` | `fix: resolve the canoes checkout dynamically; add b-model, triples-npz, cache-low and ell-bands options to the vertex build` |
| B: revision-era generator changes | `analysis3/{_plot_style.py,plot_analysis3_nlo_decomposition.py,plot_analysis3_cl_decomposition.py,plot_cl_EB_polarization.py,plot_multiz_kappa_2x5.py}`, `mc_sachs_2pt/{_plot_style.py,fig_xi_channels.py}` | `feat: 2026-08 revision figures (faint FK outside its converged range, 2-12 arcmin band, percent strips, EB parity label, real per-z FF slices, FK markers from mc_fk_complete)` |
| C: regenerated intermediates and outputs | `analysis3/outputs/{multiz_kappa_2pcf_5z.npz,multiz_kappa_xi_cl_2x5.pdf,.png}`, `mc_sachs_2pt/{outputs/appendix_mc_curve.npz,figures/xi_kappa_channels.pdf,.png}` | `data: revision outputs (perm-aware cut15360 FK planes; FK Monte-Carlo markers)` |
| D: new derivation | `mathematica/beam_width_regulator.wl` | `feat: beam-width regulator derivation (xAct, 15 of 15 checks)` |
| not committed | `equal_time_limber/_smoke_equal_time_limber.npz` | archived in Phase 2 |
