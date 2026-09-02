# Figure manifest

Every figure the manuscript renders, traced to its generator, inputs and environment.
Rewritten 2026-09-02 after the reproducibility reorganisation; each row was verified by
regenerating the figure and comparing it byte-for-byte with the deployed PDF.

**All 17 reproduce with no manual rebinding.** Run them together with

```bash
python reproduce/regen_figure.py --all && python reproduce/check_figures.py
```

The warning that used to head this file, that running a generator standalone silently
produced an FK channel a factor 4.5 too low, no longer applies: the generators default to
the manuscript's converged sweep. See `CHANGELOG_REORG_2026-09.md`.

Environments (subagent shells cannot `conda activate`; use the absolute interpreter):

* **PyCCL**: `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`
* **sft-wick**: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`
* `+canoes` marks a run that imports canoes; the callables resolve it themselves through
  `CANOES_ROOT`, an importable `canoes`, or two known locations.

Paths below are relative to this package unless they start with `figures/`, which is the
manuscript's figure folder one level up. The environment column is the one that
reproduces the deployed bytes: the two environments agree on every number but differ at
the 1 pt bounding-box level.

## The 17 figures

| # | paper figure | generator | inputs | env | time |
|---|---|---|---|---|---|
| 1 | `analysis1_O0_vs_pyccl.pdf` | `sachs_sft/analyses/analysis1/plot_analysis1_O0_vs_pyccl_fkem.py` | `sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz`; `analysis1/pyccl_xi_shear_reference_fkem.npz` | sft-wick | 2 s |
| 2 | `analysis3_NLO_FFFK.pdf` | `sachs_sft/analyses/analysis3/plot_analysis3_nlo_decomposition.py` | the O0, FF and `C_corr_op_K_limber_FK_cut15360_permfix` sweeps | PyCCL | 2 s |
| 3 | `analysis3_cl_full_vs_O0.pdf` | `analysis3/plot_analysis3_cl_decomposition.py` | the same three sweeps; pyccl only for a console validation print | PyCCL | 4 s |
| 4 | `cl_EB_polarization.pdf` | `analysis3/plot_cl_EB_polarization.py` | the FF and FK sweeps, through the row-3 module's curved-sky transform | PyCCL | 2 s |
| 5 | `multiz_kappa_xi_cl.pdf` | `analysis3/plot_multiz_kappa_2x5.py` | `analysis3/outputs/multiz_kappa_2pcf_5z.npz` and `sftwick_outputs/2PCF/multiz/multiz_ff_real_all5.npz` | PyCCL | 3 s |
| 6 | `appendix_mc_workflow.pdf` | `mc_sachs_2pt/fig_xi_channels.py --from-cache` | `mc_sachs_2pt/outputs/appendix_mc_curve.npz`; FK markers from `analyses/mc_fk_complete/_markers_pooled.npz` | PyCCL +canoes | 4 s |
| 7 | `zeta_driving_field_slices_draft.pdf` | `callables/kappa3_vertex/equal_time_limber/plot_zeta_figures.py` | `equal_time_limber/outputs/zeta_bands_cut15360_gmax85.npz` | PyCCL | 1 s |
| 8 | `corr_operator_slices_draft.pdf` | `figures/make_response_corr_operator_figures.py` | the corr_op table and `sachs_sft/scripts/D_callable.py` | PyCCL +canoes | 65 s |
| 9 | `response_operator_draft.pdf` | the same script | the same | PyCCL +canoes | (same run) |
| 10 | `screen_basis_spin2_pattern.pdf` | `figures/make_screen_basis_spin2_figure.py` | none (procedural) | PyCCL | 1 s |
| 11 | `driving_field_portrait.pdf` | `figures/plot_driving_field_portrait.py --seed 42` | CAMB linear P(k) computed at run time | PyCCL | 3 s |
| 12 | `selection_rule_schematic.pdf` | `figures/make_selection_rule_schematic.py --seed 42` | none | PyCCL | 1 s |
| 13 | `null_geodesic_congruence.pdf` | `figures/make_null_congruence_schematic.py` | none | PyCCL | 1 s |
| 14 | `cosmology_coord_maps.pdf` | `figures/make_cosmology_coord_maps.py` | none; the two canoes constants and the background helpers are vendored in the script | PyCCL | 1 s |
| 15 | `FeynDiag_2pt_order2.pdf` | `figures/make_feyndiag.py` | the sft-wick action and diagram renderer | sft-wick | 1 s |
| 16 | `FeynDiag_3pt_order1.pdf` | the same script | the same | sft-wick | (same run) |
| 17 | `fk_cutoff_convergence.pdf` | `callables/kappa3_vertex/rebuild/fig_cutoff_paper.py` | `sftwick_outputs/2PCF/cutoff_ladder/` (ten folds) and the O0 sweep | PyCCL | 1 s |

Each generator writes into its own `outputs/` folder. `reproduce/deploy.py` is the only
step that copies into the manuscript, and it records md5 sums in
`reproduce/deployed_md5.txt`.

## Live data inputs

These must not be archived. Provenance for each is in the README of the folder that holds
it, and the cost to remake it is in REPRODUCE.md section 4.

* `sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz`
* `.../C_corr_op_K_limber_FF/xi_C_corr_op_K_limber_FF.npz`
* `.../C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz`
* `.../cutoff_ladder/` (ten folds) and `.../multiz/` (the source-redshift sweeps)
* `sachs_sft/callables/C_propagator/corr_op/` (the callable and its 5.5 MB table)
* `sachs_sft/callables/kappa3_vertex/equal_time_limber_cut15360_permaware/` (the production
  vertex: callable, `table_permclosed.npz`, tests)
* `sachs_sft/callables/kappa3_vertex/equal_time_limber/` (the build script, the ell-band
  decomposition, the figure-7 generator, `outputs/zeta_bands_cut15360*.npz`, the vendored
  L2 lambda grid, and the superseded cut1000 table kept for provenance)
* `sachs_sft/callables/kappa3_vertex/rebuild/products/` (the ladder tables and the triple
  grids) and `rebuild/logs/`
* `sachs_sft/scripts/{D_callable,kappa2_callable,c_fast_gl,spin_rotation}.py` and
  `scripts/inputs/`
* `sachs_sft/analyses/analysis3/{outputs/multiz_kappa_2pcf_5z.npz, inputs/multiz_components_talk_2026-06-10.npz}`
* `sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz`
* `sachs_sft/analyses/mc_fk_complete/_markers_pooled.npz`

## Not rendered by the manuscript

* `figures/make_cumulant_hierarchy.py` produces an alternative schematic that no `.tex`
  references. The generator is kept; the PDF was removed from the manuscript's `figures/`
  (recoverable at tag `pre-reorg-2026-09`).
* `figures/appendix_3cumulant_fastnc.pdf` was referenced only by an archived section and
  was removed likewise; its analysis is `sachs_sft/analyses/shear3pcf_fastnc/`, labelled
  `BACKUP_NOT_IN_PAPER.md`.
* `driving_field_cumulants.pdf` is named only inside a commented-out block and the file
  does not exist.
