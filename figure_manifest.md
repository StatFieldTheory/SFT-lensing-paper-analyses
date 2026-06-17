# Figure manifest — STF_lensing paper

Live dependency map for every figure the **current** paper renders. Built by tracing
each `\includegraphics` backward to its generator, inputs, env, and any copy/rename step.

Environments (subagent shells cannot `conda activate`; use the absolute interpreter):
- **PyCCL**: `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`
- **sft-wick**: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`
- Paths below are relative to the repo root `/Users/zzhang/Documents/MyDrafts/STF_lensing`.

## IMPORTANT — the paper renders 16 figures, not 17

A naive `grep includegraphics` reports 17, but `figures/driving_field_cumulants.pdf`
is in a **commented-out** `figure` block (`sections/cosmology.tex:775`, line starts with `%`)
and the PDF does not exist. It is **not** a live dependency. The comment-stripped live set
(`grep includegraphics ... | grep -v '^\s*%'`) is the 16 rows below. The paper builds clean
(`main.pdf`, 30 pages, no missing-figure errors).

---

## Live figures (16)

| # | Paper figure (figures/) | Generator | Generator output → copy step | Env | Key inputs |
|---|---|---|---|---|---|
| 1 | `analysis1_O0_vs_pyccl.pdf` | `sachs_sft/analyses/analysis1/plot_analysis1_O0_vs_pyccl_fkem.py` | writes `analysis1/outputs/analysis1_O0_vs_pyccl_fkem.pdf` → **copied** (suffix `_fkem` dropped) | sft-wick | `sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz`; `analysis1/pyccl_xi_shear_reference_fkem.npz` |
| 2 | `analysis3_NLO_FFFK.pdf` | `sachs_sft/analyses/analysis3/plot_analysis3_nlo_decomposition.py` | writes `analysis3/outputs/analysis3_nlo_O0_FF_FK.pdf` → **copied** | sft-wick | `sftwick_outputs/2PCF/{C_corr_op_O0, C_corr_op_K_limber_FF, C_corr_op_K_limber_FK}/xi_*.npz` |
| 3 | `analysis3_cl_full_vs_O0.pdf` | `sachs_sft/analyses/analysis3/plot_analysis3_cl_decomposition.py` | writes `analysis3/outputs/analysis3_cl_O0_FF_FK.pdf` → **copied** | PyCCL (needs `pyccl` for validation + scipy Wigner-d) | same 3 `xi_*.npz` as #2 |
| 4 | `cl_EB_polarization.pdf` | `sachs_sft/analyses/analysis3/plot_cl_EB_polarization.py` | writes `analysis3/outputs/cl_EB_polarization.pdf` → **copied** (same name) | PyCCL | `sftwick_outputs/2PCF/{C_corr_op_K_limber_FF, C_corr_op_K_limber_FK}/xi_*.npz` |
| 5 | `multiz_kappa_xi_cl.pdf` | **2-stage**: `sachs_sft/analyses/analysis3/compute_multiz_kappa_2pcf.py` then `plot_multiz_kappa_2x5.py` | plot writes `analysis3/outputs/multiz_kappa_xi_cl_2x5.pdf` → **copied** (suffix `_2x5` dropped) | compute=sft-wick, plot=PyCCL | compute reads **external** `talk/assets/figures/_data/multiz_components.npz` + `talk/scripts/run_multiz_components.py` + the 3 live `config_L2.yaml` + corr_op & equal_time_limber callables → `analysis3/outputs/multiz_kappa_2pcf_5z.npz`; plot reads that npz |
| 6 | `appendix_mc_workflow.pdf` | `sachs_sft/analyses/mc_sachs_2pt/fig_xi_channels.py` | writes `mc_sachs_2pt/figures/xi_kappa_channels.pdf` (+ cache `mc_sachs_2pt/outputs/appendix_mc_curve.npz`) → **copied** | PyCCL | `sftwick_outputs/2PCF/<name>/xi_<name>.npz`; local modules `driver_stats.py`, `fk_analytic.py`, `sachs_mc_core.py`, `background.py`, `_plot_style.py`; `--replot` reads only the cache npz |
| 7 | `zeta_driving_field_slices_draft.pdf` | `sachs_sft/callables/kappa3_vertex/equal_time_limber/plot_zeta_figures.py` | writes **directly** to repo `figures/zeta_driving_field_slices_draft.{pdf,png}` (hardcoded path — **no copy; overwrites the paper file if rerun**) | numpy+matplotlib (either env) | `equal_time_limber/ell_band_decomp_results.npz` |
| 8 | `corr_operator_slices_draft.pdf` | `SFT-lensing-paper-analyses/figures/make_response_corr_operator_figures.py` | writes **directly** to repo `figures/corr_operator_slices_draft.{pdf,png}` (no copy; overwrites if rerun) | numpy+scipy (PyCCL fine) | imports `corr_op_C_callable` (+ its `TABLE_PATH`) and `D_callable` from `sachs_sft/callables/C_propagator/corr_op/` (added to `sys.path`) |
| 9 | `response_operator_draft.pdf` | `SFT-lensing-paper-analyses/figures/make_response_corr_operator_figures.py` (same script as #8) | writes **directly** to repo `figures/response_operator_draft.{pdf,png}` | numpy+scipy | same corr_op callable as #8 |
| 10 | `screen_basis_spin2_pattern.pdf` | `SFT-lensing-paper-analyses/figures/make_screen_basis_spin2_figure.py` | writes **directly** to repo `figures/screen_basis_spin2_pattern.{pdf,png}` | python+matplotlib (no data) | none (procedural) |
| 11 | `driving_field_portrait.pdf` | `SFT-lensing-paper-analyses/figures/plot_driving_field_portrait.py` | writes **directly** to repo `figures/driving_field_portrait.{pdf,png}` (default `--out-stem`) | PyCCL + **camb** | none external; `--seed 42`, n_pix=192, box=300 Mpc, z=0.7 |
| 12 | `selection_rule_schematic.pdf` | `SFT-lensing-paper-analyses/figures/make_selection_rule_schematic.py` | writes **directly** to repo `figures/selection_rule_schematic.{pdf,png}` | python+matplotlib | none |
| 13 | `null_geodesic_congruence.pdf` | `SFT-lensing-paper-analyses/figures/make_null_congruence_schematic.py` | writes **directly** to repo `figures/null_geodesic_congruence.{pdf,png}` | python+matplotlib | none |
| 14 | `cosmology_coord_maps.pdf` | `SFT-lensing-paper-analyses/_archive/canoes_pipeline/scripts/plot_coord_maps.py` (**generator only exists in archive**) | writes repo `figures/cosmology_coord_maps.{pdf,png}` | python; needs **external `/Users/zzhang/projects/canoes/src`** + canoes `_plot_style` | `FiducialCosmology()` from canoes; z-grid linspace(0,5,401) |
| 15 | `FeynDiag_2pt_order2.pdf` | **NO GENERATOR FOUND** in repo | — | — | hand-made vector art or deleted generator; `talk/scripts/make_feyndiag_dark.py` makes differently-named talk SVGs |
| 16 | `FeynDiag_3pt_order1.pdf` | **NO GENERATOR FOUND** in repo | — | — | same as #15 |

### One-line regeneration commands (cheap, safe-to-rerun figures)

```bash
# repo root
PY_CCL=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
SFTW=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

# Schematics / glyphs (SFT-lensing-paper-analyses/figures/ — write directly into figures/; use --out-stem or a scratch dir to avoid clobbering):
python SFT-lensing-paper-analyses/figures/make_screen_basis_spin2_figure.py
python SFT-lensing-paper-analyses/figures/make_selection_rule_schematic.py
python SFT-lensing-paper-analyses/figures/make_null_congruence_schematic.py
$PY_CCL  SFT-lensing-paper-analyses/figures/plot_driving_field_portrait.py --seed 42
$PY_CCL  SFT-lensing-paper-analyses/figures/make_response_corr_operator_figures.py     # corr_operator_slices_draft + response_operator_draft

# Analysis figures (run inside the analysis dir):
( cd SFT-lensing-paper-analyses/sachs_sft/analyses/analysis1 && $SFTW plot_analysis1_O0_vs_pyccl_fkem.py )
( cd SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3 && $SFTW plot_analysis3_nlo_decomposition.py )
( cd SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3 && $PY_CCL plot_analysis3_cl_decomposition.py )
( cd SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3 && $PY_CCL plot_cl_EB_polarization.py )
( cd SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt && $PY_CCL fig_xi_channels.py --replot )   # replot from cache

# zeta slices (callable; writes straight into figures/):
$SFTW SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber/plot_zeta_figures.py

# multi-z (2-stage; compute is expensive):
( cd SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3 && $SFTW compute_multiz_kappa_2pcf.py && $PY_CCL plot_multiz_kappa_2x5.py )
```

### Live data inputs (must NOT be archived)

- `sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz`
- `sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FF/xi_C_corr_op_K_limber_FF.npz`
- `sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK.npz`
- the three matching `config_L2.yaml` (read by the multi-z compute step)
- `sachs_sft/callables/C_propagator/corr_op/` (corr_op_C_callable.py, D_callable.py, `corr_op_table_*.npz` 5.3M)
- `sachs_sft/callables/kappa3_vertex/equal_time_limber/` (equal_time_limber_kappa3_callable.py, `*.npz` 1.3M, plot_zeta_figures.py, ell_band_decomp_results.npz)
- `sachs_sft/scripts/{D_callable.py, kappa2_callable.py, c_fast_gl.py, spin_rotation.py}` + `scripts/inputs/{F_tensor.npy, kappa2_grid_lcut0_low_empty.npz}`
- external: `talk/assets/figures/_data/multiz_components.npz`, `talk/scripts/run_multiz_components.py`

### Not a live dependency (documented, not built)

- `figures/driving_field_cumulants.pdf` — commented out (`cosmology.tex:775`); file absent; a canoes data figure dropped from the paper.
- `figures/cumulant_hierarchy.pdf` — exists + generated by `SFT-lensing-paper-analyses/figures/make_cumulant_hierarchy.py`, but **not referenced by any `.tex`** (orphan / alternative schematic). Generator kept.
