# Sachs dynamics and SFT-Wick calculations

This directory contains optical callables, SFT-Wick run configurations, saved
correlation functions and analysis drivers for the weak-lensing paper. Start
with the repository's [reproduction guide](../REPRODUCE.md) for the supported
figure and numerical-summary workflow.

## Current revision

[../reproduce/active_products.json](../reproduce/active_products.json) selects
the current R1 products and their hashes. Revision calculations live under
`analyses/r1_sft061/`, with the aligned Order-0/PyCCL calculation under
`analyses/analysis1/r1_aligned/`. Historical run folders remain available for
comparison and as dependencies of the replay drivers.

The optical dynamics use the local quadratic Sachs couplings together with
two kinds of callable input: the two-point correlation propagator `C` and the
non-local three-point driving-field vertex `K`. Numerical configurations select
specific callable implementations and data by path. The figure manifest selects
the saved products read by the plotting tools; it does not replace those run
configurations.

## Layout

| Path | Role |
|---|---|
| `callables/C_propagator/corr_op/` | Correlation-propagator implementation and stored table |
| `callables/kappa3_vertex/equal_time_limber_cut15360_permaware/` | Permutation-aware vertex callable; its historical table is archived |
| `callables/kappa3_vertex/equal_time_limber/` | Historical vertex inputs, build helpers and driving-field slice plotter |
| `callables/kappa3_vertex/rebuild/` | Vertex assembly, variant folding, cutoff studies and saved tables |
| `callables/kappa3_vertex/b_model/` | Matter-bispectrum models |
| `analyses/r1_sft061/` | Corrected-input folds, source mapping, cutoff ladder, Monte Carlo controls and figure-product assembly |
| `analyses/analysis1/r1_aligned/` | Continuous source-response Order-0 calculation and PyCCL benchmark |
| `analyses/analysis3/` | Two-point correlation and angular-spectrum figure generators |
| `analyses/mc_sachs_2pt/`, `analyses/mc_fk_complete/` | Monte Carlo implementations and historical controls |
| `sftwick_outputs/2PCF/` | Historical per-run configurations; their arrays are archived |
| `scripts/` | Shared response, covariance, geometry and sweep utilities |

A folder's name or `PRODUCTION` marker alone does not identify the current
figure inputs. Use the active manifest and its input provenance, especially when
comparing historical and revised FK amplitudes. Record the source plane,
multipole cutoff and input/runtime versions with every numerical result.

## Recomputing products

The revision drivers separate planning, preparation and execution where
supported; saved `prepared.json` records identify the accepted settings.
Some paths and guarded runtime hashes refer to the original workstation.
Inspect those dependencies before moving a calculation, and retain original
records alongside new manifests.

Use fresh run directories and keep accepted outputs unchanged. Run folds and
vertex builds sequentially; worker count is constrained by memory use.
Regenerating the figures from saved products is much cheaper than repeating the
folds, vertex integrations and Monte Carlo calculations. Neither operation
modifies the manuscript unless the separate deployment tool is invoked.
