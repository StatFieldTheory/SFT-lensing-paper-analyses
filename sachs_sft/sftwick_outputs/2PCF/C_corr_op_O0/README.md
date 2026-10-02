# 2PCF / C_corr_op_O0

> **Result arrays archived on 2026-10-02.** The arrays of this run are superseded and were moved to
> `_archive/cleanup_2026-10-02/` (index: `reproduce/archived_inputs.json`). The configuration stays here
> because the referee-revision drivers in `sachs_sft/analyses/r1_sft061/` use it as their template. The
> manuscript's products are the ones selected by `reproduce/active_products.json`. The text below
> describes the run as it stood and is kept as a record.

**Observable**: lensing 2-point correlation `ξ(γ)` (κκ, κγ±, γ±γ±) at **Order-0**
(Born/tree level), using the **corr_op** correlation propagator C = ⟨Φ Φ⟩.
This is the sachs_sft re-creation of the old `analysis_1b` corr_op-L2 Order-0 run.

## Wiring (`config_L2.yaml`)
- C propagator: `../../../callables/C_propagator/corr_op/corr_op_C_callable.py` → `C_fn`
  (`c_closed_form_only: true`, `c_closed_form_vectorized: false`).
- NO κ³ K vertex (`system.nonlocal_vertices` absent).
- `expand.orders: [0]`, `sweep.vertex_types: ["F"]` (F declared but Order-0 → no O2 diagrams).
- `sweep.integrate_over: "all"`.
- Response: `../../../scripts/D_callable.py:response_R`; kappa2: `../../../scripts/kappa2_callable.py`.

## Geometry
40-pt γ grid 0.5′–5000′, `t_final = 2313.029` Mpc (z_s=5), 6 component_pairs, `n_gauss = 24`
(the +3.086e-5 baseline geometry; identical across all sftwick_outputs runs).

## Run (when ready)
```
SFT_WICK_CONFIG_YAML=<abs path to this config_L2.yaml> \
  python ../../../scripts/run_FF_single.py \
    ../../../scripts/inputs/kappa2_grid_lcut0_low_empty.npz \
    ./xi_C_corr_op_O0.npz
```
(Output + caches land in this folder.)

## Validation target
Order-0 κκ / ξ± should match PyCCL at 5-10% (the corr_op C path's validated regime).

## Status (2026-09-02)

Live production sweep of the manuscript (unchanged since 2026-06-01). The `+3.086e-5 baseline geometry` wording refers to the 40-point gamma grid, t_final = 2313.029 and n_gauss = 24 shared by every run; the number itself is a pre-2026-06-09 stamp of the FK kappa-kappa value. `PRODUCTION` marks this folder for the runner's overwrite guard.
