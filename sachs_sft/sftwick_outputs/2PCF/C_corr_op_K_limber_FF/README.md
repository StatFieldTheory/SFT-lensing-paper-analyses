# 2PCF / C_corr_op_K_limber_FF

**Observable**: lensing 2-point ξ(γ) at **Order-[0,2]**, **FF sweep**, with the **corr_op**
C propagator and the **equal-time-Limber** κ³ vertex declared (but excluded from FF
diagrams). The F-only counterpart of [`../C_corr_op_K_limber_FK/`](../C_corr_op_K_limber_FK/).

## Wiring (`config_L2.yaml`)
Byte-identical to `C_corr_op_K_limber_FK` **except** `sweep.vertex_types: ["F"]`
(vs `["FK"]`) and the folder-local cache/output paths.
- C propagator: `../../../callables/C_propagator/corr_op/corr_op_C_callable.py` → `C_fn`.
- κ³ K vertex (declared): `../../../callables/kappa3_vertex/equal_time_limber/equal_time_limber_kappa3_callable.py`
  → `coupling_fn_batch` (vectorized), `equal_time: true`, no `already_R_contracted`
  key (→ false), `coupling_vectorized: true`.
- `expand.orders: [0, 2]`, `sweep.vertex_types: ["F"]`, `sweep.integrate_over: "all"`.

> **Token note (verified against sft-wick).** `vertex_types` entries must match the
> diagram's `_vertex_type_label`, which is the **sorted-unique** coupling names joined
> (`expansion.py:422`): a both-F diagram → `"F"`, an F×K cross diagram → `"FK"`. So the
> both-F selector is **`["F"]`, NOT `["FF"]`** — `["FF"]` matches no diagram, yielding
> an empty sweep and a `KeyError: 'x'` in `sweep.totals()`. (The sibling
> `C_corr_op_K_contracted_FF` still carries the buggy `["FF"]` and will fail the same
> way until corrected to `["F"]`.)

## FF vs FK (what this run isolates)
`vertex_types` selects the Order-2 diagram by its vertex-composition label:
- **F** (this run): only the *both-F* (local-cubic × local-cubic) post-Born diagrams.
- **FK** (sibling): only the F × K (local-cubic × κ³) cross diagrams.

The K (κ³) vertex is **declared** in `nonlocal_vertices` (so its callable is imported),
but it does **not** enter any both-F diagram. Consequently the result is **independent
of the κ³ table/wiring** — `C_corr_op_K_limber_FF` should reproduce
[`../C_corr_op_K_contracted_FF/`](../C_corr_op_K_contracted_FF/) κκ to numerical
agreement (same corr_op C, same F-only diagrams) once that sibling's `["FF"]`→`["F"]`
bug is fixed. Completes the {contracted, limber} × {F, FK} matrix.

## Expectations
The `vertex_types: ["F"]` filter keeps **only the Order-2 both-F diagram**; the
Order-0 Born term (zero vertices, label `""`) does **not** match and is excluded — so
the output carries **`order = 2` rows only** (no Order-0). To assemble the full
Order-[0,2] κκ, add `../C_corr_op_O0/` (Order-0) + this (Order-2 F) + the FK Order-2.
- **Order-2 (F)** κκ: the both-F post-Born correction; should equal the contracted-FF
  Order-2 (κ³-wiring-independent).

## Geometry
The +3.086e-5-baseline geometry verbatim, inherited from `C_corr_op_K_limber_FK`
(40-pt γ from 0.5′ to 5000′, **t_final=2313.0288751857356**, 6 component pairs,
n_gauss=24).

## Run (when ready)
The κ³ callable still loads its bundled equal_time_limber table at import (even
though FF excludes K), so that table must exist. Cap the joblib worker pool to keep
memory bounded (the config sets `n_jobs: -1`):
```
LOKY_MAX_CPU_COUNT=4 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
CANOES_SUPPRESS_METAL_WARNING=1 \
SFT_WICK_CONFIG_YAML=<abs path to this config_L2.yaml> \
  python ../../../scripts/run_FF_single.py \
    ../../../scripts/inputs/kappa2_grid_lcut0_low_empty.npz \
    ./xi_C_corr_op_K_limber_FF.npz
```
Do NOT set `SFT_WICK_EXPAND_ORDERS` / `SFT_WICK_SWEEP_N_GAUSS` — let the config's own
`orders: [0, 2]` and `n_gauss: 24` stand.

## Validation target
Order-2 (F) κκ should match `C_corr_op_K_contracted_FF`'s Order-2 (F) κκ once that
sibling is corrected to `vertex_types: ["F"]` (κ³-wiring-independent). Output is
Order-2 only (see Expectations).

## Status (2026-09-02)

Live production sweep of the manuscript (unchanged since 2026-06-01). The `+3.086e-5 baseline geometry` wording refers to the 40-point gamma grid, t_final = 2313.029 and n_gauss = 24 shared by every run; the number itself is a pre-2026-06-09 stamp of the FK kappa-kappa value. `PRODUCTION` marks this folder for the runner's overwrite guard.
