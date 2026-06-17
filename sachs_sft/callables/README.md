# callables/ — sft-wick callable objects (one folder per implementation)

Two kinds, two implementations each. Each implementation is a **self-contained
folder**: the callable module + its product(s)/build-recipe + convention notes +
validation all live together. A sft-wick config references exactly one folder per
kind.

## C_propagator/ — correlation propagator `C = ⟨Φ Φ⟩`

| folder | implementation | source (in archive) | status |
|---|---|---|---|
| `corr_op/` | canoes `corr_op` 2-pt `C_ℓ(λ',λ'')` table → propagator | `_archive/canoes_pipeline/inputs/l2_callables_2026_05_26/corr_op_*` + `corr_op_l2_callable.py` | validated |
| `limber/` | canoes Limber-version C → propagator | `_archive/canoes_pipeline/...` (Limber C path) | validated |

## kappa3_vertex/ — non-local κ³ vertex `K`

| folder | implementation | source (in archive) | status |
|---|---|---|---|
| `R_contracted/` | R-contracted-already (windowed) κ³ | `_archive/.../windowed_squeezed_kappa3_l2_callable.py` + table | validated |
| `equal_time_limber/` | (non-contracted) equal-time Limber κ³ | `_archive/.../build_kappa3_equal_time_limber_*` + zeta/sigma3 inputs | validated |

## Per-folder README must state

1. **sft-wick contract**: coordinate (λ_project physical Mpc?), measure (chi vs
   lambda density), units, `already_R_contracted` (true/false), component order,
   `t_final`/horizon, the `coupling_fn(n_list, t_list)` or `C(...)` signature.
2. **Provenance**: which archived build produced the bundled product; cosmo;
   growth; atom scope; path (A/C); the exact build command.
3. **Validation**: what it was checked against and the result.
