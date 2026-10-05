# sftwick_outputs/

Each folder here is **one concrete sft-wick run**: a self-contained L2 `config_L2.yaml`
plus (after running) its output NPZ and joblib caches. Configs are wired to the
**sachs_sft callables** (`../callables/...`, absolute paths) — lensing correlation
observables produced by feeding those callables through sft-wick's L2 expansion.

**2PCF vs 3PCF is decided by the OBSERVABLE point count** (`expand.observable`), NOT
by which vertices appear. The 5 runs below all have `observable: [phi_a(x), phi_b(y)]`
= **2-point**, so they all live under `2PCF/`. The "FF"/"FK" in a run name refers to
the sft-wick **vertex types** (F = local cubic, K = non-local κ³) that enter the
**Order-2 post-Born corrections** to the 2-point function — a κ³ vertex inside an
Order-2 diagram still contributes to ξ(γ), it does not make the observable 3-point.
`3PCF/` is reserved for future genuine 3-point runs (`observable: [phi(x),phi(y),phi(z)]`).

```
sftwick_outputs/
├── 2PCF/                            # observable = [phi_a(x), phi_b(y)] (2-point ξ(γ))
│   ├── C_corr_op_O0/                # corr_op C, Order-0  (≙ old analysis_1b corr_op L2 O0)
│   ├── C_limber_O0/                 # Born-Limber C, Order-0
│   ├── C_corr_op_K_contracted_FF/   # corr_op C + R-contracted κ³, O[0,2], F-vertex sweep
│   ├── C_corr_op_K_contracted_FK/   # corr_op C + R-contracted κ³, O[0,2], F+K sweep
│   └── C_corr_op_K_limber_FK/       # corr_op C + equal-time-Limber κ³, O[0,2], F+K sweep
└── 3PCF/                            # (reserved) genuine 3-point: observable [phi(x),phi(y),phi(z)]
```

## Historical templates and current candidate wiring

The following May/June configuration descriptions are historical templates, not
current accepted products or safe replay commands. Current figures resolve the
hash-checked selections in `reproduce/active_products.json`. The accepted T-001 saved-product selection is now published separately under
`analyses/r1_sft061/t001_accepted/`; these historical templates retain their
original input bindings and are not corrected-fold replay configurations.

Reviewed corr_op FF/FK candidate configurations bind `C_fn_batch` with
`c_closed_form_vectorized: true`; the historical scalar-only wiring below does
not describe those configurations. Use each frozen candidate's actual config
and provenance, not this prose, to establish its source endpoint and callable.
Do not run `gen_sftwick_configs.py` in place: it can overwrite hash-pinned
historical configs and targets absent folders. Prepare new configs in a fresh,
reviewed output directory; keep historical prepared/run records unchanged.

## How a historical run was wired (the L2 config contract)

Every `config_L2.yaml` has 5 top-level blocks (`system`, `expand`, `propagators`,
`sweep`, `output`). The sachs_sft-specific choices:

- **C propagator** → `propagators.c_closed_form_module` = an absolute path to a
  sachs_sft C callable (`corr_op` or `limber`), `c_closed_form_attr: C_fn`,
  `c_closed_form_only: true`. **`c_closed_form_vectorized: false`** — the sachs_sft
  `C_fn` wrappers are scalar (they do not re-export the underlying `.batch`), so a
  vectorized array-t call would crash; the C propagator is cached once anyway.
- **κ³ K vertex** (the Order-2 runs) → `system.nonlocal_vertices[0].coupling_module` =
  an absolute path to a sachs_sft κ³ callable. The load-bearing flags:
  - R_contracted: `already_R_contracted: true`, `equal_time: false`,
    `coupling_attr: coupling_fn` (scalar). sft-wick runs the full worldline integral
    and must NOT re-fold the response (it is pre-folded in the table).
  - equal_time_limber: `equal_time: true`, NO `already_R_contracted` key (→ false),
    `coupling_attr: coupling_fn_batch` + `coupling_vectorized: true` (this callable
    exposes a vectorized batch). sft-wick collapses time to `t_final` and folds R itself.
- **Response R / growth** → `system.linear.R_time_module` = `../scripts/D_callable.py`
  (`response_R`, the Jacobi D(λ)=a·χ growth). **kappa2 source** → `../scripts/kappa2_callable.py`.
  **F vertex** (Order-2) → `../scripts/inputs/F_tensor.npy`.
- **`sweep.integrate_over: "all"`** in every config (the source-plane marginalisation
  is always performed).
- **`sweep.vertex_types`**: `["F"]` (FF / local-cubic post-Born) vs `["FK"]`
  (F+K cross diagrams). `expand.orders`: `[0]` (Born only) vs `[0,2]` (+ post-Born).

## Provenance: the +3.08e-5 baseline geometry

All 5 configs are generated from the verified `+3.086e-5` FK κκ baseline geometry
(archived `config_FK_LIMBER_z5.yaml`) by `../scripts/gen_sftwick_configs.py`. The
**geometry is byte-identical across all 5** — 40-pt log-spaced γ grid (0.5′–5000′),
`t_final = 2313.0288751857356` Mpc (z_s=5), 6 `component_pairs`, `n_gauss = 24`.
The only differences between configs are the intended wiring (which C, which K,
FF vs FK). This "bit-equivalence by construction" guarantees a config diff is pure
physics intent, never an accidental geometry drift.

## Reproducing the configs

```
python ../scripts/gen_sftwick_configs.py     # rewrites all 5 config_L2.yaml
```

## Status (2026-05-31)
Configs authored + reviewed (1 hard bug fixed: `c_closed_form_vectorized` →
`false`; `integrate_over` → `"all"` per requirement; equal_time_limber uses the
vectorized batch). **No runs yet** — these are the wiring, ready to execute.
