# 2PCF / C_corr_op_K_limber_FK

**Observable**: lensing 2-point ξ(γ) at **Order-[0,2]**, **FK sweep**, with the **corr_op**
C propagator and the **equal-time-Limber** κ³ vertex. This is the run that reproduces
the validated **FK κκ ≈ +3.086e-5** baseline (γ=0.5′, z_s=5).

## Wiring (`config_L2.yaml`)
- C propagator: `../../../callables/C_propagator/corr_op/corr_op_C_callable.py` → `C_fn`.
- κ³ K vertex: `../../../callables/kappa3_vertex/equal_time_limber/equal_time_limber_kappa3_callable.py`
  → **`coupling_fn_batch`** (vectorized), `equal_time: true`, **no `already_R_contracted`
  key** (→ false), `coupling_vectorized: true`.
- `expand.orders: [0, 2]`, `sweep.vertex_types: ["FK"]`, `sweep.integrate_over: "all"`.

## equal_time_limber contract (the key difference vs R_contracted)
- **`equal_time: true`** — sft-wick collapses the worldline time to `t_final` (one
  source plane per triple) and uses the table's λ-shell grid as the radial coordinate.
- **`already_R_contracted` absent (false)** — the κ³ table is NOT response-folded;
  **sft-wick applies the response R itself**. The callable's fail-loud guard checks the
  table metadata carries `already_R_contracted=False` AND `equal_time=True`.
- Uses the vectorized `coupling_fn_batch` (`.vectorized=True`) — the only run that
  exercises a vectorized κ³ path.

## Why this is the baseline
The +3.086e-5 FK κκ was measured on this equal-time Limber κ³ (LOW exact-3j ℓ≤60 +
HIGH Born-Limber ℓ>60). It is the comparison target for the windowed R-contracted FK
(`../C_corr_op_K_contracted_FK/`): windowed = full non-Limber, this = χ=λ collapse.

## Geometry
THE +3.086e-5 baseline geometry verbatim (40-pt γ, **t_final=2313.0288751857356**,
n_gauss=24). All other configs inherit this geometry from this run's archived source.

## Run (when ready)
```
SFT_WICK_CONFIG_YAML=<abs path> SFT_WICK_EXPAND_ORDERS=0,2 SFT_WICK_SWEEP_N_GAUSS=24 \
  python ../../../scripts/run_FF_single.py \
    ../../../scripts/inputs/kappa2_grid_lcut0_low_empty.npz \
    ./xi_C_corr_op_K_limber_FK.npz
```
No `SFT_WICK_KAPPA3_NPZ` env — the equal_time_limber callable bundles its combined
LOW+HIGH table (unlike the archived run which summed two NPZs at runtime via env).

## Validation target
κκ(0.5′) ≈ **+3.086e-5**; SPT tree-level cross-check ~30%.

## SUPERSEDED (2026-09-02)

This is the cut1000 sweep of the June 2026 draft, made with the corner-frozen vertex table and the
sorting callable. It is kept byte-identical because SFT-WL-B vendored it (sha256 prefix
`c3eae8cb8515668`) and because its config is the template `run_fk_variant.py` copies. The
manuscript's FK is `../C_corr_op_K_limber_FK_cut15360_permfix/`. The `+3.086e-5` figure above is a
pre-2026-06-09 stamp (before the h^4 and (1+z)^-4 fixes); this file's own kappa-kappa value at 0.5'
is `+4.2877e-6`, and the manuscript's converged value is `+1.9480e-5`. The byte-identical twin
`xi_C_corr_op_K_limber_FK_PROD_REF_2026-06-10.npz` was removed (history at tag pre-reorg-2026-09).
