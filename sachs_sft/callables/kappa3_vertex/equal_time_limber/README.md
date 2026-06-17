# kappa3_vertex / equal_time_limber

**Kind**: non-local κ³ vertex `K` for sft-wick.
**Implementation**: NON-R-contracted equal-time Limber κ³ — LOW exact Wigner-3j
(ℓ≤60) + HIGH flat-sky Born-Limber (ℓ>60), combined. The source of the validated
FK κκ **+3.086e-5** baseline. The SECOND κ³ implementation alongside R_contracted.

## sft-wick contract
- **Callable**: `equal_time_limber_kappa3_callable.py` → `coupling_fn(n_list,
  t_list) -> (3, 3, 3)` (scalar) + `coupling_fn_batch(n_arr, t_arr) ->
  (n_samples, 3, 3, 3)` (`.vectorized=True`).
- **t_list**: three source positions in **λ_project [physical Mpc]**.
- **(3,3,3) tensor**: equal-time κ³, component order **(Φ00, Re Ψ0, Im Ψ0)**
  (Im-Ψ0 block structurally zero).
- **`already_R_contracted = False`** — sft-wick **applies the response R itself**
  (this is the key contract difference vs R_contracted, which is pre-convolved).
- **`equal_time = True`**, **`radial_measure = "lambda"`** (lambda-density).
- Indexed by COSINE-triples (cKDTree, cos γ_ij = n_i·n_j) + λ_shells, NOT
  (L,L,ℓ₃)+λ — intrinsic to the `Kappa3Output` zeta-table format.
- Explicit `TABLE_PATH` (ONE combined LOW+HIGH NPZ, no `SFT_WICK_*` env, no
  runtime LOW/HIGH summation); conventions read from `cosmo_meta`; fail-loud guards.

## Physics: how this differs from R_contracted
- **R_contracted**: windowed squeezed κ³, R-folded at build time, full non-Limber,
  (L,L,ℓ₃) harmonic; sft-wick must NOT re-convolve.
- **equal_time_limber (this)**: equal-time (single source plane per triple), Limber
  HIGH (ℓ>60 δ-collapse) + exact 3-j LOW (ℓ≤60), NOT R-folded; sft-wick applies R.
  `spt_kind="tree_phi"` (the per-leg Poisson A(a) is applied internally; the build
  feeds MATTER P_δ(z=0), not P_φ). Growth: (1+z)^4 bare driving field per leg
  (D⁴ both LOW and HIGH after the 2026-05-30 fix).

## λ-grid dependency (must match R_contracted)
The build reads the L2 λ-grid from an L2 NPZ (`--l2-npz`), so this κ³ shares the
identical 16-node source grid as R_contracted (~397–2327 Mpc). **Build
R_contracted FIRST**, then point `--l2-npz` at its table. Fallback: the archived
deployed L2 NPZ.

## MEASURE (the most error-prone convention here)
`radial_measure="lambda"` — the +3.08e-5 baseline is the LIMBER grid equal_time
(1+z)^4/leg **lambda-density** result. The callable reads `radial_density_measure`
from `cosmo_meta` and **FAILS LOUD** if it is not `lambda` (a chi-density table
would bias the K-vertex by (dχ/dλ)^3 per shell). The archived `kappa3_callable.py`
carried a 2026-05-17 note that applying the (dχ/dλ)^3 Jacobian inflates FK Order-2
25×; that note was for a different consumption path. The lambda-density build IS
the one that reproduces +3.08e-5 through sft-wick's equal-time K-vertex.

## Provenance (reproducible)
Table built by `build_equal_time_limber_table.py` (ports the archived
`analysis_3/build_kappa3_equal_time_limber_l2setup.py`):
```
LOW  = compute_kappa3_zeta_table(triples, λ_h, pk_matter=PCAMBz0, cosmo=FidCosmo,
         ell_max=60, Nmax=128, spt_kind="tree_phi", radial_measure="lambda", ...)
HIGH = compute_kappa3_sigma3_high(triples, λ_h, ell_cut=60, ell_high_max=1000,
         n_ell=96, n_phi=64, radial_measure="lambda", ...)
combined = kappa3_combine_low_high(LOW, HIGH)
```
Cosmo: `FiducialCosmology()` = canoes default = Ω_m=0.3160919980475834, h=0.6711,
n_s=0.97 (== R_contracted / corr_op). **P(k) table: `PCAMBz0.txt`** (the +3.08e-5
baseline was measured with it — a documented difference vs R_contracted/corr_op's
`PCAMB_pyccl_stf_fid_z0.txt`; the equal-time Limber baseline target pins the choice).

## Validation
- The **+3.086e-5 FK κκ baseline** (κκ at γ=0.5'): re-run sft-wick FK on the
  `config_FK_LIMBER_z5.yaml` geometry (40-pt positions_grid, t_final=2313.03,
  n_gauss=24, equal_time) and confirm κκ(0.5')≈+3.086e-5. SPT tree-level
  cross-check ~30%. (The original NEWchi FK *config* provenance was lost; this
  build re-establishes a reproducible source.)

## Files
- `equal_time_limber_kappa3_callable.py` — the sft-wick callable.
- `build_equal_time_limber_table.py` — reproducible build (`--smoke`).
- `_local_cosmo_pk.py` — vendored FiducialCosmology + load_pk_delta (self-contained).
- `equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz` — the combined LOW+HIGH table.
- `README.md` — this file.

## Status (2026-05-31)
Callable + build script + vendored deps done. Table build pointed at the
R_contracted λ-grid. FK +3.08e-5 re-validation is the end-to-end check.
