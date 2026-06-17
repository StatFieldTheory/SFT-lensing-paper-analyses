# C_propagator / corr_op

**Kind**: correlation propagator `C = ⟨Φ Φ⟩` for sft-wick.
**Implementation**: canoes `corr_op` 2-point `C_ℓ(λ', λ'')` table → propagator callable.

## sft-wick contract
- **Callable**: `corr_op_C_callable.py` → `C_fn(n1, t1, n2, t2) -> (3, 3)`.
- **t1, t2**: source positions in **λ_project [physical Mpc]** (`QUERY_COORDINATE
  = "lambda_project_mpc"`).
- **(3,3) tensor**: ⟨Φ Φ⟩, component order **(Φ00, Re Ψ0, Im Ψ0)**.
- **units**: physical Mpc; `output_h_power = 0.0` (no h-rescale; value self-consistent).
- The callable loads the table by **explicit path** (no `SFT_WICK_*` env) and reads
  all conventions **from the table** (`config.source_coord`, `config.cosmology`),
  never hand-typed.

## Coordinate (the chi-vs-lambda clarification)
corr_op's natural source coordinate is the **affine λ**, and the table IS
**λ_project-indexed** (`source_coord = lambda_project_mpc`; `lambda_grid` =
406..2328 Mpc = λ_project, 21 nodes). The inner `chi_min=50` / `n_chi_per_pair=256`
are only the **line-of-sight integration grid** for `C_ℓ(λ',λ'')`, NOT the table's
index. (The old `corr_op_l2_callable.py` mislabeled this via a triple-assigned
`NATIVE_TABLE_COORDINATE` last-winning `chi_mpc`, contradicting the table; the math
was still correct because it read `source_coord` from the table. This clean callable
reads from the table config and adds a fail-loud guard, so the contradiction cannot
recur.)

## Provenance (reproducible)
Table built by `build_corr_op_table.py` (the EXACT validated setup):
```
python -m canoes.sachs.sft_input.corr_op.build
  --pk examples/data/PCAMB_pyccl_stf_fid_z0.txt
  --omega-m 0.3160919980475834 --h 0.6711 --n-s 0.97
  --ell-max 5000 --dense-threshold 500 --n-log-per-decade 16
  --lambda-grid zs5-21pt          # → source_coord=lambda_project_mpc
  --chi-min 50.0 --n-chi-per-pair 256
  --accuracy default              # FFTlog Nmax=256
  --ell-hi-threshold 800          # S3: Limber ℓ>800, t-form ℓ≤800
  --output corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz
```
(`ell_taper_frac=0.2`, `cl_tensor_convention=bare` are CorrOpTable2DConfig defaults.)
Cosmo: Ω_m=0.3160919980475834, h=0.6711, n_s=0.97 (PyCCL STF fiducial).
Requires the canoes CLI `%`-escape fix (commit `0e24fd6`).

## Validation
- This setup is the previously-validated corr_op table (PyCCL cross-check 5-10% on
  kk / xi_+ Order-0; the +3.08e-5 FK κκ baseline derives from this C path).
- Rebuilt 2026-05-30 with the identical setup; **cross-checked vs the archived
  `_archive/.../corr_op_prod_ellmax5000_zs5_21pt_stf_fid_omega03161_h06711.npz`:
  max abs diff across cl_pp/cl_px/cl_xx = 2.9e-15 (machine precision) → PASS,
  bit-for-bit reproduction.** (The naive max *rel* diff is ~4e-2 only because of
  division by near-zero off-diagonal cells; the absolute agreement is exact.)
  Grids (lambda_grid 21, ell_grid 515) match identically.

## Files
- `corr_op_C_callable.py` — the sft-wick callable (`C_fn`).
- `build_corr_op_table.py` — reproducible build (`--smoke` for a 30s wiring check).
- `corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz` — the table (built here).
- `corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.meta.json` — build sidecar.
- `README.md` — this file.

## Status (2026-05-30)
DONE. Full production table built + cross-check PASS (abs diff 2.9e-15 vs archive);
callable validated (`C_fn → (3,3)` finite, source_coord guard passes).
