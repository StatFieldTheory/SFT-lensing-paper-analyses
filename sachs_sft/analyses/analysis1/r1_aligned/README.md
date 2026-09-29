# Referee R1: aligned Order-0 comparison

**Last updated: 2026-09-28**

This directory contains the calculation for revising the manuscript's Figure 10.
The current method integrates the continuous outer source response before
computing the Gaussian angular covariance. It evaluates the manuscript's
effective affine kernel (`eq: K affine kernel`, `eq: xi kappa folded`) without
interpolating the old 21-source C table. No amplitude is fitted.

The numerical comparison is accepted at the stated resolution, as recorded in
[source_folded_final_validation.json](source_folded_final_validation.json).
The largest withheld-grid observable change is 0.01616% of its channel peak,
below the 0.1% target. This acceptance does not certify manuscript deployment.
Existing production Order-0, FF, and FK products are separate from these
isolated calculations.

## Shared inputs and conventions

| Setting | Value |
|---|---|
| Source plane | z = 5, chi = 7971.183639083942 Mpc, lambda = 2317.9695958203943 Mpc |
| Cosmology | Omega_m = 0.3160919980475834, h = 0.6711, n_s = 0.97 |
| Power and growth | Shared CAMB z = 0 linear matter table and canoes growth/lightcone mapping |
| Physical radial support | 50 Mpc through the source distance |
| Angular band | Every integer ell from 2 through 5000 in the final spectra |
| Numerical split | Non-Limber through ell = 800 inclusive, Limber above 800 |
| High-ell power argument | (ell + 0.5) / chi |
| Angular taper | Cosine taper over the final 20% of the band, starting at ell = 4000.4 |
| Angular samples | Identical 40-sample grid from the production configuration, also saved in `outputs/order0_table_integral/order0_observables.npz` |

The four signed observables are convergence, xi-plus, xi-minus, and
convergence–tangential shear in the manuscript's pair-aligned screen basis.
Both calculations use the same finite full-sky angular sums, with the appropriate
spin factors. Neither calculation contains nonlinear optical F or K insertions.
The distinction between optical perturbative order and Limber is independent
of this implementation comparison.

## Current calculation

[source_folded_reference.py](source_folded_reference.py) numerically integrates
the continuous source response on 100001 nodes, checks it against adaptive
quadrature, and registers the resulting weighted Ricci and Weyl operators in
an isolated calculation. The field derivatives retain their original meaning.
The physical lightcone is unchanged. The old table's artificial lower source
floor at lambda = 406 Mpc is not imposed.

The accepted full-Sachs product is
[source_folded_full_boundary_n4096_k2048.npz](source_folded_full_boundary_n4096_k2048.npz),
with its [run metadata](source_folded_full_boundary_n4096_k2048.json):

- Below the split, the existing canoes direct-k helper is loaded in a separate
  namespace, using 4096 radial nodes and 2048 k nodes over 0.0001–10 Mpc^-1.
  The finite lower-radial-end integration-by-parts terms are included explicitly.
- All ell = 2–30 and 49 logarithmically spaced integer nodes from 31–800 are
  evaluated directly. Fixed-sign spectra are interpolated with a cubic spline
  in log ell and log absolute spectrum. Every ell = 801–5000 is evaluated using
  the canoes ID-basis Limber approximation.
- Bessel values use a bounded 512 MiB cache keyed by the exact input content.
  [The cache check](source_folded_cache_validation.json) found bitwise-identical
  cached and uncached spectra in the three-multipole probe. The installed
  canoes and SciPy modules are not changed.
- Full Ricci, Ricci–Weyl, and Weyl spectra are saved. Transverse-only Ricci
  spectra are retained separately using their exact angular-factor relation
  to the Weyl channel.

The convergence product
[source_folded_withheld_n2048_k1024.npz](source_folded_withheld_n2048_k1024.npz)
evaluates the
withheld midpoints between the sparse nodes with 2048 radial and 1024 k nodes.
It checks interpolation together with radial and k refinement. Its results
must be read with the other convergence probes in the final validation record.

[ccl_reference.py](ccl_reference.py) supplies the conventional transverse
lensing projection using CCL 3.3.1. The accepted benchmark
[ccl_upperpad_n8192.npz](ccl_upperpad_n8192.npz) uses
8192 FKEM radial samples, with zero padding below the physical 50 Mpc cutoff
and beyond the source through chi = 13500 Mpc. FKEM applies through ell = 800,
and CCL's Limber integration applies above it. Separate angular prefactors for
convergence and shear are checked against direct tracer calculations.
The earlier [ccl_padded_n4096.npz](ccl_padded_n4096.npz) is preserved.

## Reproduce the accepted settings

Run sequentially from this directory with the PyCCL environment. These commands
write new `_rerun` products to preserve the existing products. Choose another
unused suffix for subsequent runs. The source metadata and shared angle file
under `outputs/order0_table_integral/` are prerequisites.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 CANOES_SUPPRESS_METAL_WARNING=1
R1_PYTHON=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python

"$R1_PYTHON" source_folded_reference.py \
  --ells sparse --sparse-nodes 49 --n-chi 4096 --low-ell-n-k 2048 \
  --cached-direct --observer-boundary --ell-block 4 \
  --channels cl_pp,cl_pw,cl_ww --split 800 --limber-shift 0.5 \
  --out source_folded_full_boundary_n4096_k2048_rerun.npz

"$R1_PYTHON" source_folded_reference.py \
  --ells withheld --sparse-nodes 49 --n-chi 2048 --low-ell-n-k 1024 \
  --cached-direct --observer-boundary --ell-block 8 \
  --channels cl_pp,cl_pw,cl_ww --split 800 --limber-shift 0.5 \
  --out source_folded_withheld_n2048_k1024_rerun.npz

"$R1_PYTHON" ccl_reference.py \
  --fkem-nchi 8192 --radial-zero-pad-max 13500 \
  --z-s 5 --chi-min 50 --ell-max 5000 --l-limber 800 \
  --out ccl_upperpad_n8192_rerun.npz
```

Each run writes an NPZ plus JSON settings, input/script hashes, warnings, and
timings. The source-folded run also records the direct-k helper hash and cache
statistics. [plot_comparison.py](plot_comparison.py) accepts `--sft`,
`--reference`, and a fresh `--output-dir` for the comparison and residual plots.

## Validation and interpretation

The numerical convergence target is a maximum absolute change below 0.1% of
each observable's reference peak on the shared 40-angle grid. The final record
identifies each comparison and separately reports relative differences where
the reference exceeds 1% of its peak. Sign crossings remain visible in the
signed curves.

- Inserting the withheld nodes computed with half the radial and k resolutions
  changes the four angular statistics by at most 0.01616% of their peaks.
- Doubling CCL's FKEM radial samples from 4096 to 8192 changes the probe spectra
  by at most 0.002075%.
- The dense Ricci/Weyl covariance satisfies positivity. The transverse angular
  prefactor identities hold to a relative precision of 6.7e-16 after interpolation.
- The final source-folded run took 300 seconds with peak resident memory of
  1.051 GiB. Its saved `.source.py` snapshot matches the run's script hash.
  The saved snapshot preserves that exact run. The current script additionally
  supports source-map selection with `--z-source` and a scalar-only transform
  with `--channels cl_pp`. Source-coordinate round trips and bitwise parity of
  the scalar transform against accepted z = 5 are recorded in
  [source_folded_multiz_source_checks.json](source_folded_multiz_source_checks.json).

Keep these limits explicit when interpreting the cross-method comparison:

- The full Ricci operator contains time and radial derivative terms absent
  from conventional transverse CCL convergence. Compare the saved transverse
  channel as an additional operator-matched diagnostic. A residual between
  the full and transverse calculations is not itself a numerical error bar.
- The source-folded high-ell calculation retains the ID-basis approximation.
  It does not evaluate all full-Ricci derivative terms above ell = 800.
- FKEM integrates over a finite returned k grid. Upper radial zero padding
  reduces the omitted low-k tail without eliminating it. The lowest-ell
  discrepancy has a finite-k contribution and must not be attributed solely
  to the Ricci operator or treated as a canoes defect. CCL metadata records
  the returned k ranges and this caveat.
- The [source-window diagnostic](outputs/high_ell_window_diagnostic.json)
  traced the earlier roughly 6–8% observable bias to source-table interpolation.
  The continuous source integration removes that interpolation step.
- Earlier full-Ricci t-form probes failed covariance positivity, as recorded in
  [source_folded_validation.json](source_folded_validation.json). They are
  diagnostics, not candidate physical predictions. Unpadded CCL products in
  [ccl_archive_unpadded_20260928](ccl_archive_unpadded_20260928/) and
  [rejected comparison plots](rejected_comparisons_archive_2026-09-28/README.md)
  are retained for provenance and are not replacement figures.

Inspect the final validation record and rendered plots before deployment.
Manuscript figures are deployed only through `reproduce/deploy.py` after
review. This README does not certify unreported convergence or deployment.
