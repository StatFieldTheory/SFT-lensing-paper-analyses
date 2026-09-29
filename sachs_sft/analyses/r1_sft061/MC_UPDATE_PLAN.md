# Monte Carlo validation update at the corrected source plane

## Required scope

The author approved the corrected kappa3 inputs and consistent figure and text
updates. The source plane is **z_s = 5**, with
`lambda_source = 2317.9695958203943` physical Mpc. There is no exception retaining
the historical endpoint for this figure. Existing MC realizations must not be
relabeled or rescaled into this source plane.

This is an execution handoff, not evidence that the new computations passed.
No production files or numerical products were changed while preparing it.

All paths below are relative to `SFT-lensing-paper-analyses/` unless absolute.

## Shared inputs and process isolation

- Python: `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`.
- Background and response: `sachs_sft/analyses/r1_sft061/aligned_response.py`.
  Use its `D_at`, `chi_of_lambda`, `z_of_lambda`, and `response_R` directly.
- Corrected production K input:
  `sachs_sft/callables/kappa3_vertex/rebuild/products/pieces_nphi512/table_permclosed_np512.npz`.
- Callable:
  `sachs_sft/callables/kappa3_vertex/equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable.py`.
- FK comparison fold: the **new corrected-table, aligned-response, true-z5**
  main fold when available. Do not use `fk_same_inputs_gl24/xi.npz`, which is an
  old-endpoint version-control experiment.
- Covariance source remains the recorded corr_op table unless the root task
  explicitly changes it. Record its hash and the Sigma2 construction separately.

Use fresh output directories and one fresh subprocess per MC block. Pass table,
fold, endpoint, response module, and output paths as explicit wrapper arguments,
and write them with hashes to a sidecar manifest. Do not edit historical modules
or production defaults.

### Does this require changing the callable's module global?

The existing callable has a module-global `TABLE_PATH` and `_CACHE`, and its
factory has no table-path argument. Merely starting a subprocess does not make
its default input configurable. Two isolated routes are available:

1. **Smallest change:** after all legacy imports, set `pa.TABLE_PATH` and clear
   `pa._CACHE` in the fresh child process, then bind `driver_stats._k3 = pa`.
   This mutates only child-process memory, not a file or another running job.
   Perform no further `_bootstrap.wire()` call afterward because it resets the
   routing. Assert and record the active table path before creating any nodes.
2. **No override of `pa.TABLE_PATH`:** stage a byte-identical callable module in
   the run directory with the corrected input available there as
   `table_permclosed.npz`, import it under a unique module name, and bind
   `driver_stats._k3` to that isolated module. Its existing default resolves to
   the staged input. Record module and table hashes. This still needs the
   driver binding because `driver_stats` otherwise imports its production
   module by name. A copied input or a read-only symlink is sufficient.

The first route is preferable unless avoiding even an isolated module-global
override is a specific requirement. Avoid modifying installed canoes or the
production callable.

### Aligned background adapter

Legacy `mc_sachs_2pt/background.py::Background` samples the old `D_callable` onto
another spline. Replacing only its endpoint is insufficient. Provide a small
process-local adapter accepting its constructor keywords and exposing `D` via
the aligned `D_at`, `response_R` via the aligned function, and the aligned `z`,
`chi`, and `a` accessors if needed. Bind the background module's `Background`
before constructing grids or Sigma2 builders. The current MC integration uses
ratios of `D` directly and does not require `theta_sa`.

Pass the exact endpoint explicitly to `fk_expect_exact.build_grid` and
`sachs_mc_core.MCConfig`. Changing `background.LAM_SOURCE_BASELINE` after import
does not change the already evaluated dataclass default in `MCConfig`.

Check MC response ratios against `aligned_response.response_R` at the actual
grid nodes. Keep the physical affine units and the existing lower support
boundary of 406 Mpc explicit.

## FK markers

### Reusable implementation

- `sachs_sft/analyses/mc_fk_complete/run_sweep.py::main` already runs the required
  gamma, correlation-length, and seed sweeps, records per-seed values, and
  computes the estimator expectation and injection checks.
- `fk_expect_exact.py::build_grid` accepts `lam_source` directly.
- `fk_expect_exact.py::node_stats`, `calib_matrix`, `zeta_injected`,
  `fk_reference`, and `expectation` provide the existing validation diagnostics.
- `fk_complete_core.py::simulate_fk_complete` accepts explicit `grid` and
  `nodes`, so the simulation itself needs no change.
- `run_gate2.py::analytic_fk` reads `_bootstrap.FOLD`. A wrapper can override
  that path after imports, or supply its own strict fold loader.
- `pool_markers.py` documents the existing pooling rule, but hardcodes old
  filenames and overwrites `_markers_pooled.npz`. Reuse its calculation in the
  new wrapper with explicit inputs and a fresh output path.

`run_sweep.py` does **not** currently accept table, fold, response, or endpoint
CLI arguments. The following is the reusable argument tail for a new isolated
wrapper, not a command that updates the old script correctly by itself:

```text
--gammas 1.0 2.6 6.7 17.3 --sigmas 8.0 4.0
--n-lambda 1000 --n-real 24000 --batch 2000 --seeds 24
--calibrate smeared --seed-base <block-base> --out <fresh-block.npz>
```

Use the four distinct block bases `6010001`, `6020002`, `6030003`, and `6040004`
to reproduce the existing sampling design. Reusing seeds at changed inputs is
acceptable and improves controlled comparisons. Keep the sigma=8 and sigma=4
seed lists paired within each block. Pool the extrapolated block means and use
the scatter of independent block means for error bars, as in `pool_markers.py`.

Recommended staged execution:

1. Prepare one gamma=1, sigma=8/4 grid pair. Check routing, node response,
   injected versus target cumulant, and sampled covariance versus the covariance
   assumed by `expectation`. `simulate_fk_complete` uses `_cholesky_psd`, which
   can clip negative eigenvalues. The exact expectation is not automatically an
   exact reference for the sampled field if it still uses unclipped matrices.
2. Run one short seed at each sigma to measure current runtime and check finite
   outputs and retained-realization counts.
3. Run the four full blocks, checkpointing each gamma/sigma result rather than
   only saving at the end of an eight-row block.
4. Verify finite-correlation and lambda-step errors with the existing exact
   machinery at selected points before keeping the historical linear
   extrapolation claim. The changed K input does not guarantee the old residual
   estimates still apply.
5. Recompute marker/fold ratios using the new main FK product. Do not preserve
   historical `analytic` denominators in the new marker file.

### Runtime evidence

Historical `_blk6010001.log`, `_blk6020002.log`, `_blk6030003.log`, and
`_blk6040004.log` contain eight row durations each. Their sums are 829, 829,
819, and 829 seconds, respectively: **3306 seconds, about 55 minutes serial**
for the full four-block production design, plus preparation. This is historical
evidence, not a current-machine guarantee. Each row took roughly 97 to 109
seconds for 24 seeds. A short current probe should replace this estimate.

## FF: a consistent validation requires a matching covariance reference

The old FF MC does not use the full nonlocal cosmological C of the main FF
fold. It uses a local-in-lambda Gaussian drive obtained by differentiating the
equal-time diagonal of corr_op. This distinction must survive the update.
Changing the endpoint or increasing MC statistics does not remove it.

### Existing reusable reference

`sachs_sft/analyses/mc_fk_complete/inputmatch/ff_functional.py::ff_terms`
implements the existing FF contractions for an arbitrary state covariance and
returns the connected, disconnected, and full moment terms. Its `kind="full"`
uses corr_op. Its `kind="mc"` uses a continuum local-covariance surrogate with
the lower-bound contribution subtracted. No new symbolic derivation is needed
to reuse those contractions.

Historical `inputmatch/ff_functional.npz` already exposes the difference. At
1 arcmin it stores FF_full=1.94914e-6, FF_mc=1.71872e-6, and the historical
workflow reference=1.95357e-6. Tightening the MC errors therefore cannot justify
the old claim of agreement with the full cosmological reference.

### What must be aligned beyond stock `kind="mc"`

The stock surrogate is **not** exactly the finite process the sampler realizes.
The sampler's `Sigma2Builder` differentiates a 160-node cubic spline, clips at
the support endpoints, and `_cholesky_psd` can alter the covariance. The FF
vertex also uses the preceding time node in the Euler step.

For a strict matched reference:

1. Obtain `lam`, `dlam`, `resp`, and `L_white` from the same
   `sachs_mc_core._precompute(cfg, cos_gamma)` used by the MC.
2. Reconstruct the Gaussian state covariance by propagating the **actual**
   per-node noise covariance `L_white @ L_white.T` with those response ratios.
   This includes PSD clipping automatically. Keep the full six-component
   covariance, then extract the two within-ray blocks and the cross-ray block.
3. Supply those covariance blocks to the existing FF contraction routine.
   A continuum evaluation on a coarser grid requires a checked MC step-size
   limit. An exact finite-grid prediction additionally needs the existing
   pre-vertex-node convention and discrete response/window weights. Do not call
   a plain substitution into the continuum `ff_terms` an exact Euler reference.
4. `psd_fix/r6_exact_ff_mean.py` already implements the discrete mean and
   disconnected-channel check and is a useful independent consistency gate.
5. The MC F-on/F-off estimator contains higher F orders, while `ff_terms` is
   Order F-squared. Check that distinction at selected points or state it in
   the validation description. Do not silently equate the two objects.

The minimal practical route is the continuum matched local-covariance
reference with a small step-size study of the MC, using the existing
contractions and the actual sampled covariance. This validates this local
driving model and must not be labeled a numerical check of the full nonlocal
cosmological FF result.

### Reusable MC engine and cost

- `inputmatch/ff_anti.py::ff_moment_anti(mc, cfg, gamma)` implements the same
  F-on/F-off integration with antithetic Gaussian draws. It removes odd-noise
  fluctuations and returns the mean, pair-based standard error, pair count,
  and blowup count.
- Construct `MCConfig(lam_source=2317.9695958203943, n_lambda=..., n_real=...,
  batch_size=3000, apply_anchor=False, seed=...)` explicitly. Use the aligned
  background adapter described above.
- `inputmatch/run_nl_scan.py` is an existing step-size driver, but it hardcodes
  the old source through config defaults and reads the old full-FF reference.
  Reuse its loop, not its unmodified command.
- Historical `inputmatch/out/nl_g1_off.log` records eight seeds times 96,000
  realizations at 1 arcmin: N_lambda=1000/2000/4000/8000/16000 took
  40/81/165/351/693 seconds. The high-statistics 4000-node calculation therefore
  costs about 2.75 minutes per angle historically. All 15 existing marker
  angles at that precision would be about 41 minutes, before reference work.
  Smaller pilot statistics are appropriate for measuring current runtime.

## Figure, text, and number bookkeeping

- `mc_sachs_2pt/fig_xi_channels.py --from-cache` reads all curves and markers
  from `outputs/appendix_mc_curve.npz`. It ignores alternative marker arguments
  when the cache is used. Use an explicit new cache and generator route.
- The existing plot has one O0 line, one FF line, one FK line, and two marker
  series. If retaining the full cosmological FF line, add a clearly labeled
  local-covariance FF reference for the MC comparison. Alternatively split
  the validation channels into two compact panels. Do not place local-model
  markers over a full-C line and describe them as a matched-input test.
- A minimal fallback, if the matched FF reference cannot be validated promptly,
  is to show the newly computed FF markers as a labeled local-covariance
  comparison and remove the claim that their apparent agreement validates the
  full nonlocal FF calculation. Quantify the actual difference instead.
- The caption and Appendix validation paragraph must state true z_s=5 and
  distinguish the FK matched-cumulant test from the local-Gaussian FF test.
  Revise the related summaries in `sections/insights.tex` and
  `sections/conclusion.tex` only after the new numerical checks finish.
- Recompute the independent FK quadrature at the new endpoint, aligned D, and
  corrected table. `fk_kernel_crosscheck.py` currently hardcodes the old
  endpoint and a single ordered contraction. Prefer the existing
  `fk_expect_exact.fk_reference` with the verified symmetric corrected input,
  or review/update the old quadrature's input and ordering explicitly.
- Point `reproduce/check_numbers.py` at the new marker and reference products.
  All percent-agreement and error-bar claims must be regenerated.
- Archive historical caches and plots. Deploy the accepted figure only through
  `reproduce/deploy.py`.
