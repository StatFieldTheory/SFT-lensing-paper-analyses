# Monte Carlo wrapper handoff

## Status and ownership

Prepared and smoke-tested, not independently reviewed or accepted for final
figures. All changes are confined to `mc_update/`. No sampler, installed
package, manuscript, existing marker or production cache was modified.

All current subprocesses have finished. No heavy MC job has been launched.
Root coordinates resources and the combined code/Python review before further
production work. Run each FK block as a fresh process with one BLAS thread.

## Files to review

`runtime.py`, `run_fk.py`, `pool_fk.py`, `fk_quadrature.py`, `run_ff.py`, and
`build_cache.py`. All six pass `py_compile`. The actual `run_fk` and `run_ff`
paths passed numerical smoke tests. Pooling, fold provenance checks and cache
assembly await end-to-end checks with the completed new main fold and full
marker products.

The expected comparison file is confirmed by the actual prepared manifest:
`../true_redshift_fk_gl24_corrected_main/main/xi_main.npz`. Its manifest is
`../true_redshift_fk_gl24_corrected_main/prepared.json`, not inside `main/`.
The loader checks the recorded table and response hashes in addition to the
actual affine endpoint in every selected row.

## Sampling and reference boundary

The covariance PSD repair does not deliberately change the target K. The
wrapper constructs exactly the L factors the existing sampler uses, propagates
their covariance, and calibrates the existing `solve_Q` against that covariance.
The target remains `driver_stats.zeta6` from the corrected table. Its residual
interpolation asymmetry is automatically symmetrized by `solve_Q`, which gives
a -0.121010% target-to-injected FK difference at the checked 1-arcmin angle.
The tensor peak-relative injection error is 0.209376% there. Both are saved,
not absorbed into a fitted normalization.

At that angle the finite-correlation extrapolated exact expectation differs
from the injected local reference by +0.0924085%. This excludes finite-batch
centering bias and lambda-step bias. Other angles must be checked separately.

Suggested numerical gates for retaining a percent-level FK agreement claim:
the weighted target/injection difference and finite-correlation extrapolation
residual should each remain below 0.3%, and source-grid refinement should be
below 0.3%. These are proposed accuracy targets, not passes already established
over the requested angle range. If a gate fails, report the measured bias and
adjust the validation claim or sampling resolution, rather than renormalizing.

For FF, use the sampled local covariance reference, including the original
160-node Sigma2 differentiation, endpoint behavior and L_white PSD clipping.
No claim of a matched comparison to the full nonlocal cosmological FF is valid.
The continuum contraction and finite Euler recursion remain distinct. At
N_lambda=4000 the disconnected part alone differs by 0.553% at 1 arcmin.
The full FF contraction changes by 0.0264% between 481 and 961 reference nodes.
An N_lambda=4000/8000 comparison and selected-point higher-F-order check remain
necessary if a precise agreement percentage will be claimed.

## Runtime evidence and production commands

See `README.md` for exact commands. The full four-block FK design historically
took 55 minutes serial. The current small N_lambda=1000 probes take about
7 seconds per gamma/sigma row including expectation and one 4000-realization
seed. This short probe does not isolate enough of the fixed versus sampling
cost to promise a new full runtime. Allow roughly 1 to 2 hours serial initially
and use the first full row as the updated estimate.

Historical FF N_lambda=4000, eight times 96000 realizations took 165 seconds per
angle. Four angles would be about 11 minutes historically. The 4000-node
reference alone, including 241/481/961 contraction grids, currently takes
1.4 seconds at 1 arcmin. A dense 40-angle reference curve is inexpensive relative
to MC and can be produced in a separate `--reference-only` directory.

Use `--reference-n 481 961` for the final FF marker/reference run. Then run one
selected angle at N_lambda=8000 to quantify the remaining integration change.
The code currently reuses the original nonlinear F-on/F-off estimator, so a
perturbative-order test still needs either an existing order-controlled sampler
or an explicitly reviewed coupling-strength study.

## Plot cache schema

`pool_fk.py` produces the requested `gamma`, `fk`, `err`, `analytic` keys,
plus `blocks`, `extrapolated_expectation`, `ref_inj`, `ref_true`, `n_blocks`,
`seeds_per_block`. Its JSON records fold/input hashes and the separate error
ratios.

`build_cache.py` produces:

| Keys | Meaning |
|---|---|
| `g`, `o0`, `fk3` | New main fold angular grid, its Order-0 and FK convergence |
| `g_ff`, `ff_mom` | Continuum local-covariance FF reference on its own grid |
| `g_mc`, `ffc_mc`, `ffc_se` | FF local-driving MC means and sampling errors |
| `g_fk`, `fk_mc`, `fk_se`, `fk_analytic` | FK markers and updated denominators |
| `source_redshift`, `lambda_source` | Explicit actual z_s=5 source |
| `ff_reference_label` | `FF (local covariance)` |

The current production plot assumes the FF line uses `g`. Root must change
that line to use `g_ff` when present and identify the local-covariance reference
in its label and caption. This keeps the existing single-panel layout if desired.
A dense reference curve should span the plotted marker range. No missing-angle
values are filled or extrapolated by the cache builder.

The cache's Order-0 guide line currently comes from Order-0 rows included in the
new main FK fold. If root elects a different accepted Order-0 guide product,
its own angular and source provenance must be explicitly adapted and recorded.

## Serial continuation dispatched

The root authorized the following continuation after the reviewed four-angle
`ff_n4000` production run. It is running outside the agent's lifetime, with one
BLAS thread in each process.

1. Session **6514** waits for `ff_n4000/complete.json`, then runs gamma=1,
   N_lambda=8000, eight seeds of 96000 realizations, reference grids 481 and 961,
   seed base 7010001. Output: `ff_step_n8000/`, log `ff_step_n8000.log`.
2. Session **51973** waits for that completed N8000 run, then evaluates the
   N_lambda=4000 local-covariance continuum reference on all 40 normalized
   angular positions in the new main fold configuration. Output:
   `ff_dense_reference_n4000/`, log `ff_dense_reference_n4000.log`. Its request
   and source-configuration hash are in `dense_reference_request.json`. It
   exits without launching if its predecessor does not complete within one hour.

No half-coupling MC or FK MC has been launched by this worker. Root owns the FK
schedule. The new `run_ff_coupling.py` awaits independent general/Python review.
Its syntax and proxy identity check passed: the proxy changes only the original
vertex coefficient and retains the identical original precompute function.
After approval, run it serially with other MC jobs:

```bash
python run_ff_coupling.py --baseline ff_n4000 --gamma 1 --strength 0.5 --out ff_coupling_half
```

This uses all baseline sampling settings and identical random seed streams,
normalizes by the squared coupling and records the paired seed-difference
error. It measures finite-coupling dependence, not an exact F-squared limit.
The root explicitly authorized this selected diagnostic after review, without
additional coupling amplitudes or increased precision if the result is small.

At handoff, the first three full N4000 MC points differ from their continuum
local references by -0.6240%, -0.5802% and -0.4792%, at 1, 2.6 and 6.7 arcmin.
Their sampling standard errors are 0.266%, 0.266% and 0.273% of the references.
These remain preliminary until the selected integration-step and higher-F
diagnostics finish. The fourth angle was still running at dispatch.

## Interrupted continuation and recovery

The earlier assumption that worker sessions would survive agent handoff did not
hold for this run. At the next inspection both sessions 6514 and 51973 returned
`Unknown process id`, including from the originating agent. Host process
inspection confirmed no `run_ff.py` job remained. The exact termination cause
was not diagnosed. The N8000 run preserved six complete seeds, while the dense
reference never started. Root identified the remaining generic Python watcher
PIDs as its own queues, so none was stopped or duplicated.

The new `resume_ff.py` recovers only the two missing seeds after verifying all
original input, source and runner hashes. `--check-only` passed against the
actual checkpoint. The missing seed indices are 6 and 7, with RNG seeds
7015947 and 7016938. It preserves `seeds_g00.partial6.npz`, records original
manifest/checkpoint hashes in `recovery.json`, writes checkpoints atomically,
recomputes the inexpensive references and writes the original completion
sentinel atomically last.

No recovery computation has been launched by the worker. Root must run the
helper after independent code/Python review, then launch the dense reference
under a root-owned execution session. Review delegation was attempted but the
four-agent concurrency limit blocked spawning the reviewer.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python resume_ff.py --run ff_step_n8000
```

Use the established PyCCL interpreter and an unused recovery log. The original
manifest remains unchanged. The new reference JSON labels its elapsed time as
the recovery continuation only. The earlier six seed durations are unavailable.
