# Corrected-input Monte Carlo validation

These adapters reuse the existing samplers and contraction routines in fresh
Python processes. No historical module or MC product is modified. They use the
aligned response at actual source redshift 5, affine endpoint
2317.9695958203943 physical Mpc, and the corrected n_phi=512 K table at
ell_max=15360. The lower support remains 406 Mpc.

The numerical convention is Born/FLRW with E0=1. Each run records input hashes,
the old module hashes, sampling settings and its comparison scope. Current
source-of-truth plan: `../MC_UPDATE_PLAN.md`.

## Implementations

- `runtime.py`: imports all legacy modules before a process-local K/background
  override. The new comparison fold must have a matching prepared manifest,
  K hash, response hash and endpoint.
- `run_fk.py`: one checkpointed seed block, preserving the paired sigma=8/4
  design. Its covariance expectation uses the actual sampled L L-transpose
  covariance after PSD clipping. Existing `solve_Q` calibrates the injected
  cumulant against that same covariance. Population-centered exact expectations
  and finite-batch sample centering remain distinct.
- `pool_fk.py`: the historical linear sigma extrapolation and error from
  independent block means, with checks against duplicate seed streams. It
  recomputes denominators from the new fold and records injection, correlation
  length and grid effects separately.
- `fk_quadrature.py`: the existing target-cumulant contraction at successively
  finer source grids, independent of the framework fold and the MC sampling.
- `run_ff.py`: reconstructs the full six-component state covariance from the
  sampler's actual `L_white`, including PSD clipping. It reuses the existing
  `ff_functional.ff_terms` contractions with this covariance. The result is a
  continuum Order-F-squared reference. It is neither the nonlocal cosmological
  FF prediction nor an exact prediction for the finite Euler step. The exact
  discrete Order-F mean and disconnected term are recorded separately, using
  the existing preceding-node convention.
- `build_cache.py`: assembles a new figure-6 cache from completed FF reference,
  FF MC and pooled FK products. It validates source/provenance consistency.
  The FF curve has its own `g_ff` angular grid, and its label must identify the
  local covariance. Existing `g`, `o0`, `fk3`, `g_mc`, `ffc_mc`, `ffc_se`,
  `g_fk`, `fk_mc` and `fk_se` keys remain available.

The F-on/F-off FF sampler includes higher F orders. Step refinement and a
selected-point higher-order assessment are required before a percent-level
validation claim. FK finite-batch centering can add bias of order 1/batch.

## Probe evidence

The following are smoke/profiling products, not publication-ready markers.

- `smoke_fk_n96/`: one angle, 400 realizations, two correlation lengths.
- `smoke_ff_n96/`: one angle, 400 realizations. Its deliberately coarse Euler
  step visibly differs from the continuum contraction.
- `probe_fk_n1000/`: one angle, 4000 realizations, one seed per correlation
  length. About 7 seconds per row including the exact expectation. At 1 arcmin
  the sampled-covariance PSD change is 9.0e-6 of its peak. The local injected FK
  differs by 0.121% from the tabulated target, reflecting residual input
  interpolation asymmetry. The linear sigma extrapolation of the exact
  expectation differs from the injected local reference by 0.0924%.
- `probe_ff_n1000/`: one angle, 4000 realizations. About 1.1 seconds total.
  FF = (1.997 +/- 0.073)e-6, versus the sampled-covariance continuum reference
  1.9864e-6. Refining its reference grid from 241 to 481 nodes changes it by
  0.105%. The probe's broad statistical error does not establish final accuracy.
- `reference_ff_n4000/`: at 1 arcmin, continuum FF is 1.9630003e-6 on a
  961-node reference grid. Refining 481 to 961 nodes changes it by 0.0264%.
  The continuum disconnected term is still 0.553% above its exact Euler
  counterpart at 4000 steps. This illustrates why the two predictions are not
  silently equated.

The wrappers changed after some probes to strengthen fold provenance checks.
Each probe's sidecar records the exact implementation used at its run time.

## Commands after code review and resource authorization

Use `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python` and set
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`. Paths below assume this directory is
the working directory. Run each block in a separate Python process. Do not
reuse output directories.

```bash
python run_fk.py --seed-base 6010001 --out block6010001
python run_fk.py --seed-base 6020002 --out block6020002
python run_fk.py --seed-base 6030003 --out block6030003
python run_fk.py --seed-base 6040004 --out block6040004
python pool_fk.py --blocks block6010001 block6020002 block6030003 block6040004 --out fk_markers.npz
python fk_quadrature.py --out fk_quadrature.json
python run_ff.py --gammas 1 2.6 6.7 17.3 --n-lambda 4000 --out ff_n4000
python build_cache.py --fk-markers fk_markers.npz --ff-reference-dir ff_n4000 --ff-mc-dir ff_n4000 --out appendix_mc_curve.npz
```

The FK comparison fold must exist before these production commands. It is
`../true_redshift_fk_gl24_corrected_main/main/xi_main.npz`. Probe runs can use
`--allow-missing-fold`, but those products establish no comparison to the paper
fold until a later explicit comparison is performed.

## Figure and caption boundary

A minimal honest figure can use two panels: FK versus its matched corrected
cumulant prediction and FF versus its sampled local-covariance reference.
Four marker angles shared across panels are sufficient for the focused check.
An optional full cosmological FF curve must have a separate label and cannot
serve as the matched reference for these local-driving markers. Error bars
describe sampling uncertainty, with the measured numerical biases reported
separately in the validation text.

The root task owns plotting, manuscript changes and deployment. None has been
performed here. Heavy production runs await root resource coordination and
independent code/Python review.
