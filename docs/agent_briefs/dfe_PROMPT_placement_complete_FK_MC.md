# Task: a placement-complete, variance-controlled Monte-Carlo test of the FK channel

## What this is for

The STF_lensing paper claims that a path-integral (SFT-Wick L2) solution
reproduces the statistics of the stochastic Sachs system. That claim is
currently tested against a direct numerical solution of the stochastic
equations **only for the FF channel**. The FK channel, which is the paper's
novel object (the hierarchy-mixing diagram, a three-point driving-field
cumulant entering a two-point function through one nonlinear Sachs
correction), has **never** been compared against a direct stochastic solution.
Your job is to build that comparison.

The paper is under review and its current wording is honest about this: it
says the FF check is "independent of the diagram expansion" and the FK check
"independent of the workflow's implementation of it". If you succeed, that
sentence can be strengthened. If you fail, it stays as it is. Do not weaken
or overstate it either way.

## Why the three existing estimators do not do this

All three live in
`SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt/sachs_mc_core.py`.
Read all three before writing anything.

| routine | placement-complete | variance | verdict |
|---|---|---|---|
| `simulate_crn` (demo2 deformation, line ~258) | yes | blows up | per-realisation product amplified by Q^2, Q ~ 1e10, heavy-tailed |
| `simulate_fk_pathA` (line ~354) | yes | pathological | Q sits inside the per-realisation evolution |
| `simulate_fk_vr` (line ~444) | **no** | controlled | Q pulled outside the average, at the cost of dropping placements |

`simulate_fk_vr` is what the published figure used. An exact
second-order expectation analysis (verified against 21 probe points to about
1.5 seed standard errors) showed that it does not estimate the FK diagram:

* It injects the third cumulant **only on the observable-side leg** (its
  docstring calls this the "no s_d ODE" shortcut). Of the three placements
  making up `Zeta = Q:VV + perms` it captures one. In the eigenbasis of the
  node covariance the captured share is
  `d_q d_r / (d_p d_q + d_q d_r + d_r d_p)`, which kernel-weighted is
  **0.245 at 0.5', zero near 10', negative beyond**.
* A second, finite-`sigma_lambda` term enters because the AR(1) smearing
  decoheres the node-local `Q` from the smeared covariance, releasing the
  cancellation of the tabulated tensors' lambda-grid graininess.
* The historical "agreement" at `sigma_lambda = 8` was those two errors
  summing through unity at a value that had itself been tuned so the ratio
  came out one. The FK amplitude is close to linear in `sigma_lambda`
  (8.70, 19.8, 37.9 x 1e-6 at 0.5' for sigma = 4, 8, 16), intercept
  consistent with zero.

So the FK markers were withdrawn from the paper's validation figure. Do not
reinstate them from `simulate_fk_vr` under any circumstances.

## The design to implement

`simulate_fk_pathA` is placement-complete: its skew increment
`df = 0.5 Q (z x z - V)` is fed through the **linearised** F-vertex into a
second field `s_d`, so it reaches the F-legs as well as the observable leg.
Its problem is only that `Q` multiplies a per-realisation quantity.

But `s_d` is **linear** in its drive. So `Q` can be pulled outside the average
exactly as `simulate_fk_vr` does, without dropping any placement:

1. Replace the scalar drive `df^a = 0.5 Q^a_{ij} (z_i z_j - V_ij)` by the
   **tensor-valued** drive `(z_i z_j - V_ij)` with the `Q` contraction removed.
2. Evolve a tensor-valued increment field `s_d^{a,ij}` under the same linear
   ODE, `ds_d = resp s_d + DF(s_g).s_d dlam + (zz - V) dlam`, and integrate it
   to `kappa_d^{a,ij}`.
3. Accumulate the covariance `C^{ab,ij} = <kappa_g^a kappa_d^{b,ij}>` over
   realisations, subtracting the disconnected piece with the **matched sample**
   means (see the comment at `simulate_fk_vr` lines ~494-499 explaining why the
   sample `<zz>` must be used, not the theoretical `V[k]`).
4. Contract with `Q` **outside** the average, and symmetrise over the two rays
   exactly as Path A does (`FK = <kg x kd> + <kd x kg>`).

`(z_i z_j - V_ij)` is symmetric in `(i,j)`, so there are 21 independent
components for a 6-vector field: the linear ODE is solved 21 times per
realisation instead of once. Expect roughly a 20x cost over `simulate_fk_vr`
(which is ~17 s per (gamma, seed) at `n_real=24000`, `n_lambda=1000`), so
budget a few minutes per (gamma, seed). Vectorise over the 21 components with
numpy rather than looping if you can; the batch axis is already there.

Put the new routine in
`driver_field_emulators/code/mc_fk_complete/` as a **new file** that imports
from `sachs_mc_core`. Do NOT edit `sachs_mc_core.py`, and do not touch any
file under `SFT-lensing-paper-analyses/` or `sections/` or `figures/`.

## Stopping rule, in this order

Do not skip ahead; each gate is cheap relative to the next.

1. **Placement share.** With `sigma_lambda` fixed at 8 and a single gamma
   (1 arcmin), check that the new estimator's expectation is the FULL vertex,
   not 0.245 of it. The cheapest check is deterministic: reproduce the
   share calculation for the new placement structure and confirm it is 1.
   If it is not 1, stop and report why.
2. **One point.** Run at gamma = 1' and compare with the analytic FK,
   `driver_field_emulators/products/table_permclosed_cut15360_permfix_xi.npz`
   (kappa-kappa, order 2; the loader pattern is in
   `driver_field_emulators/code/rebuild/fk_kernel_crosscheck.py`). Report the
   ratio and its seed-to-seed scatter over at least 8 seeds. If the scatter
   swamps the signal, stop and report the variance, do not throw realisations
   at it.
3. **sigma_lambda dependence.** The analytic FK is the white-noise limit, so
   the honest test is an extrapolation, not a single value. Run
   `sigma_lambda` in {4, 8, 16} at gamma = 1' and report whether the ratio to
   the analytic value approaches 1 as sigma decreases. Watch for the known
   trap recorded in `fig_xi_channels.py`: at fixed sigma a FINER lambda grid
   makes things worse (heavy-tailed skewness), so `n_lambda` is held fixed and
   `sigma_lambda` alone is swept.
4. **Only then**, a gamma sweep over 0.5' to about 20'. Beyond roughly a
   degree the analytic FK itself is not converged (the multipole sum's low-
   and high-ell branches cancel), so do not compare there and do not treat
   disagreement there as a finding.

## Ground truth for comparison

Analytic FK (kappa-kappa, converged cutoff `ell_max = 15360`, permutation-closed
vertex), from the npz above:

| gamma | FK |
|---|---|
| 0.5' | +1.948e-05 |
| 1.0' | +1.598e-05 |
| 2.0' | +1.223e-05 |
| 5.0' | +7.038e-06 |
| 17.3' | +1.966e-06 |

Cross-checked by an independent quadrature
(`code/rebuild/fk_kernel_crosscheck.py`, agreement 0.99-1.02 over 0.5'-5'), so
these numbers are solid as a target for the *assembly*; they share the vertex
table with everything else, which is exactly why the stochastic test is wanted.

## Environment and machine

* Interpreter: `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`
* `export PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src`
  (the callables import `canoes`; without this you get `ModuleNotFoundError`)
* The vertex callable to bind is
  `driver_field_emulators/code/callable_fixed/perm_aware_kappa3_callable.py`;
  set `pa.TABLE_PATH` to
  `driver_field_emulators/products/table_permclosed_cut15360.npz`, set
  `pa._CACHE = None`, then `driver_stats._k3 = pa`. The pattern is in
  `driver_field_emulators/code/figures_corrected/regenerate_mc.py`.
* 96 GB, 28 cores, shared with other sessions. One single-threaded process per
  (gamma, seed) is the right parallelism; cap `OMP_NUM_THREADS=2` and do not
  exceed about 6 concurrent processes. An OOM reboot has happened on this
  machine before.
* **zsh does not word-split unquoted variables.** `set -- $PAIR` silently puts
  the whole string in `$1`. Use explicit 1-indexed arrays and quote every use,
  and verify background jobs with `pgrep -fl` before believing they started.

## What to report back

State plainly which gate you reached. If the estimator works, give the ratio
to the analytic FK with its uncertainty, the sigma_lambda trend, and a
recommendation for what the paper may then claim. If it does not, give the
failure mode and the number that shows it, and say so; a clean negative is a
useful result here and the paper's current wording already accommodates it.

Do not edit the manuscript, its figures, or its generators. Report, and let
the main session decide.
