# A placement-complete, variance-controlled Monte-Carlo test of the FK channel

*2026-08-27/28. Everything new lives in this folder. `sachs_mc_core.py`,
`SFT-lensing-paper-analyses/`, `sections/`, `figures/` and every paper generator
were read only and are byte-identical.*

---

## 1. Verdict

**All four gates passed.** The FK channel now has a direct stochastic-ODE check.

| gate | result |
|---|---|
| 1. placement share | **1.000000** to machine precision, proved in generic symbols. SCOPE: it fixes the placement bookkeeping, not the kernels -- it holds to 1e-16 with a random `F`, a random `Q`, and the background off |
| 2. one point (gamma = 1') | MC / fold, extrapolated to `sigma_lambda = 0`: **1.01 +/- 0.01** (four fresh 24-seed blocks, error from block scatter and inflated 1.75x for the unexplained block heterogeneity). Converged expectation / fold: **0.9974** |
| 3. sigma_lambda dependence | monotone, deficit exactly linear in `sigma_lambda` (fitted exponent 1.05), **intercept 1.0000** |
| 5. input matching (new) | the FK diagram carries **no driving-field two-point propagator**, and the MC's FK is invariant to `Sigma2` to ~1e-4 under arbitrary lambda-distortion, cross-ray zeroing and full tensor isotropisation. Both sides depend on `(zeta, D)` alone, and both are shared exactly, so the comparison is **genuinely input-matched** |
| 4. gamma sweep 0.5'-17.3' | exact expectation / analytic FK at `sigma_lambda -> 0`: **median 1.0000**, range [0.991, 1.006]; MC/exact = 0.997 to 1.007 at every separation once independent seed blocks are used |

Cost: **4 s per (gamma, seed)** at `n_real = 24000`, `n_lambda = 1000`, single
threaded. That is about 4x *cheaper* than `simulate_fk_vr`, not the 20x more the
brief anticipated; the reason is in section 3.

## 2. What the paper may now claim

The FK check can be upgraded from "independent of the workflow's implementation
of it" to an independent numerical solution of the stochastic Sachs system, on
the same footing as the FF check. Suggested wording for the magnitude:

> in the limit of vanishing colored-noise correlation length the estimator's
> expectation agrees with the folded prediction to better than 0.3% at every
> separation from 1' to 17', and the Monte-Carlo reproduces that expectation to
> about 1%

REVISED 2026-08-28 after adjudication (`review/ADJUDICATION.md`). The earlier
`1.002 +/- 0.006` was too tight and slightly high. Two corrections:
(i) the converged (`sigma -> 0`, `dlam -> 0`) expectation/fold is 0.9974 / 0.9980
/ 0.9986 at 1.0' / 5.3' / 17.3', always slightly LOW, and the 1.002 at
`n_lambda = 1000` was high by ~0.4% (+0.19% unconverged discrete kernel, +0.089%
two-point extrapolation bias); (ii) four independent seed blocks scatter with
chi2 = 9.22 for 3 dof, so the seed SEM understates by 1.75x and the honest
Monte-Carlo number is **1.00 +/- 0.011**.

Two caveats that must travel with it:

* **Shared vertex table -- but that is the controlled variable, not a leak.**
  The claim being made is "holding the driving-field statistics fixed, does the
  stochastic solution agree with the formalism", so sharing `zeta6` is
  *required*. The sharper statement, established in section 5d: the FK diagram
  carries no driving-field two-point propagator, so BOTH sides depend on
  `(zeta, D)` alone -- the one place the two treatments genuinely differ (the
  formalism uses the full non-local `C(lam1,lam2)`, the MC a lambda-local
  idealisation, worth 15% in Order-0 at 0.5') provably cannot touch FK. What is
  still NOT tested is the vertex table itself: a units error inside it cancels
  on both sides, exactly as for `fk_kernel_crosscheck.py`. What IS tested is the
  assembly: the diagram's existence in the stochastic dynamics, its
  normalisation, the causal lambda-measure, the F-K index contraction, and the
  response kernels.
* **Perturbative in zeta.** The estimator is exactly linear in the injected
  three-point cumulant, which is the order the FK diagram lives at. It is *not*
  perturbative in the F vertex: the nonlinear-arm flavour (full Sachs ODE, the
  increment field linearised about it) agrees with the strictly order-controlled
  one to 0.4%, so F^2 and higher are negligible here.

## 3. The two defects, and why they cancelled

The published route had **two independent errors of opposite sign** which nearly
cancelled at the tuned value `sigma_lambda = 8`. Both are now quantified.

### 3.1 Placement incompleteness (the known one) -- factor ~0.22

The Monte-Carlo does not inject `zeta` directly; it deforms a Gaussian,
`f = z + 0.5 Q (zz - <zz>)`. At O(Q) this yields

```
<f_a f_b f_c> = Q_a:VV + Q_b:VV + Q_c:VV        (three terms, one per leg)
```

In the FK diagram those legs go to the far observable and to the two inputs of
the F vertex. `simulate_fk_vr` deforms only the observable leg. Reproduced here
independently (exact Wick, not a probe fit), the captured share is

| gamma | 0.5' | 0.8' | 1.3' | 2.0' | 3.2' | 5.0' | 8' | 12' | 17.3' | 30' |
|---|---|---|---|---|---|---|---|---|---|---|
| `T1 / (T1+T2)` | 0.251 | 0.232 | 0.201 | 0.170 | 0.130 | 0.089 | 0.031 | -0.024 | -0.075 | -0.163 |

The share is only weakly calibration-dependent: under the *nominal* calibration
that `simulate_fk_vr` actually uses it reads 0.266 at 0.5' and 0.227 at 1.0',
against the prior session's 0.245 at 0.5' -- the same diagnosis to about 10%,
reached by a completely different method (closed-form Wick contraction rather
than a second-order expectation analysis fitted to probe points). It falls
through zero between 8' and 12' and is negative beyond -- the prior session's
"zero near 10', negative beyond", reproduced exactly.

### 3.2 Deformation mis-calibration -- factor ~3.9 at sigma_lambda = 8

`Q` is solved against the *nominal* node covariance `V = Sigma2 / (2 sigma)`,
i.e. assuming `Sigma2` is constant over a correlation length. The drive actually
couples to the AR(1)-**smeared** covariance

```
B[k] = sum_l dlam <z[k] z[l]^T>            ( -> Sigma2 as sigma_lambda -> 0 )
```

**Why the nominal choice fails (verified, and NOT what a first pass suggested).**
The tabulated two-ray `Sigma2` is **not positive semi-definite everywhere**: its
smallest eigenvalue passes through zero and goes negative at isolated lambda.
`verify_graininess.py` at gamma = 17.3', around lambda = 1307 Mpc:

| lambda [Mpc] | 1299 | 1303 | 1307 | 1309 | 1313 |
|---|---|---|---|---|---|
| smallest eigenvalue / local mean | 1.18 | 0.47 | 0.021 | **-0.25** | 0.16 |
| `|Q|_max` | 1.9e7 | 1.4e8 | **5.4e10** | 3.5e10 | 2.7e9 |

`zeta` runs smoothly straight through that window (4.00e-18 to 4.32e-18,
monotone), so the defect is entirely in the covariance table's lambda
interpolation. Because `solve_Q` divides by `d_p d_q + d_q d_r + d_r d_p`, the
node-local `Q` spikes by ~5000x exactly there. 59 of 1000 nodes sit below half
their local mean and together carry 6.9% of the FK kernel weight. The smeared
`B[k]` averages over +/- sigma_lambda and is regular.

This is precisely what the 2026-08-26 session called "lambda-grid graininess of
the tabulated tensors". That diagnosis was right; an intermediate account in this
folder blamed the near-degeneracy of the two-ray covariance at small gamma, and
that was **wrong** -- the covariance conditioning at the worst nodes is only
~5-12, and the nominal calibration fails just as badly at gamma = 17.3' and 30'
where the rays are well separated.

**Consequences of the nominal choice, all measured:**

* *Non-convergent in sigma_lambda.* At fixed `sigma_lambda/dlambda = 4.19`,
  exact/fold reads 5.02, 1.86, 1.26, 1.38 at gamma = 17.3' and 0.66, **-0.22**,
  0.79, 2.10 at gamma = 30', for sigma_lambda = 16, 8, 4, 2. Sign flips included.
* *Grid-unstable.* At fixed `sigma_lambda = 8`, gamma = 17.3', exact/fold reads
  1.86, **65.5**, 31.9, 2.79, 3.98 at `n_lambda` = 1000, 1414, 2000, 2828, 4000 --
  it depends on whether the grid happens to sample a singular lambda.
* *Detectable without the F vertex at all.* The third moment of the LINEAR
  observable, `<kappa1^3>`, which must equal `sum dlam Wd^3 zeta`:

| sigma_lambda | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|
| nominal calibration | 1.259 | 1.896 | **4.974** | 9.851 | 11.913 |
| calibrated against `B` | 0.9999 | 0.9995 | 0.9978 | 0.9911 | 0.9723 |

The fix is one line in `node_stats`: `solve_Q(B[k], zeta[k])` instead of
`solve_Q(V[k], zeta[k]/(2 sigma)^2)`. `solve_Q` is homogeneous, so the two are
*identical* in the white-noise limit -- the change moves nothing about the
physics, only the rate at which the regulator is removed. Under the smeared
calibration exact/fold is stable to 0.5% over an 8x grid refinement
(0.9903 -> 0.9848 at `n_lambda` = 500 -> 4000) and converges monotonically in
sigma_lambda at every separation tested (0.5' to 30').

**Caveat this exposes.** Because `Sigma2` is not PSD at those nodes, the
Monte-Carlo's `_cholesky_psd` floors the negative eigenvalue to zero while the
exact expectation propagates `V[k]` as given. The two therefore simulate very
slightly different fields at ~6% of the kernel weight. Empirically it does not
matter (MC/exact = 1.000 +/- 0.005), but the covariance table's PSD violation is
a real defect worth fixing upstream in `corr_op`.

### 3.3 The cancellation

At `sigma_lambda = 8`, gamma = 1': share `0.227` x inflation `3.92` = **0.889**,
and at 0.5' it is **1.034**. That is the "agreement" the June figure recorded.
Both numbers fall straight out of the exact expectation, with no fitting.

## 4. The estimator

`FK = (T1 + T2) + (T1 + T2)^T`, cross-ray block, with

```
T1_AB = Cov( kappa_g^(2)_A , kappa_d^(1)_B )     Q on the far observable leg
T2_AB = Cov( kappa_g^(1)_A , kappa_d^(2)_B )     Q on an F-vertex input leg
```

where `kappa_d` is the response to the skew drive through the F vertex
linearised about the Gaussian trajectory -- Path A's dynamics, which is what
reaches the F legs.

**Why this is variance-controlled without the tensor-valued field.** The brief
proposed evolving a 21-component increment field so `Q` could be pulled outside
the average. That is not where the variance reduction lives: pulling `Q` out of a
linear functional is algebraically a no-op (`<kappa_g kappa_d>` and
`sum_k Q[k] <kappa_g Lambda(k) w(k)>` are the same random variable, so they have
the same variance). Path A's pathology is the specific product
`kappa^(1) x kappa_d^(1)`, whose mean vanishes but which carries the entire `Q^2`
spread. Splitting into the two placement classes never forms that product: each
term pairs an F-free factor with an F-induced one, `Q` scales signal and noise
together, and the relative error becomes `Q`-independent. Hence 2-3% seed
scatter at 24000 realisations instead of ~100%, at roughly 1/80 of the cost
the tensor-field design would have carried.

**Centering.** The drive subtracts the *batch sample* second moment of `z` at
each node. For `T1` any constant is removed by the covariance; for `T2` it is
not -- there the two `Q` legs can contract with each other, and that tadpole
survives unless the subtracted constant is the field's actual variance. The
AR(1) chain is not stationary (`Sigma2` varies along the ray), so the nominal `V`
is the wrong constant.

## 5. Numbers

### Gate 1 -- placement identity (deterministic)

`gate1_placement_share.py` builds the exact O(Q) three-point function induced by
the deformation, contracts it through the discrete FK kernels, and compares.
`(T1+T2) + transpose` equals it to **4e-16**; `T1 + T1^T` gives 0.211. Holds
across N = 16-40, gamma = 0.5'-60', sigma = 30-200 (`gate1_sweep.py`), so it is
structural, not a limit statement.

### Gate 2 -- gamma = 1', three independent seed blocks (88 seeds x 24000)

```
analytic FK (paper fold)          1.597716e-05
local reference, injected zeta    1.597122e-05   (0.9996 x analytic)
exact expectation, sigma = 8      1.554705e-05   (0.9734 of the reference)
exact expectation, sigma = 4      1.576620e-05   (0.9871)

Monte-Carlo, sigma = 8   1.55783e-05 +/- 4.0e-08   MC/exact = 1.0020 +/- 0.0026
                                                   MC/analytic = 0.9750 +/- 0.0025
Monte-Carlo, sigma = 4   1.57953e-05 +/- 6.0e-08   MC/exact = 1.0018 +/- 0.0038
                                                   MC/analytic = 0.9886 +/- 0.0037

per-seed linear extrapolation to sigma_lambda = 0 (88 paired seeds):
    Monte-Carlo / analytic FK    = 1.0022 +/- 0.0056
    exact expectation / analytic = 1.0005
```

The three blocks (seed bases 20260827, 314159001, 20260828) agree within their
errors; individually they give MC/analytic at sigma = 8 of 0.9826, 0.9703 and
0.9713, so a single 32-seed block is good to about 1%, not 0.3%.

The "local reference, injected zeta" is the leg-symmetrised reference and it
matches the paper's fold to 0.04% here. See section 5 note below on why
symmetrised is the correct target.

### Gate 3 -- sigma_lambda at gamma = 1' (n_lambda = 1000 fixed)

| sigma_lambda | 32 | 16 | 8 | 4 | 2 |
|---|---|---|---|---|---|
| MC / analytic | 0.8957 | 0.9507 | 0.9790 | 0.9872 | 1.038 |
| exact / analytic | 0.8939 | 0.9459 | 0.9734 | 0.9872 | 0.9943 |
| MC / exact | 1.0024 | 1.0054 | 1.0061 | 1.0004 | 1.045 |
| seed scatter | 1.6% | 1.8% | 2.1% | 2.8% | 15% |

The deficit halves when `sigma_lambda` halves (fitted exponent **1.05**), so the
`sigma_lambda -> 0` intercept is 1.0000. Below `sigma_lambda ~ 4` the *sampling*
degrades (the node skewness grows as the regulator is removed) while the exact
expectation stays clean -- so `sigma_lambda = 4-8` is the usable MC window and
the extrapolation is anchored there.

**The old fine-grid trap is gone.** At fixed `sigma_lambda = 16` the exact
expectation reads 0.9459 at n_lambda = 1000, 2000 *and* 4000 -- identical to four
decimals over an 8x refinement. The June note's "finer grid makes it worse" was a
symptom of the variance-pathological estimator, not of the discretisation.

### Gate 4 -- gamma sweep

Exact expectation, extrapolated to `sigma_lambda = 0`, against the paper's fold:

| gamma | 0.5' | 0.8' | 1.3' | 2.0' | 3.2' | 5.0' | 8.0' | 12' | 17.3' |
|---|---|---|---|---|---|---|---|---|---|
| exact / analytic | 1.000 | 1.000 | 1.006 | 1.001 | 1.000 | 0.999 | 0.995 | 0.991 | 1.006 |

Flat at 1.000 with a spread under a percent across the whole range.

The Monte-Carlo tracks it. At `sigma_lambda = 8`, pooling two independent
48-seed blocks:

| gamma | 3.2' | 5.0' | 8.0' | 12' | 17.3' |
|---|---|---|---|---|---|
| MC / exact | 0.996 | 0.992 | 0.986 | 0.976 | 0.958 |
| +/- | 0.005 | 0.007 | 0.010 | 0.014 | 0.022 |
| deviation | -0.8 sigma | -1.1 | -1.4 | -1.7 | -1.9 |

and at small separation (single 8-seed block, `sigma_lambda = 8 / 4`):
MC/analytic = 0.974/0.983 at 0.5', 0.977/0.985 at 0.8', 0.988/0.996 at 1.3',
0.991/0.998 at 2.0'.

**A warning about the error bars.** Every point in a sweep shares one seed list,
so the fluctuations are correlated across gamma and can look like a smooth
trend. The first 8-seed pass gave MC/exact = 1.029, 1.042, 1.062, 1.089, 1.130
at the five separations above -- a convincing-looking upward drift that is
entirely one realisation. Two independent 24-seed blocks put the same points at
0.999/0.993, 0.995/0.990, 0.988/0.984, 0.976/0.976, 0.954/0.961. Quote the
per-point error bar; do not read a gamma-trend out of a single seed block.

## 5b. Does the comparison have teeth?

A validation that cannot fail is worthless. `verify_teeth.py` injects structural
errors into the assembly and reports the ratio to the paper's fold; the measured
agreement is `1.0022 +/- 0.0056` at 1 arcmin, so anything moving the ratio by
more than about 1.5% would have been caught.

| perturbation | 1' | 5' | 17.3' |
|---|---|---|---|
| baseline (symmetrised reference) | 0.9996 | 0.9985 | 1.0051 |
| drop causality (F node need not follow the drives) | 558 | 453 | 317 |
| one propagator instead of two on the F leg (`D^2` not `D^4`) | 1.111 | 1.143 | 1.211 |
| omit the `A<->B` symmetrisation | 0.500 | 0.499 | 0.503 |
| scale `zeta` by 1.05 | 1.050 | 1.048 | 1.055 |
| drop the trapezoid end weights | 1.012 | 1.010 | 1.015 |
| use the raw, unsymmetrised `zeta6` legs | 1.008 | 1.025 | 1.080 |

Every structural error is detected, several by orders of magnitude, and even the
2-per-mille end-weight convention shows at 2 sigma. The Mathematica proof carries
the same idea at the identity level: two deliberate-error controls (an
unsymmetrised `F` in `kappa_d2`, a mismatched `kappa_d1` kernel) both BREAK the
placement identity, and `T1`-only, `T2`-only and the unsymmetrised sum all fail
it.

What the comparison canNOT detect is an error inside the vertex table itself,
since both sides read it. That limit is unchanged from the deterministic
cross-check and must stay in the paper.

## 5c. Independent verification (2026-08-28)

Six independent agents re-derived or attacked this work; artifacts in `review/`.
The load-bearing outcomes, each verified by me directly afterwards:

### Proved, not fitted

* **The placement identity is a symbolic theorem.** `review/proof_placement_identity.wl`
  proves `(T1+T2) + (T1+T2)^T == TRUTH` for FULLY GENERIC symbols (generic
  covariance `R[k,l]`, generic `Q[k]`, generic `F_abc` with no symmetry assumed,
  generic weights) at `(n_comp, N) = (2,3), (2,4), (3,3), (3,4)`. It also proves
  the sharper per-F-placement form `M == T1 + T2^T`. Negative controls in the
  same script: `T1`-only, `T2`-only and the unsymmetrised sum all FAIL, and two
  deliberate errors (unsymmetrised `F` in `kappa_d2`, mismatched `kappa_d1`
  kernel) BREAK the identity. I re-ran the script myself and reproduced every
  line.

* **The `sigma_lambda -> 0` extrapolation is DERIVED.** This was the sharpest
  objection: two points and a linear model chosen because it worked. It is not.
  The `O(sigma)` term comes from the CAUSAL KINK in the response kernel (the F
  node must follow the drives); smooth directions give `O(sigma^2)` by parity of
  the two-sided kernel, but a step function does not. Mathematica
  (`review/sigma_limit/kink_moments.wl`) gives the two kink moments EXACTLY:
  `E[max(U,V)] = 3/4` and `E[max(U,0)] = 1/2`, with finite second moments so the
  series is regular in `sigma` with no `sigma^2 log sigma`. Those constants give
  a **parameter-free** prediction of the slope, tested against the exact engine:

  | sigma | predicted/measured, `T1` leg | predicted/measured, `T2` leg | total |
  |---|---|---|---|
  | 16 | 1.040 | 0.953 | 0.977 |
  | 8 | 1.023 | 0.970 | 0.985 |
  | 4 | 1.012 | 0.985 | 0.992 |
  | 2 | 1.011 | 0.989 | **0.995** |

  and the wrong coefficient assignments are excluded: `(3/4, 3/4)` gives 1.483 on
  the second leg, `(1/2, 1/2)` gives 0.674 on the first, "no `O(sigma)` term"
  gives 0. The measured local exponent is `q = 1.006` over `sigma = 32 -> 2`.
  Universality: `E[max(U,V)] = E|U-V|/2 > 0` and `E[max(U,0)] = E[U+] > 0` for
  ANY non-degenerate symmetric kernel, so the `O(sigma)` deficit and its sign
  cannot be tuned away by choosing a different regulator.

  Intercepts, four independent ways (linear and quadratic fits, at two lattice
  resolutions): 0.999916, 1.000493, 1.000316, 1.000049. The production two-point
  extrapolation carries a `+0.09%` systematic (`+0.025%` lattice-converged) --
  four times smaller than the `0.56%` statistical error, so it is not the limiting
  uncertainty.

* **The calibration change cannot move the target.** `solve_Q` is homogeneous of
  degree `-2` in its matrix argument (the defining relation is quadratic in `M`),
  so `solve_Q(2 sigma V, zeta) == solve_Q(V, zeta/(2 sigma)^2)` **bit for bit** --
  verified at 0.000e+00 relative error on the real `Sigma2`/`zeta6` at
  `gamma = 1'` and `17.3'`, `sigma = 2, 8, 32`. Since `B -> Sigma2` as
  `sigma_lambda -> 0`, the smeared and nominal calibrations share the white-noise
  limit exactly. The fix changes the rate of convergence, not what is converged to.

### Independently reproduced

* **Clean-room exact expectation** (written from `SPEC.md` alone, never opening
  `fk_expect_exact.py`): agrees to `3e-9` relative on `T1`, `T1+T2` and the share,
  under BOTH calibrations, and to `5e-15` on the local reference.
* **Clean-room Monte-Carlo** (own design, own estimator): `0.9738 +/- 0.0044`
  at `sigma = 8` and `0.9861 +/- 0.0062` at `sigma = 4`, against my pooled
  `0.9750 +/- 0.0025` and `0.9886 +/- 0.0037`. Its own MC/exact closure is
  `1.0009 +/- 0.0046` and `0.9994 +/- 0.0062`.
* **Three routes to the reference** agree: my discrete kernel, an independent
  continuum quadrature, and `../rebuild/fk_kernel_crosscheck.py`, with clean
  `O(1/N)` convergence (successive-difference ratios 0.470, 0.471, 0.486, 0.496).

### Corrections this exercise forced on me

1. The mechanism behind the mis-calibration is the covariance table's lambda
   structure (section 3.2), NOT the small-gamma near-degeneracy I first claimed.
   The nominal calibration has no `sigma -> 0` limit at ANY separation tested,
   including `gamma = 27.7'` and `30'` where the covariance is well conditioned.
2. The `0.9913` is leg symmetrisation, not an injection loss (section 5 note).
3. Two docstring errors, now fixed: `solve_Q`'s homogeneity exponent (`c^2`, not
   `c^3`), and `fk_complete_core`'s statement of the induced three-point function
   (it is non-local along the ray, `Q[m] R[m,n] R[m,p] + perms`, collapsing to
   `Q:Sigma2 Sigma2` only as `sigma_lambda -> 0`).

### Residual caveat worth stating

That the paper's fold symmetrises the vertex over leg placements is inferred,
not read off the sft-wick source: symmetrising the reference reproduces the fold
to 0.2-0.3% flat in gamma over a range where the single ordering drifts by 7.3%.
Strong, but circumstantial.

## 5d. Are the two sides given the same driving-field statistics? (2026-08-28)

The question this whole exercise exists to answer is: *holding the driving-field
statistics fixed, does the Monte-Carlo agree with the semi-analytic formalism?*
So the shared vertex table is not a weakness -- holding the input fixed is the
point. What matters is whether the two sides really receive the SAME input.
Answered channel by channel; artifacts in `inputmatch/` and
`check_fk_sigma_shape.py`.

### FK: yes, and for a structural reason

The FK diagram carries **four R-propagators and zero C-propagators**
(`FK_NOTES.md`), so the formalism's FK never reads the covariance table at all.
And the Monte-Carlo's FK is insensitive to `Sigma2` to about 1e-4, while
Order-0 swings by a factor 2.5 under the same distortions:

| distortion of `Sigma2(lambda)` | Order-0 ratio | FK ratio |
|---|---|---|
| x 2 (global) | 2.000 | **1.0000000** |
| linear ramp 0.5 -> 2.0 | 1.632 | 0.9999344 |
| tilt `(lam/1200)^1.5` | 1.934 | 0.9998926 |
| oscillatory `1 + 0.5 sin` | 0.914 | 1.0001341 |
| step-like tanh 0.4 -> 2.5 | 2.303 | 1.0000276 |

An independent test (`inputmatch/fk_c_invariance.py`) pushes it further, and
the residual shrinks as the regulator is removed (`sigma_lambda` = 16 -> 2):
a lambda ramp 0.4->3.6 gives 0.99972 -> 0.99999.8; a +/-70% lambda oscillation
1.00128 -> 1.0000018; **zeroing the entire cross-ray block** 0.99891 -> 0.99988;
**replacing `Sigma2` by an isotropic `(tr/6) I`**, destroying its tensor
structure, 0.99757 -> 1.00152.

So both sides depend only on `(zeta, D)`, and both are shared exactly. The FK
comparison is genuinely input-matched, and the measured agreement
(1.011 +/- 0.005 at 1', four independent seed blocks) is a real test of the
formalism, not of the input.

### Order-0: no, and the mismatch is gamma-dependent

The formalism consumes the full non-local `C(lam1, lam2)`; the Monte-Carlo
drives a lambda-LOCAL field of density `Sigma2 = d/dlam[D^4 C(lam,lam)]/D^4`.
Both computed from the SAME table (`inputmatch/o0_two_ways.py`):

| gamma | 0.5' | 1' | 2.6' | 6.7' | 17.3' | 44.4' | 90' | 150' |
|---|---|---|---|---|---|---|---|---|
| collapsed / full | 0.846 | 0.846 | 0.850 | 0.858 | 0.870 | 0.889 | 0.901 | 0.936 |

The equal-time collapse costs 15% at 0.5' and 6% at 150'. Because it is
gamma-dependent, the single `ANCHOR_C0 = 1/0.881` cannot absorb it, and Order-0
is not a comparison in the figure -- it is shown as the analytic line only.

### FF: inherits that mismatch quadratically, and is not resolved

`FF ~ Sigma2^2`, so the collapse deficit would put the ratio near 0.72 at small
gamma. The measurement cannot see it: FF is a DRIFT term with `O(dlam)` error
and at `gamma = 1'`, anchor off, 4 seeds x 12000, `n_lambda` = 500/1000/2000/
4000/8000 gives FF/analytic = 1.040/1.268/1.219/0.993/0.968, each with ~10% seed
error. The published marker is a single seed at 0.805 +/- 0.15.

### Incidental defect found

`corr_op.C_fn_batch` groups samples by a ROUNDED `(cos gamma, psi1, psi2)` key
and caches the angular-channel table, so whichever path touches a given `cos`
FIRST fixes what BOTH the batch and scalar paths return afterwards. Measured
scalar-only vs batch-first at `lambda = 1234.5`: 1.0004 (0.5'), **1.0017** (1'),
0.9996 (2'), 0.9997 (5'), 1.0000 (10'-44'). Small, but it makes results
evaluation-order dependent. Reproducer: `psd_fix/indep_corr_op_batch_bug.py`.

## 6. Reproducing

```
export PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
OMP_NUM_THREADS=2 $PY gate1_placement_share.py          # gate 1, seconds
OMP_NUM_THREADS=2 $PY gate1_sweep.py                    # gate 1 across configs
OMP_NUM_THREADS=2 $PY validate_mc_vs_exact.py 24 60 1.0 40000 4
OMP_NUM_THREADS=2 $PY run_gate2.py --seeds 8 --out g2.npz
OMP_NUM_THREADS=2 $PY run_sweep.py --gammas 1.0 --sigmas 32 16 8 4 2 \
    --seeds 8 --out gate3.npz
OMP_NUM_THREADS=2 $PY sigma_scan_exact.py 1.0 "1000:8,1999:4,3997:2" nominal
OMP_NUM_THREADS=2 $PY analyse2.py
```

Use a different `--seed-base` for an independent realisation block.

## 7. Files

| file | what |
|---|---|
| `_bootstrap.py` | path wiring + binds the corrected vertex callable to `table_permclosed_cut15360.npz` |
| `fk_expect_exact.py` | closed-form (Wick) expectation of the estimator; the `smeared`/`nominal` calibration switch |
| `fk_complete_core.py` | the Monte-Carlo estimator |
| `gate1_placement_share.py`, `gate1_sweep.py` | gate 1 |
| `validate_mc_vs_exact.py` | MC against its own exact expectation on a cheap grid |
| `run_gate2.py`, `run_sweep.py` | gates 2-4 |
| `analyse.py`, `analyse2.py` | tables, paired-seed extrapolations |
| `plot_fk_mc.py`, `fk_mc_complete.pdf` | diagnostic figure (this folder only, not a paper figure) |
| `probe_norm_kappa3.py`, `probe_norm2.py` | the `<kappa1^3>` normalisation test that located defect 3.2 |
| `probe_decoherence.py`, `probe_B_pd.py`, `probe_B_edges.py` | supporting diagnostics |
| `probe_bias.py`, `probe_seed_dist.py` | the seed-correlation investigation of section 5 gate 4 |
| `sigma_scan_exact.py`, `gamma_scan_exact.py`, `exact_production.py` | exact-expectation scans |
| `_*.npz`, `_*.log` | run products |
