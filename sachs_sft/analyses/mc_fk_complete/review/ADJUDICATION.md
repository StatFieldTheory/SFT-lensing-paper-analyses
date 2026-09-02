# Adjudication of the FK Monte-Carlo validation exercise

*Adjudicator's report, 2026-08-28. I did not do the work under review. Six
first-principles claims and six adversarial attacks were submitted; this weighs
them. Everything I ran myself is under `review/adjudication/`; nothing outside
`review/` was created or modified.*

---

## 0. What I verified myself

I did not adjudicate on the reports. Six independent checks were run here, four
of them on the load-bearing numbers. All used the PyCCL env, `OMP_NUM_THREADS=2`,
at most two processes, `n_real = 24000`, `batch = 2000`.

| # | script | what it settles |
|---|---|---|
| V1 | `gate1_placement_share.py` (unmodified) | the placement identity and its 0.211/0.789 split |
| V2 | `adjudication/v1_sigma_limit.py` | the whole `sigma_lambda` ladder, incl. switching the regulator OFF |
| V3 | `adjudication/v2_reference_convergence.py` | leg symmetrisation, and `dlam -> 0` of the reference |
| V4 | `adjudication/v3_gate1_blindness.py` | whether gate 1 can fail |
| V5 | `adjudication/v4_kappa1_cubed.py` | the F-free calibration discriminator |
| V6 | `run_sweep.py` fresh seed base + `adjudication/v5_blocks.py` | a fourth independent MC block |

**V1.** Reproduced exactly: truth `1.589134e-05`; `T1+T1^T` = `3.353740e-06`
(0.2110); `T2+T2^T` = `1.253760e-05` (0.7890); `(T1+T2)+transpose` matches the
truth to `4.015e-16` on the full 6x6.

**V2.** `gamma = 1'`, `N = 1000`, smeared calibration. Exact expectation and its
ratio to the interpolated fold `1.5977160e-05`:

| sigma | exact | exact/fold | exact/ref_inj | ref_inj/fold | injection fidelity |
|---|---|---|---|---|---|
| 16 | 1.5107761e-05 | 0.9455849 | 0.9459367 | 0.999628 | 8.7e-12 |
| 8 | 1.5547050e-05 | 0.9730797 | 0.9734418 | 0.999628 | 8.6e-12 |
| 4 | 1.5766214e-05 | 0.9867970 | 0.9871642 | 0.999628 | 8.6e-12 |
| 2 | 1.5880611e-05 | 0.9939570 | 0.9942786 | 0.999677 | 1.2e-05 |
| 1 | 1.5917482e-05 | 0.9962648 | 0.9979906 | 0.998271 | 8.9e-05 |
| 0.5 | 1.5944728e-05 | 0.9979701 | 0.9996990 | 0.998271 | 8.9e-05 |
| 0.25 | 1.5949422e-05 | 0.9982639 | **0.9999933** | 0.998271 | 8.9e-05 |

`rho = 4.83e-04` at `sigma = 0.25`: the regulator is off. NOTES' gate-2 exact
numbers (`1.554705e-05`, `1.576620e-05`) and its `0.9734 / 0.9871` are
reproduced to every digit printed. Published two-point extrapolation
`2r(4) - r(8) = 1.000514`.

**V3.** The cumulant the deformation injects is *exactly* the 6-fold symmetrised
tabulated one: `ref_inj / ref_sym = 1.000000000`. The raw table's leg asymmetry
is `ref_raw/ref_sym` = **1.0088** (1.015'), **1.0284** (5.304'), **1.0762**
(17.276'). Convergence of the symmetrised discrete reference against the fold, at
exact fold nodes:

| gamma | N=500 | 1000 | 2000 | 4000 | 8000 | 1/N Richardson |
|---|---|---|---|---|---|---|
| 1.01546' | 1.00146 | 0.99936 | 0.99837 | 0.99789 | 0.99766 | **0.99742** |
| 5.30409' | 1.00261 | 1.00022 | 0.99908 | 0.99852 | 0.99825 | **0.99797** |
| 17.27554' | 1.00418 | 1.00133 | 0.99995 | 0.99927 | 0.99894 | **0.99860** |

**V4.** The gate-1 identity, re-run inside deliberately wrong but
self-consistent systems (`resp`/`Wd`/`Hd` rebuilt from a substituted `D`;
`F`, `Q`, `A` monkeypatched):

| world | truth (kk) | identity error |
|---|---|---|
| true system | 1.5891e-05 | 4.015e-16 |
| propagator exponent 1 | 1.9340e-05 | 6.750e-16 |
| propagator exponent 6 | 1.9867e-05 | 4.959e-16 |
| D = 1 (background off) | 3.2442e-05 | 6.720e-16 |
| F -> random (6,6,6) | 1.9455e-05 | 3.773e-16 |
| Q -> random symmetric | -5.1364e-12 | 4.407e-16 |
| A -> random SPD | -1.5961e+26 | 2.986e-16 |
| all wrong at once | +1.9559e+18 | 3.038e-16 |

**V5.** `<kappa1^3>` at O(Q), no F vertex, no FK kernel, no fold, against
`sum_j dlam Wd_j^3 zeta_sym`, `N = 1000`, channels (000)/(003)/(033)/(011):

| calibration | sigma | ratios |
|---|---|---|
| nominal, 1' | 8 | 5.0986 / 4.9677 / 4.9677 / 3.0867 |
| nominal, 1' | 2 | 1.2602 / 1.2577 / 1.2577 / 1.2129 |
| smeared, 1' | 8 | 0.9976 / 0.9978 / 0.9978 / 0.9975 |
| smeared, 1' | 2 | 0.9998 / 0.9999 / 0.9999 / 0.9998 |
| nominal, 17.3' | 8 | 5.3860 / 1.5786 / 1.5786 / 3.3254 |
| smeared, 17.3' | 8 | 0.9976 / 0.9989 / 0.9989 / 0.9976 |

**V6.** A fourth independent seed block (base 5150271, 16 seeds x 24000,
`gamma = 1'`, `N = 1000`). Paired-seed `sigma -> 0` extrapolation of
MC/analytic, across every block I could read:

| block | n | MC/ana (8) | MC/ana (4) | extrapolated |
|---|---|---|---|---|
| **adjudicator, base 5150271** | 16 | 0.9605 +/- 0.0061 | 0.9670 +/- 0.0081 | **0.9735 +/- 0.0112** |
| existing `_gate3_sigma` | 8 | 0.9790 +/- 0.0073 | 0.9872 +/- 0.0098 | 0.9953 +/- 0.0150 |
| existing `_papergrid` | 24 | 0.9713 +/- 0.0053 | 0.9936 +/- 0.0080 | 1.0158 +/- 0.0119 |
| `adv_sigma` block D | 24 | 0.9742 +/- 0.0062 | 0.9960 +/- 0.0093 | 1.0178 +/- 0.0133 |

Inverse-variance mean **0.9992 +/- 0.0063**, **chi2 = 9.22 for 3 dof**
(p ~ 0.027), error scale factor **1.75**.

---

## 1. Which claims are proved, which are supported, which are unproven

### PROVED

**(P1) The placement identity.** `(T1+T2) + (T1+T2)^T == TRUTH`, and the sharper
`M == T1 + T2^T`, hold identically in generic symbols. The Mathematica proof is
real work: a generic Isserlis engine, generic `R[k,l]` with only
`R[k,l] = R[l,k]^T` imposed, generic `Q`, generic `F` with *no* symmetry, generic
weights, symbolic `dlam`, with the FK kernel read off the recursion rather than
assumed. Five sizes closed at exact zero, and the index argument that carries it
to general `(n_comp, N)` is stated and is sound: it uses only Q's last-two-index
symmetry, `C[k] = R[k,k]`, the two F-vertex legs sharing one time argument and
one kernel, and `kappa_d` reusing `kappa`'s kernels. Negative controls fail as
they must. I did not re-run the `.wl`, but V1 and V4 confirm its numerical
shadow, including the leg-by-leg split.

**Scope, which is essential and is currently mis-stated in `PAPER_PROPOSAL.md`.**
V4 settles this: the identity holds at 3-7e-16 in a world where the background is
switched off, where `F` is a random tensor, where `Q` is random, where the
covariance is random, and where the truth is `1.96e+18`. It is a theorem about
`expectation()` being a faithful Wick evaluation of *whatever* `(grid, A, Q, F)`
it is handed. It certifies **how the deformation's three legs are distributed**,
not **which diagram the kernels encode**. The proposal's bullet "the estimator
captures the whole vertex, by a Mathematica identity" reads as the second and
must be qualified.

**(P2) `solve_Q` is homogeneous of degree -2 in its matrix argument**, hence
`solve_Q(2 sigma V, zeta) == solve_Q(V, zeta/(2 sigma)^2)` bit for bit, and the
nominal and smeared calibrations share the white-noise limit exactly. Proved
step-by-step through the routine (the ridge scales with the matrix; the `rcond`
mask is relative; `U` enters an even number of times so sign ambiguity cancels),
and verified at `0.000e+00` on the real tensors by two agents independently. This
is what makes the calibration change non-negotiable rather than a tuning knob.

**(P3) The `sigma_lambda -> 0` limit of the estimator's expectation is the local
FK integral, and nothing survives at `sigma = 0`.** This is proved, not
asymptotic: the estimator's expectation is *exactly* the FK diagram of the
colored driving field with a `sigma`-independent kernel, so all `sigma`-dependence
sits in the drive's induced third-cumulant density, whose weak limit is the
injected cumulant. The enabling identity was verified to `6e-16..8e-17` over
eight configurations. I confirm the consequence directly: at `sigma = 0.25`
(`rho = 4.8e-04`) `exact/ref_inj = 0.9999933`, and the approach is monotone from
0.9459.

**(P4) The leading correction is O(sigma^1), from the causal kink.** Proved: the
smooth directions have zero first moment against the AR(1) kernel, but the
product of the two retarded legs, `theta(max(m,n) < j)`, has a kink whose mean is
not zero. The two coefficients are exact rationals, `E[max(U,V)] = 3/4` and
`E[max(U,0)] = 1/2`, with finite second moments, so the expansion is regular with
no `sigma^2 log sigma`. The universality argument (`E|U-V|/2 > 0`, `E[U+] > 0` for
any non-degenerate symmetric kernel) is correct and means the deficit's sign
cannot be regulator-tuned away. **What is proved is the order and the
coefficients; that the measured slope *is* that term is supported, not proved**
(agreement 0.5% at 1', 1.0% at 5', residual shrinking with sigma). The
normalisation-free discriminator, the ratio of the two placement slopes
(measured 1.4114 vs derived 1.4433, equal-coefficient 0.9622, reversed 0.6415),
is the strongest part of that support.

**(P5) The smeared calibration is forced, not fitted.** The derivation is
F-free: `<kappa1^3>` at O(Q) equals
`sum_j dlam Wd_j [Q[j]:Y[j]Y[j] + 2 perms]` with
`Y[j] = sum_p dlam Wd_p R[j,p]`; since `Wd` is slowly varying on the scale
`sigma` (verified, `1e-4..7e-3`), `Y[j] -> Wd_j B[j]`, and because `Wd_j^3 > 0`
and the identity must hold for any window, the condition is node-wise
`cum3_from_Q(Q[j], B[j]) = zeta[j]`, whose unique solution is
`Q = solve_Q(B, zeta)`. No F vertex, no FK kernel, no fold enters. I verified the
consequence myself (V5): smeared reads 0.9975-0.9999 in every channel at both
separations; nominal reads 1.26 to 5.39. Strictly this is a derivation with one
controlled approximation (`Y -> Wd B`), not a theorem, so I record it as **proved
modulo a numerically bounded step**.

### SUPPORTED NUMERICALLY

* **The clean-room exact expectation.** Twelve quantities agreeing to 3e-9 and
  the local reference to 5e-15, from a re-derivation written against `SPEC.md`
  alone. This is a strong implementation check. It is not an independent physics
  check: both sides consume the same tables.
* **The clean-room Monte-Carlo.** Independent design (pathwise derivative in the
  deformation amplitude rather than the author's split), landing within errors,
  with its own MC/exact closure at 1.0009 +/- 0.0046. Its independent
  rediscovery of the 0.22 placement share and the 3.91 nominal inflation, from a
  different construction, is the most valuable part.
* **Three routes to the reference kernel** with clean O(1/N) convergence; my V3
  is a fourth and agrees.
* **Gate 4 across gamma**, and its converged form (V3).
* **The identification of the measured `sigma`-slope with the derived kink term**
  (see P4).
* **MC/exact closure** at 1.00 to about 1% (V6 pooled, error inflated).

### UNPROVEN

* **That the paper's fold symmetrises the vertex over leg placements.** This is
  inferred from a fingerprint, not read from the `sft-wick` source. The
  fingerprint is strong and I reproduced it (V3): the symmetrised reference sits
  at 0.997-0.999 of the fold flat across 1'-17', while the raw ordering drifts
  1.0088 -> 1.0762 over the same range. But it remains circumstantial, and it
  matters, because the injected cumulant is the symmetrised one *exactly*
  (`ref_inj/ref_sym = 1.000000000`), so the entire comparison is blind to the
  table's leg asymmetry, which is 0.9% at 1' and 7.6% at 17.3'.
* **The vertex table's units and normalisation.** Shared; cancels. Already
  stated.
* **Three further shared inputs the caveat does not name**: the callable that
  *reads* the table (the fold's vertex module is a verbatim copy of
  `perm_aware_kappa3_callable.py`), the background `D(lambda)`
  (`background.py` wraps the same `D_callable.py` the sft-wick config names), and
  the `F` tensor (`F_tensor.npy` identical to `sachs_mc_core._F` at 0.0).
  Quantified blind spots: a `(1+z)^-4` radial-measure error is a factor 11.5 and
  invisible; a 1% tilt in `D` biases FK by 0.8% and is invisible; a pure rescaling
  of `D` is invisible to machine precision.
* **That the Monte-Carlo estimates its own expectation without bias below 0.1%.**
  A measured, one-sided `O(1/batch_size)` bias of `-0.11% +/- 0.07%` at
  `batch = 2000` is real and unbudgeted (`cov()` normalises by `mb` with batch
  sample means; `c_node` is a batch sample second moment).
* **Why the seed blocks are heterogeneous.** Two independent detections now
  (the statistics agent's chi2 = 11.57/5 on six blocks; my chi2 = 9.22/3 on four,
  V6). No mechanism found; RNG stream correlation is ruled out, and a careful
  93-dof variance decomposition finds no excess *within*-seed variance, which
  means the excess must be a per-block common offset that nobody has explained.
* **Everything outside the kappa-kappa entry.** Only `(a,b) = (0,0)` was tested.
  The `F` contraction reaches only row 0: `F_101 = F_202 = -2` are never touched,
  and halving them changes nothing. Any statement covering `xi_+`, `xi_-` or
  `kappa-gamma` is unsupported.

---

## 2. Which attacks survived

Four of the six attacks were reported refuted and one was reported as landing.
My reading differs from the submitted verdicts in two places.

### SURVIVED (real residual doubt)

**A. Gate 1 is blind to the physics — the strongest surviving objection.**
Reported `refuted=false`, and I confirm it, having reproduced the blindness
myself (V4). The attack's *strong* form fails: the identity has teeth against
wrong *estimators* (T1-only 0.2110, T2-only 0.7890, no `A<->B` symmetrisation
0.500, off-by-one causality 1.120, `R -> delta` 0.490), and it catches the
actual historical bug at 0.211. But it is exactly invariant under wrong
*physics*. Consequence: gate 1 is not evidence that the estimator computes the
paper's FK diagram, and the proposal's phrasing must change.

**B. MC/exact is not an independent gate.** The narrow form of the shared-bug
attack is correct and was demonstrated: in a world whose response exponent is 1
instead of 2 and whose FK is 16% wrong, MC/exact still reads 1.0012 +/- 0.0371.
MC/exact is a sampling-vs-Wick check on one shared model. Presenting four gates
as a chain of independent confirmations overstates the evidence; the
non-tautological load sits entirely on the comparison to the `sft-wick` fold,
plus the external anchors (the Order-0 comparison against analysis-3 over a 290x
gamma lever arm, which pins the response exponent to p = 2 with a spread of
0.029 against 1.29 and 1.82 for p = 1 and 3).

**C. The converged agreement is 0.997-0.999, not 1.002.** This is the residual
from the extrapolation and circularity attacks, and my V2+V3 settle it
independently. Decompose

```
exact/fold  =  (exact/ref_inj)  x  (ref_inj/fold)
```

The first factor goes to 1 exactly as `sigma -> 0` (proved; measured 0.9999933).
The second is `sigma`-independent while the injection is exact, equals the
symmetrised discrete reference over the fold, and **drifts with the lattice**:
0.99936 at `N = 1000` but 0.99742 in the `dlam -> 0` limit at `gamma = 1.015'`
(0.99797 at 5.30', 0.99860 at 17.28'). So the honest joint
(`sigma -> 0`, `dlam -> 0`) statement is **0.9974 / 0.9980 / 0.9986**, always
slightly low, monotone in gamma. The published `1.002` at `N = 1000` is high by
about 0.4%, of which roughly `+0.19%` is the discrete kernel not being converged
and `+0.089%` is the two-point extrapolation model's own bias.

I also record a correction to the *attacker's* number here. `adv_sigma` reports
the direct switch-off limit as 0.998264 and treats it as the truth. It is not:
at `sigma <= 2` the injection degrades (V2, fidelity 8.6e-12 -> 8.9e-05) as `B`
loses positivity and `solve_Q`'s `rcond` mask engages, and `ref_inj/fold` steps
down by 0.136%. So the direct switch-off is biased **low** by that step, and the
published extrapolation is biased **high** by 0.089%. The clean number is neither:
it is `ref_sym/fold` at converged `dlam`.

**D. The error bar is optimistic by 1.5-1.75x.** Reported by the statistics
agent and independently reproduced by me on a fresh block (V6): four independent
blocks of the extrapolated quantity give chi2 = 9.22 for 3 dof, scale 1.75, and
my own block is the lowest at 0.9735 +/- 0.0112 against another at
1.0178 +/- 0.0133. The attack's stated *mechanism* (per-realisation kurtosis ~92)
is genuinely dead — the batch-level distribution is essentially Gaussian
(excess kurtosis +0.42, Delta AIC +51.7 against a Student-t, coverage 0.975 on
disjoint groups). But the *conclusion* survives on different evidence. A quoted
`+/- 0.006` from 88 seeds is not defensible; `+/- 0.011` is.

**E. The paper text quotes three different numbers for one quantity.**
`PAPER_PROPOSAL.md` line 20 and 205: `1.002 +/- 0.006`. Line 156 (the proposed
LaTeX): `0.997 +/- 0.006`, which no agent and no pooling I tried can reproduce.
Line 291 (the marker table): `1.016 +/- 0.012`, which is one 24-seed block alone
and which I reproduce exactly as 1.0158 +/- 0.0119. This must be resolved before
anything reaches the manuscript.

**F. An undocumented requirement on the extrapolation protocol.** From the
calibration attack, and it is the most useful thing that attack produced: the
fixed-`sigma/dlam` ladder is only valid if the smearing sum is a
second-order-accurate quadrature. A naive one-sided smearing produces a
beautifully convergent-looking sequence (0.8123, 0.8080, 0.8037, 0.8007) whose
linear extrapolation is 20% wrong, with nothing in the sequence to warn you. It
is exposed only by refining `dlam` at fixed `sigma`, or by the F-free test.

### DID NOT SURVIVE

* **"The calibration is a fudge."** Dead. Thirteen alternative calibrations were
  built and every one is rejected by the F-free `<kappa1^3>` criterion, which
  cannot have been tuned to FK; and each rejected calibration's FK error is
  *predicted* by its F-free score (causal2 0.816 -> 0.804; trunc3.0 1.091 ->
  1.075; geo 1.089 -> 1.071; nominal 2.049 -> 1.748). The legitimate family
  spreads by 0.4% in the limit while differing by 5.6% at `sigma = 8`. I
  reproduced the discriminator myself (V5). Also: the brief's "use nominal `V`
  instead of the AR(1) variance `A`" is not an alternative, it is an error —
  `Cov` from the explicit linear map matches `rho^|k-l| A[min]` to 1.07e-15 and
  `V` to 4e-02.
* **"The `sigma -> 0` extrapolation is an unjustified model."** Dead twice over:
  the order is derived with exact coefficients (P4), and the regulator can simply
  be switched off on a fixed lattice, so no model is needed at all (V2).
* **"The comparison is circular."** Dead in its strong form. Six distinct
  structural errors injected into the *estimator* are detected by factors 0.47 to
  2.0; the MC arm carries the same teeth; and the seven-source-distance sweep
  (fold amplitude spanning 196x, integrand median 784 -> 1750 Mpc, same table) is
  flat to 0.005 peak-to-peak where a wrong response exponent swings 2.011 ->
  1.263. That is a lambda-resolved constraint no table error can fake, and it is
  the single most persuasive piece of evidence in the whole exercise. It is
  currently not cited anywhere in `NOTES.md` or `PAPER_PROPOSAL.md`; it should be.
* **"Seed statistics are unsound because the estimator is heavy-tailed."** The
  mechanism is dead (see D).

---

## 3. What the paper may honestly claim

### The defensible statement

> With the input statistics held fixed, a direct Monte-Carlo of the stochastic
> Sachs equation reproduces the FK channel. In the limit of vanishing
> colored-noise correlation length the estimator's expectation agrees with the
> folded prediction to better than 0.3% at every separation tested from 1' to
> 17', and the Monte-Carlo reproduces that expectation to about 1%.

That is stronger and safer than what is currently proposed, because it separates
the two things that are known to different precision.

### Concrete numbers that are defensible

| quantity | value | status |
|---|---|---|
| exact expectation / fold, joint (`sigma -> 0`, `dlam -> 0`) | 0.9974 (1.015'), 0.9980 (5.30'), 0.9986 (17.28') | verified here, V3 |
| MC / exact expectation | 1.00 +/- 0.01 | verified here, V6, error inflated for block chi2 |
| MC / fold at `gamma = 1'`, extrapolated | **1.00 +/- 0.011** | V6, four independent blocks |
| placement share captured by the old estimator | 0.211 at 1' | V1 |
| nominal-calibration inflation at `sigma = 8` | 3.92 | NOTES, reproduced by three agents |

### Wording that must change

1. **Delete `0.997 +/- 0.006` (proposal line 156).** It is not reproducible from
   any seed block on disk.
2. **Replace `1.002 +/- 0.006` with `1.00 +/- 0.01`** at `gamma = 1'`, or drop
   the Monte-Carlo central value from the sentence entirely and quote the
   converged *expectation* instead, which is the tighter and better-understood
   number.
3. **Delete "better than a percent"** from the proposed `insights.tex`
   paragraph. It is true of the converged expectation; it is not true of the
   Monte-Carlo at the statistics actually collected.
4. **Delete "to within a few percent out to 17'"** as a claim about the
   Monte-Carlo. The 17.3' marker is `0.970 +/- 0.069`; the correct statement is
   "consistent with the folded prediction within the Monte-Carlo error, which
   grows to about 7% at 17'".
5. **Qualify the Mathematica bullet.** Not "the estimator captures the whole
   vertex, by a Mathematica identity", but "the deformation's three vertex legs
   are all reached by the estimator, by an identity proved in generic symbols;
   the identity fixes the placement bookkeeping, not the kernels".
6. **Do not present the gates as four independent confirmations.** Gate 1 is a
   structural identity about the estimator; gate 2's MC/exact is a sampling
   check on a shared model. The independent content is gate 4 against the fold,
   plus the external Order-0 anchor and the source-distance sweep.
7. **The marker table's inverse-variance mean of 1.012** is one seed block. If
   the markers go in the figure, either pool blocks or state that the error bars
   are single-block and correlated across gamma.

### Caveats that must travel with it

1. **Shared inputs — four, not one.** The tabulated `zeta_abc`, the callable that
   reads it, the background `D(lambda)`, and the `F` tensor are common to both
   sides. Name them. Quantified blind spots worth one clause: a radial-measure
   error of the `(1+z)^-4` kind would be a factor 11.5 and invisible; a 1% tilt
   in `D` biases FK by 0.8% and is invisible.
2. **The comparison is against the *symmetrised* vertex.** The deformation
   injects the 6-fold symmetrised cumulant exactly
   (`ref_inj/ref_sym = 1.000000000`), so the check is structurally blind to the
   table's leg asymmetry, which is 0.9% at 1' and 7.6% at 17.3'. That the fold
   behaves like the symmetrised object is inferred from a fingerprint, not read
   from source.
3. **Regulator.** State that the number is a `sigma_lambda -> 0` limit, and that
   the deficit at simulated values is a few percent. Already in the proposal;
   keep it.
4. **Perturbative in `zeta`, first order in `F`.** Already stated; keep it.
5. **Only the kappa-kappa entry.** Nothing here supports `xi_+`, `xi_-` or
   `kappa-gamma`, and the `F` contraction never touches `F_101`, `F_202`.

### What must NOT be claimed

Everything already on the proposal's "must not" list, plus: that the agreement
is better than half a percent *for the Monte-Carlo*; that four independent gates
confirm the result; that the identity proves the estimator computes the FK
diagram.

---

## 4. What must be fixed or re-run before the claim is made

**Blocking.**

1. **Resolve the three central values** (`1.002`, `0.997`, `1.016`). Decide which
   block or pool is quoted and make the LaTeX, the marker table and the summary
   agree.
2. **Requote the error bar** with the block-heterogeneity scale factor. Four
   independent blocks give chi2 = 9.22/3 (V6); six give 11.57/5. Use ~1.5-1.75x,
   i.e. `+/- 0.011` at `gamma = 1'`.
3. **Requote the central value at converged `dlam`, or state `n_lambda = 1000`
   explicitly.** The `1.002` is a property of that lattice. The converged number
   is 0.9974-0.9986 (V3).
4. **Fix the gate-1 wording** in `PAPER_PROPOSAL.md` section 0 (see 3.5 above).
   V4 shows the current phrasing claims something demonstrably false.

**Cheap and worth doing.**

5. **Drop the extrapolation.** The regulator can be switched off exactly on the
   same lattice (V2: `sigma = 0.25` gives `rho = 4.8e-04`). Doing so removes the
   `+0.089%` model bias for free. If the extrapolation is kept, match
   `sigma/dlam >= 8` at both anchor points (`N = 2000` at `sigma = 8`,
   `N = 4000` at `sigma = 4`), which drops the bias to `+0.026%`. **Caution:** do
   not read the direct switch-off blindly either — below `sigma ~ 2` the
   injection fidelity degrades from 8.6e-12 to 8.9e-05 and `ref_inj/fold` steps
   by 0.136% (V2). Quote the switch-off against `ref_sym`, or check the fidelity
   at whatever `sigma` is used.
6. **Add the source-distance sweep to the evidence.** It is already run
   (`review/circularity/_c4.log`, `_c5.log`), it costs minutes, and it is the
   only part of the check that is not degenerate in `lambda`. Its absence is the
   biggest gap between what was established and what is written down.
7. **Add the F-free `<kappa1^3>` criterion as the justification of the
   calibration**, with the alternative-calibration rejections. `NOTES.md`
   currently gives only nominal vs smeared. What makes the calibration *forced*
   is that this FK-independent criterion rejects every other smearing recipe and
   predicts each one's FK error. That is a much stronger argument than "the
   smeared one converges", and it is the answer to any referee who suspects
   tuning.
8. **State the second-order-accuracy requirement on the smearing quadrature**
   (residual F above), or note that `dlam` was refined at fixed `sigma` to check
   it.

**Corrections to `NOTES.md` itself.**

9. Section 3.2's PSD/graininess claim: "59 of 1000 nodes ... carry 6.9% of the FK
   kernel weight" overstates the consequence. The actual MC-vs-exact model
   mismatch from `_cholesky_psd` flooring is `+0.014%` at `sigma = 8` and
   `+0.034%` at `sigma = 4`. Fix the number.
10. Gate 3's `MC/exact = 1.045` at `sigma = 2` is attributed to sampling. It is
    `+2.87%` deterministic, from that same flooring. Relabel it. (The conclusion
    that `sigma = 4-8` is the usable window is unaffected and in fact
    strengthened.)
11. Gate 3's table row "exact / analytic 0.9734" is `exact/ref_inj`; against the
    fold it is 0.97308. Minor, but the two references are used interchangeably in
    places.
12. The `1.0084` "discretisation" in gate 2: only `1.0019` of it is
    discretisation. The rest is the vertex table's leg asymmetry (V3:
    `ref_raw/ref_sym = 1.0088`). NOTES' own correction list already says this
    (item 2); the gate-2 table should say it too.
13. Record the unbudgeted `O(1/batch_size)` bias of `-0.11% +/- 0.07%`, or remove
    it by using `ddof`-corrected covariances and an out-of-batch `c_node`.

**Not blocking, but flag for follow-up.**

14. The `Sigma2` zero-crossing near `lambda = 1309 Mpc` is a resolved feature of
    the tabulated object (visible on a 0.5 Mpc scan), not a one-node glitch, and
    it means the tabulated `Sigma2` is not a valid covariance density there. That
    is an upstream `corr_op` defect, independent of this exercise, and it is the
    reason the nominal route is singular rather than merely biased.
15. Nobody explained the seed-block heterogeneity. It is now detected twice
    independently. Adding seeds will not help until it is understood, and may
    make the quoted precision worse rather than better.

---

## 5. Bottom line

The work is sound and the exercise was unusually well conducted: the placement
identity is a genuine theorem, the regulator limit is genuinely derived, the
calibration is genuinely forced by an FK-independent criterion, and the
adversarial phase found real defects rather than rehearsing the defence. Two
independently written clean-room re-implementations landing on the same numbers
is strong evidence that the implementation is what it is claimed to be.

What the exercise does *not* support is the precision currently written down.
The `1.002 +/- 0.006` is a single-lattice, single-pooling number carrying a
`+0.4%` known-sign systematic and an error bar too small by a factor of about
1.7. Replace it with `1.00 +/- 0.01` for the Monte-Carlo, and quote the
converged expectation, `0.997-0.999` across `1'-17'`, as the sharp result. Fix
the gate-1 scoping sentence, name the four shared inputs rather than one, and
publish the source-distance sweep, which is the best evidence in the folder and
is currently invisible.
