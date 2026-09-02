# Proposed revision: the FK numerical check

*2026-08-28. Nothing applied. Three text blocks and one figure. Every number is
from `NOTES.md`; the section references there carry the evidence.*

---

## 0. What changed, and what may now be said

The FK channel previously had no stochastic check: `simulate_fk_vr` was shown
(2026-08-26) to measure only a placement share of the vertex, so its markers
were withdrawn and the appendix was rewritten to say the FK check is
deterministic. A placement-complete, variance-controlled estimator
(`driver_field_emulators/code/mc_fk_complete/`) now reproduces the analytic FK
from a direct numerical solution of the stochastic Sachs system.

**The claim, and why it is well posed.** The question is: *holding the
driving-field statistics fixed, does the stochastic solution agree with the
semi-analytic formalism?* Sharing the vertex table is therefore not a leak, it
is the controlled variable. And for FK the input matching is exact in a strong
sense (`NOTES.md` section 5d):

* the FK diagram carries four `R` propagators and **zero `C` propagators**, so
  the formalism's FK never reads the driving-field covariance;
* the Monte-Carlo's FK is invariant to `Sigma2` to about `1e-4` under an
  arbitrary lambda-distortion (ramps, tilts, oscillations, steps that swing
  Order-0 by a factor 2.5), under zeroing the entire cross-ray block, and under
  replacing `Sigma2` by an isotropic `(tr/6) I`.

So both sides depend on `(zeta, D)` alone and both are shared exactly. This
matters because the one place the two treatments genuinely differ *is* the
two-point sector: the formalism uses the full non-local `C(lam1, lam2)`, the
Monte-Carlo a lambda-local idealisation, and that costs **15% in Order-0 at
0.5'**, falling to 6% at 150'. It provably cannot touch FK.

**Numbers.**

| quantity | value |
|---|---|
| MC / fold at `gamma = 1'`, `sigma_lambda -> 0` | **1.01 +/- 0.01** |
| converged expectation / fold (`sigma -> 0`, `dlam -> 0`) | 0.9974 (1'), 0.9980 (5.3'), 0.9986 (17.3') |
| placement share of the withdrawn estimator | 0.211 at 1' |

**What must NOT be claimed.** That the vertex table itself is validated; that
the Monte-Carlo agreement is better than about a percent; that the Monte-Carlo
is non-perturbative in the three-point cumulant (it is exactly first order in
it, which is the order the diagram lives at); or anything about `xi_+`, `xi_-`
or `kappa-gamma` -- only the kappa-kappa entry was tested, and the `F`
contraction there never touches `F_101` or `F_202`.

**Do not present the checks as a chain of independent confirmations.** The
placement identity is a structural theorem about the estimator, and the
Monte-Carlo-versus-its-own-expectation closure is a sampling check on a shared
model: in a world with the response exponent set to 1, where FK is 18% wrong,
that closure still reads 0.987. The non-tautological content is the comparison
to the fold, plus the Order-0 lever arm and the source-distance sweep.

---

## 1. `sections/insights.tex`, the Validation paragraph (l. 313)

### Current

```latex
\paragraph*{Validation.} \edited{The diagrammatic prediction is checked in two
ways: with the input statistics held fixed, a Monte-Carlo integration of the
Sachs equation reproduces the analytic FF channel, and a direct contraction of
the tabulated three-point cumulant with the FK kernel reproduces the folded FK
channel to about a percent at the separations that carry the signal
(Fig.~\ref{fig: val mc}; Appendix~\ref{append: validation}); the first check
is independent of the diagram expansion, the second of the workflow's
implementation of it.}
```

### Proposed

```latex
\paragraph*{Validation.} \edited{With the input statistics held fixed, a
Monte-Carlo integration of the stochastic Sachs equation reproduces both
Order-2 channels (Fig.~\ref{fig: val mc}; Appendix~\ref{append: validation}).
The FK comparison is the sharper of the two: its diagram carries no
driving-field two-point propagator, so both the simulation and the
diagrammatic fold depend on the prescribed three-point cumulant alone, and they
agree to about a percent. A direct contraction of that cumulant with the FK
kernel, independent of the workflow code, reproduces the same channel.}
```

*Same length, and it says more.* The old "first check / second check" split no
longer maps onto FF / FK, since both now have both kinds of check; the new
sentence instead states what makes the FK comparison clean.

---

## 2. `sections/conclusion.tex`, the numerical-side sentence (l. 137)

### Current

```latex
\edited{On the numerical side, the FK leakage reported here has been validated
for a prescribed analytic three-point cumulant: an independent contraction of
that cumulant with the FK kernel reproduces the folded channel to about a
percent at the separations that carry the signal, and a direct Monte-Carlo of
the stochastic Sachs equation reproduces the companion FF channel
(Fig.~\ref{fig: val mc}).}
```

### Proposed

```latex
\edited{On the numerical side, the FK leakage reported here has been validated
for a prescribed analytic three-point cumulant, in two independent ways: a
direct Monte-Carlo of the stochastic Sachs equation reproduces it to about a
percent, and a contraction of the cumulant with the FK kernel, independent of
the workflow code, reproduces it as well (Fig.~\ref{fig: val mc}).}
```

*Two lines shorter.*

---

## 3. `sections/appendix.tex`, "Workflow validation" (ll. 570-612)

### 3a. Opening paragraph (l. 573), one clause changed

Replace

```latex
The three-point
channel FK is checked deterministically: contracting the tabulated cumulant
$\zeta_{abc}$ with the FK kernel by plain quadrature, independently of the
workflow code, reproduces the folded FK channel to about a percent at the
arcminute separations that carry the signal, and to within $7\%$ out to
$17'$.}
```

with

```latex
The three-point
channel FK is checked twice: the same Monte-Carlo reproduces the folded channel
to about a percent at $\gamma=1'$, and contracting the tabulated cumulant
$\zeta_{abc}$ with the FK kernel by plain quadrature, independently of the
workflow code, reproduces it to about a percent at the arcminute separations
that carry the signal and to within $7\%$ out to $17'$.}
```

### 3b. The paragraph at l. 588: three clauses became false

Three statements in this paragraph are no longer true. Everything else in it is
correct, reads well, and should stay untouched.

1. **"The FK channel is not given Monte-Carlo points"** -- it now is.
2. **"the recovered FK amplitude scales roughly linearly with that length"** --
   that was the behaviour of the withdrawn estimator, whose amplitude was
   proportional to `sigma_lambda` with an intercept consistent with zero. The
   corrected estimator's amplitude converges; it is the *deficit* that is linear
   in `sigma_lambda`, and it extrapolates to the folded value.
3. **"a stochastic estimate cannot fix the FK amplitude on its own. The
   deterministic contraction above can, and does."** -- it can, and does.

#### Current (13 lines)

```latex
\edited{The FK channel is not given Monte-Carlo points, for a reason worth
stating. The FF vertex is a term in the equation of motion: it can be switched
on and off at a fixed random draw, and differencing the two matched runs
cancels the common fluctuation, so the noise scales with the FF effect itself.
The FK signal is instead a property of the input statistics, the response of
$\langle\kappa\kappa\rangle$ to the three-point cumulant of the driving field; a
third moment measured from realisations scatters according to the Gaussian
fluctuations of the field, however small the cumulant. A simulable skewed
field must also carry a finite correlation length along the ray (the pointwise
skewness diverges as $\sigma_\lambda^{-1/2}$ in the white-noise limit), and
the recovered FK amplitude scales roughly linearly with that length, so a
stochastic estimate cannot fix the FK amplitude on its own. The deterministic
contraction above can, and does.}
```

#### Proposed (12 lines): two edits, nothing added

```latex
\edited{The two channels are estimated differently, for a reason worth
stating. The FF vertex is a term in the equation of motion: it can be switched
on and off at a fixed random draw, and differencing the two matched runs
cancels the common fluctuation, so the noise scales with the FF effect itself.
The FK signal is instead a property of the input statistics, the response of
$\langle\kappa\kappa\rangle$ to the three-point cumulant of the driving field; a
third moment measured from realisations scatters according to the Gaussian
fluctuations of the field, however small the cumulant. A simulable skewed
field must also carry a finite correlation length along the ray (the pointwise
skewness diverges as $\sigma_\lambda^{-1/2}$ in the white-noise limit), so the
FK estimate is made at several correlation lengths and extrapolated to zero,
the residual being linear in that length.}
```

Diff, in words: the opening clause changes from "The FK channel is not given
Monte-Carlo points" to "The two channels are estimated differently"; and the
final one-and-a-half sentences are replaced by the clause "so the FK estimate is
made at several correlation lengths and extrapolated to zero, the residual being
linear in that length". Every surviving word was already in the paragraph.

*An earlier draft of this proposal added a sentence here about the FK diagram
carrying no driving-field two-point propagator. That has been dropped: the
paragraph is about why the two channels are estimated differently, and a remark
about what the comparison is blind to answers a different question in language
the paragraph has not set up.*

---

## 3c. One imprecision this work exposed, flagged not proposed

The appendix and the figure caption both say the Monte-Carlo is "driven by the
same input statistics". At the level of the three-point cumulant that is exact.
At the level of the two-point function it is not: the workflow consumes the full
non-local `C(lam1, lam2)`, while the Monte-Carlo drives a lambda-local field of
density `Sigma2 = d/dlam[D^4 C(lam,lam)]/D^4`. Computed both ways from the same
table, the collapse costs 15% at `0.5'` and 6% at `150'`.

This does not touch FK, which carries no driving-field two-point propagator and
is measurably invariant to `Sigma2` (`NOTES.md` section 5d). It does bear on the
FF markers, which scale as `Sigma2^2`. Whether the phrase is worth qualifying is
a judgement call, and it is a pre-existing wording question rather than anything
this work changed; raising it here so it is on the record.

---

## 4. Figure `fig: val mc` (`figures/appendix_mc_workflow.pdf`)

**What the figure can and cannot carry.** Mocked up in
`mock_val_figure.pdf` (this folder; nothing written to `figures/`). The panel is
log-log over 6.14 decades in 4.26 inch, so 0.694 inch per decade. On that axis:

| quantity | height on the page |
|---|---|
| the FK agreement, 1% | **0.24 pt** |
| the FK error bar, 0.5% | **0.11 pt** |
| the FF offset at small gamma, ~20% | 3.96 pt |
| (marker diameter, for scale) | 6.5 pt |

So the FK markers sit on the line with no visible deviation and no visible error
bar: the figure can show that the check EXISTS and is consistent, but it cannot
show how good it is. The percent-level number has to live in the appendix text.
That is not a reason to leave the markers out -- "markers on the line" is the
normal idiom for a validation figure, and the FF markers already do exactly that
-- but it should be understood, and the caption should not imply the reader can
read the precision off the plot.

A ratio sub-panel WOULD make the 1% visible. It is not proposed: it would also
put the FF markers' ~15-20% scatter on display at full height, and it is a
structural change to a figure in a manuscript under review.

**Change.** Restore Monte-Carlo markers on the FK series, from
`driver_field_emulators/code/mc_fk_complete/` only -- never from
`sachs_mc_core.simulate_fk_vr`.

**Values.** Four independent 24-seed blocks (bases 6010001/6020002/6030003/
6040004), `n_lambda = 1000`, `n_real = 24000`, paired-seed linear
`sigma_lambda -> 0` extrapolation from `sigma_lambda = 8` and `4`. Error bars
from the **block-to-block scatter**, not the within-block seed error: independent
blocks scatter by 1.75x the seed SEM, for reasons not yet understood, so seed
errors would understate by that factor.

| gamma | FK (Monte-Carlo) | error | analytic fold | ratio |
|---|---|---|---|---|
| 1.0' | 1.6161e-05 | 8.5e-08 | 1.5977e-05 | 1.012 +/- 0.005 |
| 2.6' | 1.0746e-05 | 8.2e-08 | 1.0617e-05 | 1.012 +/- 0.008 |
| 6.7' | 5.5762e-06 | 8.0e-08 | 5.5344e-06 | 1.008 +/- 0.014 |
| 17.3' | 1.9649e-06 | 6.6e-08 | 1.9664e-06 | 0.999 +/- 0.033 |

Cached in `_markers_pooled.npz`.

**Grid.** Only the four published `GAMMA_MC` nodes inside the range where the
analytic FK is converged. No FK markers beyond `gamma ~ 1 degree`, where the
figure already draws the FK line faint.

**Consistency.** Both sides must read `products/table_permclosed_cut15360*`. The
deployed figure's FK line is already bound to the corrected fold by
`code/figures_corrected/make_val_figure.py --fk-xi`, and the Monte-Carlo binds
the same table through `perm_aware_kappa3_callable`; that pairing is verified
bit-identical.

### Caption

#### Current

```latex
  \caption{\emph{Workflow validation.} The convergence two-point function
  $\xi_\kappa(\gamma)$ split into its Order-$0$ (O0), nonlinear-propagation (FF)
  and three-point cumulant (FK) channels: the analytic workflow (lines) against a direct
  Monte-Carlo of the stochastic Sachs equation driven by the same input
  statistics \edited{(markers, FF channel). The FK channel is validated
  separately, by the deterministic contraction described in the text. The FK series is drawn faint beyond $\gamma\simeq1^\circ$, where its multipole sum no longer converges (Section~\ref{subsec: cutoff}) and no value is quoted.} Filled and
  hollow markers denote positive and negative values.}
```

#### Proposed

```latex
  \caption{\emph{Workflow validation.} The convergence two-point function
  $\xi_\kappa(\gamma)$ split into its Order-$0$ (O0), nonlinear-propagation (FF)
  and three-point cumulant (FK) channels: the analytic workflow (lines) against a
  direct Monte-Carlo of the stochastic Sachs equation driven by the same input
  statistics (markers). \edited{FK markers are extrapolated to vanishing
  colored-noise correlation length, with error bars from the scatter of
  independent realisation blocks. The FK series is drawn faint beyond
  $\gamma\simeq1^\circ$, where its multipole sum no longer converges
  (Section~\ref{subsec: cutoff}) and no value is quoted, and no markers are
  placed there.} Filled and hollow markers denote positive and negative values.}
```

Written timelessly: it does not say the markers are "now" restored.

---

## 5. Not proposed, and why

**Order-0 is not a comparison and should not be presented as one.** The
Monte-Carlo drives a lambda-local field of density
`Sigma2 = d/dlam[D^4 C(lam,lam)]/D^4`; the formalism uses the full non-local
`C(lam1, lam2)`. Computed both ways from the same table, the collapse costs 15%
at `0.5'`, 13% at `6.7'`, and 6% at `150'` -- gamma-dependent, so the single
`ANCHOR_C0 = 1/0.881` cannot absorb it. The figure already shows O0 as the
analytic line only; keep it that way.

**FF markers stay as they are.** FF scales as `Sigma2^2` and so inherits that
mismatch quadratically, and the measurement cannot currently resolve it: FF is a
drift term with `O(dlam)` error, and at `gamma = 1'` the ratio to the analytic
value reads 1.04 / 1.27 / 1.22 / 0.99 / 0.97 for `n_lambda` = 500 to 8000, each
with about 10% seed error. Making FF a few-percent test is a separate piece of
work.

**Nothing about `xi_+`, `xi_-` or `kappa-gamma`.** Untested; the `F` contraction
in the kappa-kappa entry never reaches `F_101` or `F_202`.

---

## 6. Also worth updating

`driver_field_emulators/CHANGES_OUTSIDE_THIS_FOLDER.md` item 7 records the
removal of the FK markers and its reasoning. That reasoning (the placement
share) was correct and should be preserved as history, with the resolution
appended rather than overwritten.
