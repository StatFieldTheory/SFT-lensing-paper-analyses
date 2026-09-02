# The lambda_min floor of the line-of-sight integral

> **Correction, 2026-09-02.** The absolute FK amplitudes in this note are
> `ell_max = 1000` values obtained with the sorting callable and the production
> linear-in-lambda fold. The manuscript's FK is evaluated at `ell_max = 15360`
> with the permutation-aware callable, where the same quantity is `+1.9480e-5`
> at `gamma = 0.5'` (2.31% of Order-0), a factor 4.5 larger. Ratios between
> variants at a fixed cutoff are unaffected. See
> [`../FK_BASELINE_NUMBERS.md`](../FK_BASELINE_NUMBERS.md) for the full key and
> `sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/` for the
> current product.

Raised as a question on 2026-08-25: the LOS integral has a lambda_min
cutoff, which might not matter at high source redshift but could at low
redshift. It does, and the redshift dependence is the right instinct.

## What is cut

The kappa3 table's shells start at lambda = 396.63 Mpc and the callable
returns exactly zero outside them (scalar path lines 155-156 and 170-175,
batch path 229-230), so anything nearer than that contributes nothing.

The production fold uses a 24-node Gauss-Legendre rule on
[0, t_final = 2313.03], whose nodes span [5.6, 2307.5] Mpc. Six of the 24
fall below the floor: lambda = 5.6, 29.2, 71.4, 131.4, 208.2, 300.6 Mpc.
By plain measure that is 15% of the integration range.

lambda_min = 396.63 Mpc always corresponds to **z = 0.100**, independently
of the source redshift, because lambda is measured from the observer. So
what is excised is always the z < 0.1 foreground.

## How much it matters, by source redshift

Path fraction is the wrong measure; the lensing efficiency
K = a^2 chi (chi_s - chi) / chi_s vanishes toward the observer. The FK
vertex carries four response lines, so K^2 and K^4 bracket the weighting:

| z_s | path fraction | K^2 weight | K^4 weight |
|---|---|---|---|
| 5.0 | 0.055 | 1.46% | 0.34% |
| 3.2 | 0.065 | 1.76% | 0.44% |
| 2.5 | 0.073 | 2.03% | 0.54% |
| 1.7 | 0.090 | 2.68% | 0.81% |
| 1.0 | 0.128 | 4.64% | 1.79% |
| 0.5 | 0.223 | 12.74% | 7.77% |

At the paper's z_s = 5 the excised foreground is worth 0.3 to 1.5% of the
kernel weight, and the true figure is smaller still: the vertex itself is
about five orders of magnitude smaller at the innermost covered shell
(396.6 Mpc, zeta_TTT ~ -6.3e-19) than at the outermost (2326.6 Mpc,
~ -2.7e-14), so the near-observer region is doubly suppressed.

At z_s = 0.5, the lowest slice of the multi-redshift figure, the same cut
is worth 8 to 13% of the kernel weight. That is not negligible, and it is
systematic in one direction: it removes signal, and it removes MORE of it
the lower the source redshift.

## Why this matters for the redshift-trend statement

The draft states that the FK leakage "strengthens steadily with source
redshift". Two separate effects push in exactly that direction as an
artefact:

1. this lambda_min cut, which suppresses low-z_s slices by up to ~10% and
   high-z_s slices by well under 1%
2. the cosine-grid freezing documented in
   `finding_grid_artifact_measured.md`, whose severity also depends on the
   angular range each slice covers

Neither is evidence that the trend is wrong; both mean the deployed
figure does not measure it cleanly. The clean measurement exists:
`theory_redshift_dependence.md` section 3 folds the corrected tables at
each lambda_f(z_s) and finds the tree FK/O0 fraction flat in z_s to
-16%/+5%, with this cut worth -6.6% at z_s = 0.5 against -0.1% at
z_s = 5 (measured, section 4 item 3 there).

## A structural floor, separate from the inherited one

Building a table with shells at lambda = 6 Mpc fails: the Limber branch
needs P(k) at k = 2 ell_max / chi_h, which at chi_h = 4 Mpc/h is
500 h/Mpc, beyond the table's 206. The hard floor from that constraint is
chi = 2 ell_max / (k_max h) = 14.5 Mpc at ell_max = 1000. The deployed
floor of 396.6 Mpc is far above it and is simply inherited from the shared
L2 lambda grid, so it can be lowered without touching the P(k) table.
Raising ell_max, however, raises the structural floor in proportion.

## Empirical check: measured, and it confirms the estimate

A tree table was rebuilt on the same dense cosine grid with four extra
shells at lambda = 50, 90, 160 and 260 Mpc, covering four of the six
dropped nodes. The innermost two nodes (5.6 and 29.2 Mpc) cannot be
covered at all, see the structural floor below, and carry K^4 weight of
order 1e-11, so they are irrelevant.

Control first: on the sixteen shared shells the two tables agree to
machine zero (max |ratio - 1| = 0.00e+00), so the extra shells perturbed
nothing. The vertex on the new shells is well behaved and does not grow
toward the observer: zeta_TTT at the smallest separation runs
-2.07e-19 (lambda = 50 Mpc), -2.87e-19 (90), -3.82e-19 (160),
-4.87e-19 (260), against -6.26e-19 at the old innermost shell and
-2.70e-14 at the outermost.

Folding both through the production pipeline, xi_FK for kk (absolute
values here carry the production lambda-interpolation inflation, +11% at
this gamma, see `theory_redshift_dependence.md` section 2.2; the floor
CHANGES in the last column are unaffected, both columns sharing the
interpolation):

| gamma ['] | floor 397 Mpc | floor 50 Mpc | change |
|---|---|---|---|
| 0.50 | +4.28149e-06 | +4.28624e-06 | +0.111% |
| 5.30 | +3.72698e-06 | +3.73147e-06 | +0.120% |
| 17.28 | +1.22382e-06 | +1.22689e-06 | +0.251% |
| 56.27 | +2.22997e-07 | +2.24376e-07 | +0.619% |
| 183.26 | +4.19107e-09 | +4.54075e-09 | +8.344% |
| 596.89 | -1.86608e-09 | -1.80757e-09 | -3.135% |

Median change over gamma < 700 arcmin: +0.12%. At the separations that
carry the signal the cut is worth about a tenth of a percent, consistent
with the K^4 estimate of 0.34% (which is an upper bound, since it ignores
the vertex being five orders smaller near the observer).

The change grows to a few percent only beyond about 150 arcmin, where the
FK curve is passing through its sign change and the absolute values are
1e-9 or below, i.e. where a percentage is not meaningful.

**Conclusion for z_s = 5: the lambda_min floor is a 0.1% effect and does
not affect any statement in this study.** The redshift dependence of the
kernel weighting is a statement about the multi-redshift figure, not
about the main analysis, and it is measured across source redshifts in
`theory_redshift_dependence.md` (-6.6% at z_s = 0.5 down to -0.1% at
z_s = 5 at 0.5 arcmin).
