# SUPERSEDED: the HIGH-branch cutoff behaviour

> **Superseded on 2026-08-25 by `finding_ell_max_resolved.md`.** The
> octave-band measurement shows the vertex integral CONVERGES for both the
> tree and the nonlinear bispectrum; the cutoff is a numerical convergence
> parameter, and the deployed value of 1000 sits below the integrand's
> peak. The reasoning below correctly identified the mechanism, an
> un-damped coincident hard pair, but drew the wrong conclusion from a
> partial sum measured over a single decade. Kept for the record.

# Original note: the HIGH-branch cutoff is a physical UV scale

Measured 2026-08-25 with `code/sweep_ell_high_max.py` on the collapsed
family (1, cos gamma, cos gamma) that the FK two-point fold actually
queries. Affects the DEPLOYED table and hence the published FK amplitude,
not only the nonlinear rebuild.

## The control: the deployed quadrature is fine

Sweeping `ell_high_max` at fixed `n_ell = 96` conflates two effects, since
the same 96 Gauss-Legendre nodes then cover a wider interval. Separating
them:

* At the deployed cutoff of 1000, the 96-node rule agrees with a 384-node
  rule to 0.1-1.1% across gammas and shells, for both tree and BiHalofit.
  The production quadrature is converged.
* An earlier naive sweep showed a 3.5x discrepancy between 96 and 192 nodes
  at cutoff 4000. That was the node-dilution artifact, not a defect of the
  deployed rule.

## The finding: the cutoff itself sets the vertex value

At CONSTANT node density (96 nodes per 1000 multipoles), zeta_TTT(cutoff)
relative to zeta_TTT(4000):

| model | gamma | lambda [Mpc] | cut 1000 | cut 2000 | cut 4000 |
|---|---|---|---|---|---|
| tree | 0.5' | 396.6 | 0.744 | 0.898 | 1.000 |
| tree | 0.5' | 1200 | 0.487 | 0.760 | 1.000 |
| tree | 0.5' | 2326.6 | 0.147 | 0.464 | 1.000 |
| tree | 42' | 2326.6 | 0.118 | 0.489 | 1.000 |
| BiHalofit | 0.5' | 1200 | 0.073 | 0.287 | 1.000 |
| BiHalofit | 42' | 2326.6 | 0.116 | 0.483 | 1.000 |
| BiHalofit | 117' | 1200 | 0.384 | 0.548 | 1.000 |

Reading the doubling behaviour: the tree vertex grows with a decelerating
exponent (d ln zeta / d ln ell_max falls from about 1.7 to 1.1 at the worst
cell, and from 0.6 to 0.4 at a typical one), so it is slowly convergent.
The BiHalofit vertex grows with an almost constant exponent near 1.8-2.0,
i.e. it shows no sign of converging over the tested decade.

## Why this is expected physics, not a bug

The Ricci focusing field is essentially the density contrast: under Poisson
plus the transverse-Laplacian screen projection, Phi_00 is proportional to
A(a) delta. The collapsed family puts TWO legs on the same direction, so
zeta(1, c, c) contains a coincident-pair moment, the three-point analogue
of sigma^2(0). For a CDM spectrum that object is UV-divergent: logarithmic
with linear P(k) at high k, power-law once the nonlinear P is used. The
alpha kernel supplies no angular damping for the coincident pair
(P_l(1) = 1), so nothing in the integrand cuts the hard-pair sum off.

`ell_high_max` therefore acts as a smoothing scale on the driving-field
three-point cumulant, and the deployed value of 1000 is a choice, not a
converged limit. A nonlinear bispectrum makes the sensitivity stronger
because it adds small-scale power exactly where the divergence lives.

## What this does and does not imply

It does NOT immediately imply the published FK number is wrong. The vertex
is folded with response propagators over lambda before it becomes an
observable, and the fold could suppress or cancel the UV-sensitive part.
What it does imply:

1. The FK amplitude in the deployed pipeline carries an undocumented
   dependence on ell_high_max, which no existing test bounds.
2. Comparing a nonlinear rebuild against the tree baseline at a FIXED
   cutoff is still meaningful (both are smoothed the same way), but the
   absolute FK amplitude of either is cutoff-conditioned.
3. The Monte-Carlo validation of the workflow (appendix figure
   `fig: val mc`) does not test this: it is seeded with the SAME input
   statistics, so it inherits the same cutoff.

## The decisive test (not yet run)

Build two otherwise identical tables at ell_high_max = 1000 and 2000, run
the production FK fold on both, and compare xi_FK(gamma) and the FK kk
value at 0.5'. Three outcomes:

* xi_FK stable to a few percent: the fold suppresses the UV-sensitive part,
  the cutoff is harmless, and the result should be stated with that
  evidence.
* xi_FK moves roughly like the vertex: the FK amplitude is UV-conditioned
  and the paper must state ell_high_max as a physical smoothing scale, with
  a justified choice (for example matched to the angular resolution of the
  observable) rather than an inherited default.
* Intermediate: quantify and report the residual as a systematic.

Until that test is run, quote FK amplitudes as conditional on
ell_high_max = 1000, and keep the cutoff fixed when comparing models so the
comparison is not contaminated.
