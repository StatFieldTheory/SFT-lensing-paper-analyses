# Does a correction at this level matter for cosmology?

Assessment, 2026-08-25. The numbers come from the measurements in this
folder; the placement against survey requirements is judgement, and is
labelled as such.

## Where it sits

After the grid correction, FK / Order-0 is flat at one to two percent
across three decades in angle (2.06% at 0.5 arcmin for tree level at
ell_max = 8192, less at lower cutoffs, more with a nonlinear bispectrum,
and not yet converged). For orientation against the systematics a
Stage-IV cosmic-shear analysis already handles:

| effect | typical size | treatment |
|---|---|---|
| baryonic feedback | 1 to 20% for k > 1 h/Mpc | modelled and marginalised |
| intrinsic alignments | 1 to 10% | modelled |
| **FK, corrected** | **1 to 2%, not converged** | not currently considered |
| reduced-shear correction | about 1% | included by Euclid |
| post-Born / lens-lens | below 1% for galaxy lensing | usually neglected |

So by size it lands where a survey would at least have to check it, and
plausibly include it, but well below the dominant systematics.

## The size is not the deciding factor

The more useful observation is that the corrected FK / Order-0 is nearly
CONSTANT with angle. A flat multiplicative offset is close to degenerate
with the amplitude parameters: a 1% shift in xi corresponds to about 0.5%
in sigma_8. Stage-IV precision on sigma_8 is of order 0.3 to 0.5%, so an
unmodelled flat 1% would bias the amplitude by roughly one sigma, and
would show up as a shifted sigma_8 rather than as a distorted shape.

That reframes the question. What decides whether this term matters is
whether it has a distinctive REDSHIFT or SCALE dependence, not whether it
is 1% or 3%:

* a flat offset is largely absorbed by sigma_8 or S_8
* a term that grows with source redshift cannot be absorbed by a single
  amplitude in a tomographic analysis, and becomes a systematic that must
  be modelled

The draft claims exactly the second behaviour, that the leakage
strengthens steadily with source redshift. That claim is the ONE thing
this session did not test, and it is contaminated by two effects that both
push in the same direction: the cosine-grid freezing, and the lambda_min
floor, which removes 8 to 13% of the kernel weight at z_s = 0.5 against
under 1.5% at z_s = 5 (`finding_lambda_min_cut.md`).

**So the cosmological priority is the opposite of the intuitive one: the
redshift trend matters more than the amplitude.**

## A five-fold discrepancy that needs resolving

An earlier internal check (2026-06-11) reconciled the FK result with the
CMB-lensing literature by noting that at ell = 100 to 1000 the ratio
(FF + FK) / Order-0 was 0.14 to 0.24%, consistent with the post-Born
corrections of Pratten and Lewis at 0.2% or below. On the dense grid the
FK ratio over the same multipole range is about 1%
(`finding_cl_and_credibility.md`), roughly five times larger.

That reconciliation therefore has to be redone. Two readings are open:
the corrected FK is genuinely larger than standard post-Born estimates,
which would be a real result since FK is not the same object (it is the
coupling of nonlinear propagation to driving-field non-Gaussianity, not a
perturbed-path correction); or the two sides are not comparing the same
quantity. Nothing here decides between them.

## Summary

At present numbers this is a correction at the threshold of Stage-IV
relevance which would be largely absorbed by sigma_8 if it is as flat as
measured. The redshift dependence that decides the question is now
measured at ell <= 1000 (`theory_redshift_dependence.md` section 3): at
tree level the FK/O0 fraction is flat in source redshift to -16%/+5%,
so sigma_8 absorbs it there; with nonlinear input the fraction FALLS
with source redshift (~2.7% at z_s = 0.5 against ~0.8% at z_s = 5), a
z_s-DEPENDENT term after all, but one that matters most for the
low-redshift bins rather than growing toward high z_s as the draft's
figure suggested. What still gates a final verdict is the amplitude
settling, which requires the convergence work at a k-fixed cutoff. The
1% figure by itself does not settle it either way.
