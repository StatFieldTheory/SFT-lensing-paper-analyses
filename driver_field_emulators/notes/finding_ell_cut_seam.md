# The ell_cut = 60 seam: keeping LOW at tree level is safe

Measured 2026-08-25 with `code/probe_ell_cut_seam.py`, plus the
absolute-share follow-up that the first reading required.

## Why the first reading was misleading

The probe compares, on the shared band (40, 60], a LOW band-difference
(exact Wigner-3j at ell_max = 60 minus ell_max = 40) against a HIGH window
integral over the same band, with the dispatcher bypassed on both sides.
Read with the usual relative gate it returns BAD in every cell: the two
methods disagree by about 100%.

That reading is wrong, because the relative gate is the wrong criterion
here. The band values are 1e-20 to 1e-27 while the full vertex is 1e-14 to
1e-19, so the gate was comparing two quantities that barely contribute. The
band (40, 60] carries between 1e-6 and 2.5e-2 of the total vertex,
depending on the cell. A 100% disagreement inside a band worth 0.001% of
the answer is not a defect worth acting on.

The boundary-validation rule says to evaluate both methods directly at the
threshold, which the probe does correctly. What it does not say, and what
this case adds, is that the verdict must be normalised to the quantity the
result feeds, not to the band itself.

## The question that actually mattered

Once HIGH takes a nonlinear bispectrum while LOW stays tree-level, the two
branches model different physics. The cost of that, expressed as a
fraction of the TOTAL vertex:

| triple | lambda [Mpc] | \\|B_nl - B_tree\\| on the band / \\|total\\| |
|---|---|---|
| gamma = 0.50' | 1200 | 1.2e-04 |
| gamma = 42.2' | 396.6 | 9.6e-03 |
| gamma = 1000' | 396.6 | 1.5e-02 |
| gamma = 1000' | 1200 | 6.1e-03 |
| open triangle | 396.6 | 9.4e-04 |

Worst cell: 1.5%. So keeping the exact-3j branch at tree level, which its
FFTlog machinery requires anyway, costs at most about 1.5% of the vertex,
and far less at the small angles and far shells that dominate the FK fold.
That is an order of magnitude below the BiHalofit model error and two
orders below the cutoff sensitivity. **The LOW-stays-tree decision is
validated.**

## A structural fact worth recording

The share table also exposes something not previously noted: at the
farthest shell the two branches largely CANCEL. At lambda = 2326.6 Mpc the
LOW branch is 2.90 times the total for gamma = 42', and 1.02 times it for
the widest angles, so HIGH carries the opposite sign at nearly the same
magnitude. The total there is a difference of two larger numbers, which
makes it more sensitive to errors in either branch than the raw magnitudes
suggest. At the near shells LOW is instead negligible (1e-5 to 1e-3 of the
total) and HIGH carries everything.

This is a property of the existing design, not of the nonlinear swap, but
it is worth knowing when interpreting far-shell numbers.

## Dispatch error, for completeness

Expressed against the total, the exact-3j versus flat-sky-Limber
disagreement on the band ranges from 1.4e-5 (far shell, small angle) to
1.7e-1 (near shells, wide angles, where the total is five orders smaller
than at the far shell). It is a pre-existing property of the LOW/HIGH
split, unchanged by anything in this work, and it is largest exactly where
the vertex contributes least.
