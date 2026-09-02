# Task A: the blind audit, and how it compares with the prior findings

An independent reviewer was given the production path of the FK figures
and nothing else. It was not told what the previous session had concluded,
and the words grid, interpolation, corner and artifact were kept out of its
brief. It was told explicitly that a clean audit would be a useful result.
It was blocked from reading `notes/`, `SUMMARY.md` and `HANDOVER.md`.

It worked forward from `figure_manifest.md` through the generator scripts,
the `xi_*.npz` files, the sft-wick run config, the callable and the table,
and it tested rather than assumed at every step.

## Verdict

It reached the same conclusion. Quoting its bottom line: the plotted FK
curve is not a faithful evaluation of the quantity the paper defines.

It also cleared a control the previous session had not run: it re-ran the
production FK sweep from a redirected copy of the config and reproduced the
deployed NPZ bit for bit, `max|rel diff| = 0`. So every comparison it then
made isolates exactly one change.

## Agreements

| claim | prior note | blind audit |
|---|---|---|
| the fold queries only the collapsed family (1, cos g, cos g) | instrumented probe | instrumented probe, 1990656 logged evaluations, one sorted triple per separation, three leg times equal to 2.3e-13 |
| the lookup is frozen at the (1,1,1) corner over the plotted range | deployed column constant 0.5' to 600' | same, with the neighbour weights tabulated |
| the small-angle amplitude survives | 0.15% at 0.5' | 0.1% at 0.5', and correct to 6% out to 3' |
| the corrected FK/Order-0 ratio is flat at a few tenths of a percent to 1% | 0.4 to 1% | 0.5 to 1% |
| no FK / Order-0 crossover once the vertex is resolved | none in the physical range | never, over the whole plotted range |
| the sign reverses at wide angles | beyond about 183' | beyond about 230' |
| ell_high_max = 1000 is far from converged | tree x13.8 at the source shell for 16000 over 1000 | tree x3.1 at 2000, x10.4 at 8000 |
| the lambda_min floor is negligible at z_s = 5 | +0.12% median on the folded curve | vertex at the z = 0.1 shell is 4.6e-5 of the source shell |
| n_ell = 96, n_phi = 64 is adequate at fixed window | (not tested this way) | doubling both moves zeta by at most 0.07% |

The two independent estimates of the overstatement factor agree
numerically where they overlap. The prior note's dense/deployed ratio at
gamma = 5.30' is 0.8692; the audit's deployed/resolved at the same
separation is 1.150, and 1/1.150 = 0.870.

## What the audit adds

Five things the prior notes did not have.

1. **The plotted curve IS the interpolation weight.** Normalising the
   deployed FK curve by its value at 0.5' reproduces the inverse-distance
   weight on the (1,1,1) node to better than 1 part in 10^4, from 0.5' to
   1535'. This is stronger than "the deployed column is essentially
   constant", and it is internal to the deployed pipeline: it needs no
   second table to establish, which is exactly what a referee would ask
   for.

2. **It explains a specific caption in the draft.** The Fig. 15 caption in
   `sections/appendix.tex` describes "the sign change at gamma = 32 deg".
   That is where the four-neighbour set drops the (1,1,1) node, at
   cos = 6/7, gamma = 31.0 deg. The caption reads an interpolation switch
   as physics.

3. **A permutation-symmetry defect in the callable, separate from the
   grid.** The query triple is sorted before lookup and the resulting
   (3,3,3) tensor is filled with one channel value in every index
   permutation. That is right for zeta_TTT, which is symmetric, and wrong
   for the spin channels, which assign specific spins to specific legs. In
   the deployed table, rows sharing a sorted triple differ by a median of
   120 to 154% of their magnitude for the spin channels against 6.6% for
   zeta_TTT. It does not touch the deployed numbers, because the corner
   (1,1,1) has no permutation duplicates, but it does limit how precisely
   the corrected xi_- and xi_kappa_gamma can be quoted. The prior note had
   noticed the symptom (permutation duplicates differing by up to a factor
   15 at wide angles) without identifying it as a callable defect.

4. **The 16 lambda shells carry their own interpolation bias.** Rebuilt on
   a 31-node lambda grid, the deployed grid's linear interpolant
   overestimates zeta by 10 to 23% at mid lambda and 1 to 4% near the
   source. This is a few percent on the folded amplitude, so it is not a
   driver, but it is a distinct effect from the lambda_min floor, which is
   all the prior notes had examined.

5. **The corrected spin-2 panels are not at the numerical floor.** With a
   resolved vertex, xi_- reaches about +2.9e-7 at 28' and
   xi_kappa_gamma about -8.9e-7 at 17', i.e. 0.5 to 0.8% of Order-0,
   where the deployed run puts them at 1e-16 to 1e-9. The exact
   cancellation is a gamma to 0 statement, not a statement about the
   plotted range.

The audit also observed the internal contradiction plainly: Fig. 7,
`figures/zeta_driving_field_slices_draft.pdf`, is built from
`ell_band_decomp_results.npz`, which IS gamma-resolved, and it shows the
vertex decaying and changing sign inside the range over which the FK
figures show it flat. The draft already contains its own counterexample.

## Disagreements

None of substance. Two differences of emphasis are worth recording.

* The audit reports the sign reversal setting in beyond about 230' where
  the prior note put it beyond about 183'. Both are inside the region where
  the vertex passes through zero, and both used a different refinement of
  the collapsed grid, so the disagreement is the expected resolution of a
  zero crossing rather than a conflict.
* The audit calls `n_ell = 96` adequate; the prior addendum reports the
  deployed single-window build differing from an octave-banded build by
  10% at gamma = 1535'. These are compatible: the audit doubled the nodes
  inside one fixed window, the addendum changed the window partition. Node
  count is adequate; the linear node mapping over a wide window is not.

## Consequence

The gate in the handover is met: the independent audit confirms the
substance, so Task B proceeds.
