# Three-effect decomposition of the FK two-point amplitude

> **Correction, 2026-09-02.** The absolute FK amplitudes in this note are
> `ell_max = 1000` values obtained with the sorting callable and the production
> linear-in-lambda fold. The manuscript's FK is evaluated at `ell_max = 15360`
> with the permutation-aware callable, where the same quantity is `+1.9480e-5`
> at `gamma = 0.5'` (2.31% of Order-0), a factor 4.5 larger. Ratios between
> variants at a fixed cutoff are unaffected. See
> [`../FK_BASELINE_NUMBERS.md`](../FK_BASELINE_NUMBERS.md) for the full key and
> `sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/` for the
> current product.

Measured 2026-08-25 at the observable level. Each pair of runs changes ONE
variable and holds everything else fixed: the same 80-row collapsed cosine
grid, the same 16 lambda shells, the same LOW cache, the same cosmology,
and the same production sft-wick geometry (n_gauss = 24, t_final =
2313.0289). Products and their generating scripts are in
`driver_field_emulators/{products,code}`.

| variant | table | xi_FK kk at 0.5' |
|---|---|---|
| deployed | covgrid16, tree, cutoff 1000 | +4.287723e-06 |
| grid fixed | dense (1,c,c), tree, cutoff 1000 | +4.281487e-06 |
| cutoff doubled | dense (1,c,c), tree, cutoff 2000, n_ell 192 | +8.432787e-06 |
| nonlinear | dense (1,c,c), BiHalofit, cutoff 1000 | +7.089177e-06 |

(Caveat added later the same day: all folded values in this note read
the shell grid through the production's linear-in-lambda interpolation,
which a truth probe subsequently showed inflates the z_s = 5 fold by
+11% at 0.5 arcmin and up to +30% at tens of arcmin; the grid-converged
tree baseline is +3.84e-6. The RATIOS between variants are unaffected,
since every variant shares the interpolation. See
`theory_redshift_dependence.md` section 2.2.)

So, relative to the grid-fixed tree baseline at 0.5 arcmin:

* cosine-grid fix: x0.9985 (negligible HERE; see below for wide angles)
* ell_high_max 1000 -> 2000 at constant node density: **x1.97**
* tree -> BiHalofit at fixed cutoff: **x1.66**

## The three effects, in the order they matter

**1. The cutoff dominates.** Doubling ell_high_max nearly doubles the FK
amplitude. This settles the question left open by
`finding_ell_max_uv_sensitivity.md`: the fold does NOT suppress the
vertex's UV sensitivity, it passes it through. The FK amplitude is
conditioned on an arbitrary numerical cutoff, and the dependence is of
order unity per doubling, with no sign of converging. It is not a constant
rescaling either: the cutoff ratio varies with separation (1.97 at 0.5',
1.14 at 8.5', 1.43 at 21.9', 1.59 at 145'), so it changes the shape as
well as the amplitude.

**2. The nonlinear bispectrum matters less than the cutoff.** The upgrade
this whole programme set out to make, tree to BiHalofit, raises FK by 1.66
at 0.5 arcmin, i.e. less than the effect of an unjustified numerical
parameter. That ordering is the headline: physics accuracy in B is not the
limiting factor today.

**3. The grid governs the angular shape.** At 0.5 arcmin the grid fix is a
0.15% effect, but at wide angles it is two to three orders of magnitude
(see `finding_grid_artifact_measured.md`). With the grid corrected, the
nonlinear FK tracks Order-0 at a roughly constant 0.8 to 1.5% across the
whole range, and no crossover occurs in the physical range.

| gamma ['] | O0 | \\|FK_nonlinear / O0\\| |
|---|---|---|
| 0.50 | +8.443e-04 | 0.0084 |
| 3.31 | +5.547e-04 | 0.0120 |
| 8.51 | +3.132e-04 | 0.0146 |
| 21.88 | +1.253e-04 | 0.0079 |
| 56.27 | +2.854e-05 | 0.0112 |
| 144.71 | +2.195e-06 | 0.0133 |
| 957.24 | -8.642e-08 | 0.0023 |

The formal crossing at 1535 arcmin sits between curves already at the
numerical floor and changing sign, so it carries no information.

## What this means for the programme

The nonlinear-bispectrum upgrade is worth having, but it cannot be
presented as an accuracy improvement while a factor-two numerical cutoff
sensitivity sits underneath it. The order of business is now:

1. Decide what ell_high_max physically means for this observable and fix it
   on a stated basis, or demonstrate a regime where the answer stops
   depending on it. Until then every FK amplitude is conditional.
2. Re-derive the FK numbers on the corrected grid, since the wide-angle
   behaviour changes qualitatively.
3. Only then compare tree against BiHalofit as a physics statement.

## Implications for the draft (recorded here, not acted on)

The draft has not been touched and will not be in this session. For the
record, the statements a corrected analysis would need to revisit are:

* the wide-angle FK behaviour and the reported crossover with Order-0 near
  2 degrees, which the grid analysis attributes to the frozen corner value
* any FK amplitude quoted without an accompanying ell_high_max statement
* the figures built on the FK curves (the Order-0/FF/FK decomposition and
  the C_ell version), whose wide-angle FK branch would change shape
* the arcminute-scale FK amplitude, which survives the grid fix at the
  0.15% level and is the one number that does NOT change

None of this is a claim that the physics of the FK term is wrong. The
diagrammatics, the measure bookkeeping and the small-angle amplitude all
stand. What changes is the angular shape at wide separations and the
status of the amplitude as a cutoff-conditioned number.

## Converged-tree run (2026-08-25): the shape conclusion is cutoff-robust

A tree table was rebuilt on the dense collapsed grid with the HIGH branch
summed over seven octave bands up to ell = 8192 (617 s; each band carries
its own 96-node quadrature so the density in ln ell stays uniform, which
is what makes a high cutoff affordable at all). Folded through the
production pipeline and cross-checked against the independent diagram
assembly: median ratio 0.9982, spread 0.0043.

xi_FK for kk at gamma = 0.5 arcmin:

| table | xi_FK | vs dense@1000 |
|---|---|---|
| deployed grid, ell_max 1000 | +4.288e-06 | 1.00 |
| dense grid, ell_max 1000 | +4.281e-06 | 1.00 |
| dense grid, ell_max 2000 | +8.433e-06 | 1.97 |
| dense grid, ell_max 8192 | +1.738e-05 | **4.06** |

The growth is decelerating, 1.97 for the first doubling and 1.43 per
doubling from 2000 to 8192, consistent with the band slopes; ell = 8192 is
therefore closer to convergence than 1000 but not converged, and the tree
amplitude will keep climbing by tens of percent per doubling for a while.

**The important result is that the angular SHAPE does not move.**
|FK / Order-0| in percent:

| gamma ['] | deployed | dense@1000 | dense@2000 | dense@8192 |
|---|---|---|---|---|
| 0.50 | 0.508 | 0.507 | 0.999 | 2.059 |
| 3.31 | 0.773 | 0.732 | 1.287 | 1.695 |
| 21.88 | 3.422 | 0.697 | 0.996 | 1.423 |
| 56.27 | 14.990 | 0.781 | 1.001 | 1.340 |
| 144.71 | 192.487 | 0.526 | 0.835 | 1.002 |

Every dense column is flat at the one-to-two percent level across three
decades in angle, while the deployed column climbs by three orders of
magnitude. Raising the cutoff rescales the FK amplitude but does not
restore the wide-angle rise, because the cutoff acts on the hard pair and
the angular dependence comes from the soft leg.

So the two questions separate cleanly:

* **amplitude** is cutoff-dependent and currently underconverged; the
  converged tree value at 0.5 arcmin is at least 4x the published one, and
  a nonlinear spectrum multiplies that further
* **shape**, and with it the claim that FK overtakes Order-0 near
  2 degrees, is a property of the cosine grid alone, and is unchanged by
  the cutoff. FK tracks Order-0 at one to two percent and never crosses it.
