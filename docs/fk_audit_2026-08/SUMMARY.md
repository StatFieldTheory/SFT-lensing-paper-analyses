# Nonlinear driving-field statistics for the SFT lensing FK term

Session of 2026-08-25. Start at `README.md` for the folder map; this file
is the conclusion of the work.

The question was whether a published emulator can supply the three-point
input statistics of the SFT lensing driving fields, and then to build the
code that redoes the three-point leakage into the two-point function with
that input. Both were done. Along the way the study measured two
properties of the existing pipeline that matter more than the upgrade
itself.

## 1. The emulator landscape

No emulator of the absolute gravity-only matter bispectrum exists. The
only genuine three-point emulator in the Euclid ecosystem is the baryonic
BOOST emulator of Burger et al. (arXiv:2506.18974), shipped as
`baccoemu.Matter_bispectrum`, which multiplies a gravity-only model rather
than replacing it. The community backbone is still the BiHalofit fitting
formula, which is also the best-performing formula in the squeezed
configurations this vertex samples. There is no Weyl or Ricci field
emulator, and none is needed: the driving fields are per-multipole linear
transfers of one scalar potential, so a matter bispectrum converts exactly.
Details and the full negative results: `notes/emulator_candidates.md`.

## 2. The conversion, specified and implemented

`notes/input_ready_b_to_zeta_spec.md` is the source-verified map from
B(k1,k2,k3;z) to the configuration-space six-channel zeta table and on to
the sft-wick coupling. Two facts made the implementation simpler than
expected: under Limber the screen-Hessian factor cancels the Poisson
factor per leg, so the kernel is a scalar; and B is evaluated once per
shell on the quadrature grid, shared across all cosine triples, so the
whole table costs about 9.4 million bispectrum calls, which is 8 seconds.

Implemented as `code/b_model` (tree, BiHalofit, BiHalofit x baryon ratio)
plus an optional `b_delta_fn` on the two canoes HIGH-branch functions. The
tree branch reproduces the canoes internal assembly to 6.6e-13 on its own
quadrature geometry; through the full vertex the residual is 3.7e-6, all
of it the growth-factor interpolation route, removable via `growth_fn`.
59 tests here, and the canoes suite passes 341 with no failures.

## 3. What the measurements found

Four table variants were folded through the production FK pipeline, each
changing one variable. At gamma = 0.5 arcmin, xi_FK for the kk pair:

| variant | xi_FK | vs dense tree |
|---|---|---|
| deployed grid, tree, cutoff 1000 | +4.287723e-06 | x1.0015 |
| dense grid, tree, cutoff 1000 | +4.281487e-06 | 1 |
| dense grid, tree, cutoff 2000 | +8.432787e-06 | **x1.97** |
| dense grid, BiHalofit, cutoff 1000 | +7.089177e-06 | **x1.66** |

(The absolute values in this table read the shell grid through the
production lambda interpolation; section 4d records the correction. The
ratios between variants are unaffected.)

**The deployed vertex is undercounted because its multipole cutoff sits
below the integrand's peak.** Doubling ell_high_max nearly doubles the FK
amplitude, more than the tree-to-BiHalofit upgrade delivers (1.66).
Resolving the differential contribution octave by octave shows the
integral does CONVERGE, for both spectra: the band slope rises
monotonically and turns positive near ell ~ 1000 to 8000 for tree and
ell ~ 8000 to beyond 16000 for BiHalofit, the delay being halofit's
one-halo term flattening the spectrum so that n_eff crosses -2 only at
k of several h/Mpc. Summing the bands, the deployed cutoff of 1000
undercounts the vertex by factors of 1.5 to 20 (tree) and 2 to 57
(nonlinear). So the cutoff is a convergence parameter, not a physical
smoothing scale, and the fix is mechanical: raise it until the bands have
turned over. The physical caveat that remains is sharper: the converged
vertex is dominated by k of several h/Mpc, where the bispectrum model is
least trustworthy and baryons are not a small correction.
(`notes/finding_ell_max_resolved.md`.)

**The deployed cosine grid freezes the vertex at its corner.** The FK fold
queries only the collapsed family, and on the 15-value uniform grid the
interpolation puts at least 96% of its weight on the (1,1,1) cell for 27
of the 40 production separations. The deployed FK curve is consequently
flat from 0.5 to about 600 arcmin, while the true vertex decays by three
orders of magnitude and changes sign near 183 arcmin. The arcminute
amplitude is unaffected at the 0.15% level.
(`notes/finding_grid_artifact_measured.md`.)

The grid conclusion does not depend on the second table. The deployed
cosine spacing is 1/7, so its first two collapsed-family samples are at
gamma = 0 and gamma = 31 degrees, while the paper's FK figure spans 0.008
to 33 degrees: the entire curve lies in the first grid gap. Interpolating
the deployed table's own rows in one dimension, with the k-nearest-neighbour
lookup bypassed, reproduces the same flat curve, and its own next sample is
three orders of magnitude below the corner value. Both tables share the
same Limber code and the same measure; only the angular sampling differs.

**Keeping the exact-3j branch at tree level is safe.** Its FFTlog
machinery needs a separable bispectrum, so it cannot take a general model.
Measured on the shared band, that choice costs at most 1.5% of the vertex,
an order of magnitude below the BiHalofit model error.
(`notes/finding_ell_cut_seam.md`.)

## 4. Consequences for the draft, recorded not acted on

The draft was not modified and no draft figure was regenerated.
`notes/impact_on_draft_discussion.md` gives the assessment, verified
line by line against the text. In brief: the small-angle numbers the draft
quotes are confirmed (FK is 0.507% of Order-0 and 2.13 times FF at
0.5 arcmin); the FF statements are untouched, since FF never reads this
table (FF overtakes Order-0 at 183 arcmin, as the draft says); the
mechanism claim that the squeezed configuration routes small-scale
structure into the two-point function is confirmed and in fact stronger
than stated; but the three statements that FK overtakes Order-0 beyond
about 2 degrees do not survive the grid fix, because they read the frozen
corner value against a genuinely decaying Order-0. The redshift-trend
statement is scored in section 4d.

`figures/fk_variants_kappa.pdf` shows this in the draft's own visual
language.

## 4b. The three follow-up runs (2026-08-25, later)

**Dense grid at a converged cutoff.** The HIGH branch now supports octave
banding (`--ell-bands`): the cutoff range is split into disjoint windows in
max-leg, each with its own 96-node quadrature, so node density in ln(ell)
stays uniform. That is what makes a high cutoff affordable; one run to
ell = 8192 at constant density would need n_ell ~ 790 and a
(790, 790, 64) grid. Seven bands, 617 s. Result at gamma = 0.5 arcmin:
xi_FK = +1.738e-05, which is 4.06x the deployed value, still climbing at
about 1.43 per doubling, so ell = 8192 is closer to converged than 1000
but not converged.

**The shape conclusion is robust to the cutoff.** |FK / Order-0| in
percent stays flat at one to two across three decades in angle for EVERY
dense variant (1000, 2000, 8192), while the deployed column climbs from
0.5 to 192. Raising the cutoff rescales the amplitude; it does not restore
the wide-angle rise, because the cutoff acts on the hard pair while the
angular dependence comes from the soft leg. So the finding about the
2-degree crossover does not depend on the convergence question.

**Monte-Carlo cross-check, with a control that settles its scope.** The
independent diagram assembly in `analyses/mc_sachs_2pt/fk_analytic.py`
reproduces the production fold on the dense table to a median ratio of
0.9988 (spread 0.0012) across three orders of magnitude of decay, and on
the converged table to 0.9982. Pointed at the DEPLOYED table it agrees
equally well, 0.9987 with zero spread, on a curve that is flat instead of
decaying. The cross-check therefore validates the assembly and the fold,
and has no power to judge the input statistics, because both sides read
the same table. (`notes/finding_mc_crosscheck.md`.)

**The lambda_min floor is a 0.1% effect at z_s = 5.** Six of the 24
Gauss-Legendre nodes fall below the table's innermost shell at 396.6 Mpc
and are silently zeroed. Rebuilding with shells down to 50 Mpc changes
xi_FK by a median +0.12% below 700 arcmin, consistent with the K^4 kernel
estimate of 0.34%. The excised region is always z < 0.1, so its kernel
weight grows as the source redshift falls: 1.5% (K^2) at z_s = 5 against
12.7% at z_s = 0.5. That is a reason to distrust the multi-redshift
figure's trend, not the main analysis. (`notes/finding_lambda_min_cut.md`.)

## 4c. Harmonic space, and the credibility of an underconverged loop

Transformed with the draft's own curved-sky Wigner-d projection, the ell
dependence reverses rather than merely rescaling. |C_l^FK / C_l^O0| for the
deployed table is 17.95 at ell = 3, crosses unity near ell = 9 and collapses
to 1e-4 by ell = 200. Every dense variant instead rises gently with ell,
sitting at 1 to 2% from ell = 40 to 1500 and never crossing unity. The
corrected behaviour is the physically expected one for a term sourced by
the squeezed bispectrum's hard modes. The lowest multipoles are unstable
between variants and should not be read.

On credibility: the FK term is a loop, and whether it converges is a
property of the input spectrum, the same question that arises for post-Born
and lens-lens couplings. The measurement says it DOES converge, slowly,
because the hard pair contributes as INT dk k P(k) and the nonlinear
spectrum only steepens past n_eff = -2 near k = 10 h/Mpc. So the amplitude
must be quoted with its k_max, as a loop integral should be. The shape, the
grid diagnosis and the machinery checks are all independent of that.
(`notes/finding_cl_and_credibility.md`.)

## 4d. Redshift dependence, measured (2026-08-25, latest)

The redshift dependence is now explicit at both levels
(`notes/theory_redshift_dependence.md`). Per shell, the tree table
factorizes exactly into growth times geometry, and the full explicit
chain, three legs of A(a)(1+z)^4, the (1+z)^-4 collapse Jacobian,
chi^-4, D^4 and the unit conversions, is certified against the real
table to 0.5% across shells at z = 0.10 to 5.7. The nonlinear
enhancement of the vertex is a low-redshift-shell phenomenon (21.6x at
the z = 0.1 shell, 1.00 at z >= 3), and even the z_s = 5 fold takes half
its weight from shells below z = 0.88, so the late universe controls the
observable at every source redshift.

Folding the same dense tables at lambda_f(z_s) measures the source-
redshift trends directly: the tree FK/O0 fraction is flat in z_s to
-16%/+5% (0.40-0.50% at 0.5 arcmin at every z_s), while with nonlinear
input the fraction FALLS with source redshift (~2.7% at z_s = 0.5
against ~0.8% at z_s = 5, at ell <= 1000: the direction is robust, the
absolute values are cutoff-conditioned), the opposite direction from the
deployed
multi-z figure's impression. Four distortions of the deployed trend are
now measured rather than estimated: the cosine-grid artifact (dominant
at wide angles at every z_s), the fixed-ell cutoff (fold-weighted
undercount 1.39 at z_s = 0.5 vs 2.14 at z_s = 5, a 1.5x tilt against
high z_s), the lambda_min floor (-6.6% at z_s = 0.5, -0.1% at z_s = 5),
and a NEW one: the production fold's linear-in-lambda interpolation of
the sparse mid-path shell grid inflates xi_FK at z_s = 5 by +11% at
arcminutes, decaying to +4% at ~1 degree and -10% at 145 arcmin
(truth-probed at mid-gap shells; the grid-converged tree baseline is
+3.84e-6, not +4.28e-6; all variant RATIOS in sections 3 and 4b are
unaffected). The symbolic layer
(`code/verify_redshift_scalings.wl`, 15 checks PASS) pins the power-law
exponents, including the low-redshift fold anchor
FK/O0 ~ chi_s^(-(n+1)), gamma-independent.

## 5. What should happen next, and why not now

The obvious next step, rebuilding the production 1671-triple table with
the collapsed family densified, is deliberately NOT done. Its output would
be an amplitude conditioned on a cutoff that changes it by a factor of two
per doubling, so the expensive build would have to be repeated once that
question is settled. The order that makes sense:

1. Raise the cutoff until the octave bands have turned over and decayed,
   which the slopes put at ell of order 30000 for the nonlinear model at
   the far shell, and FIX IT IN k RATHER THAN ELL: section 4d shows a
   fixed-ell cutoff undercounts the shells z_s-dependently (fold-weighted
   1.39 at z_s = 0.5 against 2.14 at z_s = 5) and so redshift-biases any
   multi-z result. This also needs a longer P(k) table or an explicit
   high-k extrapolation: the fiducial CAMB table stops at k = 206 h/Mpc
   and the third leg reaches 2 x ell_max, so the tree branch runs out of
   table above ell ~ 27000 at the nearest shell.
2. Rebuild the production table on the densified collapsed family,
   keeping the open triangles the order-1 three-point vertex needs, and
   fix the SECOND grid dimension found in section 4d: densify the
   lambda shells across the mid-path gaps (or make the callable
   interpolate in log-lambda), since the current linear interpolation
   inflates the z_s = 5 fold by +11% at arcminutes.
   `code/probe_lambda_interp.py` is the acceptance test.
3. Then, and only then, compare tree against BiHalofit as a physics
   statement, and revisit the draft's wide-angle passages with numbers
   that are no longer cutoff-conditioned.

The redshift trend no longer blocks on any of this: it is measured at
ell <= 1000 (section 4d), and only its converged-cutoff repeat waits for
steps 1-2. `HANDOVER.md` carries the updated per-task status.

Everything needed for steps 2 and 3 exists here: the build takes
`--b-model`, `--triples-npz` and `--cache-low`, the variant runner is
isolated from production artifacts, and the comparison and plotting
scripts consume the outputs directly.

## 6. Reproducing

```bash
CAN=/Users/zzhang/projects/angular_statistics/canoes
cd driver_field_emulators/code
PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python -m pytest b_model/tests -q
PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python sweep_ell_high_max.py
PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python probe_ell_cut_seam.py
```

Table builds, variant folds and figures are documented in `code/README.md`.
Everything this session touched outside the folder, and how to revert it,
is in `CHANGES_OUTSIDE_THIS_FOLDER.md`.
