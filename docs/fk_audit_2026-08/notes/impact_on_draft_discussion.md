# Does this change the draft's qualitative discussion?

Assessment only. The draft was not modified, and no draft figure was
regenerated. Line references are to the current text, read read-only.

Measured inputs: the four FK variants of
`finding_three_effect_decomposition.md`, plus the deployed Order-0 and FF
runs, all on the same 40-point gamma grid.

## The numbers the draft quotes at small angles are CONFIRMED

At gamma = 0.5 arcmin, with the corrected (dense) grid and tree-level B:

| quantity | measured | draft's wording |
|---|---|---|
| xi_FK | +4.2815e-06 | "a few times 1e-6" (insights) |
| xi_FK / Order-0 | 0.507% | "about 0.5% of Order-0" |
| xi_FK / FF | 2.13 | "roughly twice the FF term" |

(These confirmations are at the production linear-in-lambda fold. The
grid-converged values are +3.84e-6, ~0.45% and ~1.9x, which still
satisfy the draft's wording; FF never reads this table and is not
rescaled. See `theory_redshift_dependence.md` section 2.2.)

All three survive the grid fix (which is a 0.15% effect there). With
BiHalofit instead of tree the same numbers become 7.09e-06, 0.84% and
3.53.

## The FF statements are untouched

FF carries no K vertex and never reads the kappa3 table, so nothing in
this study bears on it. Measured independently: FF overtakes Order-0 in
absolute amplitude at 183 arcmin, matching the draft's "beyond about 200
arcmin" (insights.tex:110). FF's own shape is genuinely a slowly varying
plateau: relative to its 0.5-arcmin value it falls to 0.29 by 145 arcmin
and then flattens (0.2886 at 145', 0.2874 at 372'). The draft's
description of FF is accurate.

## What does NOT survive: the wide-angle FK claims

Three passages state that FK overtakes Order-0 beyond about 2 degrees:

* insights.tex:121 (caption of Fig. `NLO 2pt FFFK`): "overtakes Order-0 in
  absolute amplitude beyond gamma ~ 2 degrees"
* insights.tex:195: "on a z_s = 5 source plane it overtakes the falling
  Order-0 signal beyond ~2 degrees"
* conclusion.tex:36: "overtakes the linear signal in absolute amplitude
  beyond roughly 2 degrees"

Measured with the corrected grid, FK/Order-0 stays flat at 0.4 to 1.5%
across the entire range and never crosses unity in the physical range. The
crossing at 145 arcmin (2.4 degrees) is reproduced exactly by the deployed
table, and disappears when the cosine grid is densified.

The mechanism is visible in the shape. Normalised to its own 0.5-arcmin
value, the FK curve behaves as:

| gamma ['] | FF | FK deployed | FK dense (true) |
|---|---|---|---|
| 0.50 | 1.0000 | 1.0000 | 1.00000 |
| 21.88 | 0.3814 | 0.9997 | 0.20399 |
| 56.27 | 0.3070 | 0.9978 | 0.05208 |
| 144.71 | 0.2886 | 0.9853 | 0.00270 |
| 372.19 | 0.2874 | 0.9080 | -0.00077 |

The deployed FK is flatter than FF, which is what let it overtake a
decaying Order-0. That flatness is the frozen (1,1,1) corner value, not
physics. The true FK decays in step with Order-0, keeping their ratio
nearly constant.

Two further draft statements follow from the same passages and would need
revisiting on the same evidence:

* "where it exceeds FF" (insights.tex:121) holds at small angles (2.13x at
  0.5 arcmin) but reverses at large ones: at 145 arcmin FF is about 50
  times FK.
* "it strengthens steadily with source redshift (Fig. multiz kappa)"
  (conclusion.tex:37). NOW MEASURED (`theory_redshift_dependence.md`,
  sections 3-4): the absolute FK amplitude does grow steeply with z_s,
  but so does Order-0; the tree-level FRACTION FK/O0 is flat in z_s to
  -16%/+5% (0.40-0.50% at 0.5 arcmin at every source redshift), and with
  nonlinear input the fraction FALLS with z_s (~2.7% at z_s = 0.5
  against ~0.8% at z_s = 5). The figure's wide-angle impression that the
  leakage matters more at higher z_s is the cosine-grid artifact; the
  low-z_s deployed slice also overstates the tree FK by ~3.4x.

## What survives, and is in fact strengthened

The physical mechanism the draft describes is not in question. The
conclusion states that the three-point cumulant enters the two-point
function in its squeezed configuration, so that "small-scale structure has
a route into large-angle lensing statistics" and the result is "sensitive
to small-scale modes of the matter bispectrum" (conclusion.tex:30-35).

That sensitivity is real, and this study makes it sharper than the draft
does: the coincident hard pair of the collapsed configuration receives no
angular damping, so the vertex is UV-sensitive to the point where the
multipole cutoff sets its value (a factor 1.97 in xi_FK per doubling of
ell_high_max, `finding_ell_max_uv_sensitivity.md`). The route from small
scales into the two-point function exists exactly as claimed.

What does not follow from it is the observable consequence the draft
attaches: that the leakage dominates the signal at wide separations. The
sensitivity shows up as an amplitude that depends on the small-scale
cutoff, not as a term that outgrows Order-0 with angle.

Likewise untouched: the diagrammatic derivation, the selection rule that
sends FK into xi_kappa and xi_+ while cancelling it in xi_- and
xi_kappa_gamma, the measure and units bookkeeping, and the arcminute-scale
amplitude that all the quantitative statements rest on.

## Summary

| draft statement | status |
|---|---|
| FK ~ 0.5% of Order-0 at small angles | confirmed (0.507% at the production fold; ~0.45% grid-converged) |
| FK ~ twice FF at small angles | confirmed (2.13x at the production fold; ~1.9x grid-converged) |
| FF overtakes Order-0 beyond ~200' | confirmed (183'), unaffected |
| FK overtakes Order-0 beyond ~2 degrees | NOT confirmed: artifact of the cosine grid |
| FK exceeds FF | only at small angles; reverses by ~145' |
| squeezed configuration routes small scales into the 2PCF | confirmed, and stronger than stated |
| FK strengthens with source redshift | amplitude yes, fraction no: FK/O0 flat in z_s at tree level, falls with z_s for nonlinear input (`theory_redshift_dependence.md`) |
| FK amplitude as a prediction | conditional on ell_high_max (x1.97 per doubling) |

## Why the appendix Monte-Carlo agreed with the paper's FK

The appendix figure `fig: val mc` compares a direct Monte-Carlo of the
stochastic Sachs equation against the analytic workflow, and the two agree
on the FK channel. That agreement is not evidence that the FK curve is
right, because the Monte-Carlo reads the SAME table through the SAME
callable.

The cross-check's own code says so. `analyses/mc_sachs_2pt/fk_analytic.py`
states in its docstring that it "uses the SAME bare equal-time zeta
callable (equal_time_limber), the SAME scalar response and the SAME
F-tensor", so it "cross-checks the sft-wick FK diagram assembly (the
equal-time time-collapse, the causal measure, the F-K index contraction and
the cosine-triple geometry)", and explicitly that it "does NOT
independently re-derive the zeta table's internal units (both this
quadrature and analysis-3 read the same table)". `driver_stats.py`, which
supplies the driving statistics to the Monte-Carlo, imports
`equal_time_limber_kappa3_callable` directly.

So the validated scope is the machinery downstream of the table: the
diagram assembly and the fold. The input statistics are common to both
sides, so any property of the table, including the interpolation studied
here, cancels in the comparison. A Monte-Carlo seeded from a densified
table would agree just as well, at different values.

The appendix caption is consistent with this reading: it reports the
Monte-Carlo tracking "FK its decay toward the sign change at gamma = 32
degrees". That is the deployed curve's shallow decay, a factor of about
2.5 across the whole range with a sign change near 1944 arcmin, which is
what a lookup drifting slowly off the corner produces. The independently
computed vertex falls by three orders of magnitude over the same range.
