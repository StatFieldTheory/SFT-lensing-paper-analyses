# Harmonic space, and how far to trust an underconverged loop

Two things this note settles: what the corrected FK looks like as an
angular power spectrum, and what the ell_max non-convergence does and does
not undermine.

## The C_ell picture inverts as well

Transformed with the draft's own curved-sky Wigner-d projection
(`analysis3/plot_analysis3_cl_decomposition.py`, kappa-kappa through
d^l_{0,0} = P_l, same apodisation and DC subtraction), |C_l^FK / C_l^O0|:

| ell | deployed | dense@1000 | dense@2000 | BiHalofit@1000 | dense@8192 |
|---|---|---|---|---|---|
| 3 | **17.95** | 0.197 | 0.450 | 0.278 | 0.085 |
| 7 | 3.81 | 0.084 | 0.189 | 0.142 | 0.108 |
| 17 | 0.241 | 0.029 | 0.013 | 0.050 | 0.034 |
| 39 | 0.0008 | 0.0088 | 0.0082 | 0.0126 | 0.0007 |
| 93 | 0.0022 | 0.0079 | 0.0107 | 0.0103 | 0.0127 |
| 514 | 0.0001 | 0.0088 | 0.0126 | 0.0150 | 0.0171 |
| 1211 | 0.0003 | 0.0003 | 0.0144 | 0.0007 | 0.0204 |

The ell dependence is not merely rescaled, it is REVERSED. The deployed FK
band power falls with ell and is concentrated at ell < 10, where it exceeds
Order-0 by up to a factor 18 and crosses unity near ell = 9. Every dense
variant instead RISES gently with ell, sitting at about 1 to 2% of Order-0
from ell = 40 out to 1500, and never crossing unity anywhere.

The corrected behaviour is the physically expected one. FK is sourced by
the squeezed bispectrum's hard modes, so its power should live at small
angular scales, not be concentrated in the first few multipoles. A flat
band power at low ell followed by a collapse is the harmonic-space
signature of a real-space curve that is constant out to 10 degrees, which
is what the frozen corner value produces.

Caveat on the lowest multipoles: the ell < 20 points are unstable, bouncing
between variants (0.197, 0.450, 0.278, 0.085 at ell = 3). They come from
transforming the wide-angle tail of xi, beyond 10 degrees, which
`finding_grid_artifact_measured.md` already flags as oscillatory and
cancellation-dominated. Read the ell >= 40 plateau, not the first points.

This bears on the recorded low-ell result from 2026-06-11, where the FK
excess over Order-0 at ell < 10 was measured at 18x with a crossover near
ell = 9 to 10, and reconciled with the CMB-lensing literature by noting
that at ell = 100 to 1000 the ratio was only 0.14 to 0.24%. On the dense
grid the 18x becomes about 0.2, i.e. no crossover, while the ell = 100 to
1000 ratio rises from ~0.2% to ~1%. Both halves of that comparison move.

## Is an underconverged loop still credible?

The FK term is a loop: the K vertex supplies the three-point cumulant and
the fold integrates over the internal multipole of the coincident hard
pair. Whether such a loop converges is a property of the input spectrum,
not of this formalism, and the same question arises for the post-Born and
lens-lens couplings computed routinely in weak lensing.

What the measurement established:

1. The integral does converge. The octave-band slopes turn positive for
   both spectra (`finding_ell_max_resolved.md`). It is not divergent.
2. It converges slowly, because the hard pair contributes as
   INT dk k P(k) and the nonlinear spectrum only steepens past n_eff = -2
   near k = 10 h/Mpc. Ninety percent of the integral needs k ~ 8 h/Mpc for
   the linear spectrum and k ~ 80 h/Mpc for halofit.
3. The deployed ell_max = 1000 therefore undercounts. The tree amplitude at
   0.5 arcmin rises by 4.06 from ell_max 1000 to 8192 and is still climbing
   at about 1.43 per doubling.

So the amplitude is a number that must be quoted with its k_max, exactly as
a loop integral should be. That is a limitation of the prediction, and it
should be stated, but it is not a reason to disbelieve the calculation.

What it does NOT undermine, and this is the part that matters for the
conclusions drawn here:

* **The shape.** |FK / Order-0| is flat at one to two percent across three
  decades in angle for cutoff 1000, 2000 and 8192 alike. The cutoff acts on
  the hard pair; the angular dependence comes from the soft leg. So the
  finding that FK does not overtake Order-0 near 2 degrees is independent
  of where the cutoff sits.
* **The grid diagnosis.** That rests on the deployed table's own two
  samples at gamma = 0 and 31 degrees, and on a like-for-like comparison at
  fixed cutoff.
* **The machinery.** The independent diagram assembly agrees with the
  production fold to 0.1 to 0.4% on every table tested.

The honest summary is that this study has made the FK prediction MORE
credible in shape and LESS certain in amplitude: the angular behaviour is
now measured rather than interpolated, while the amplitude is revealed to
be an underconverged loop whose converged value depends on the
one-halo-regime bispectrum. Both statements should travel together.

## Consistency point worth checking

The two-point side of the pipeline (corr_op) is built to ell_max = 5000,
while the three-point vertex stops at 1000. Any FK-versus-FF comparison in
the draft therefore compares terms truncated at different scales. This has
not been quantified here.
