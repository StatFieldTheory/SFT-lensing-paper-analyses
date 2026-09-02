# ell_high_max resolved: a convergence parameter, not a physical scale

> **Correction, 2026-09-02.** The absolute FK amplitudes in this note are
> `ell_max = 1000` values obtained with the sorting callable and the production
> linear-in-lambda fold. The manuscript's FK is evaluated at `ell_max = 15360`
> with the permutation-aware callable, where the same quantity is `+1.9480e-5`
> at `gamma = 0.5'` (2.31% of Order-0), a factor 4.5 larger. Ratios between
> variants at a fixed cutoff are unaffected. See
> [`../FK_BASELINE_NUMBERS.md`](../FK_BASELINE_NUMBERS.md) for the full key and
> `sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/` for the
> current product.

Measured 2026-08-25 with `code/probe_uv_convergence.py`; raw output in
`products/uv_convergence_bands.txt`. This supersedes the open question left
by `finding_ell_max_uv_sensitivity.md`, and corrects the guess recorded
there that the nonlinear case might be genuinely divergent.

## Why the earlier evidence was ambiguous

Sweeping the cutoff measures a partial sum, and over one decade a slowly
climbing convergent sum is indistinguishable from a divergent one. Both
looked like growth. The differential contribution settles it: evaluate the
HIGH branch on successive octave bands (L, 2L], each with its own
quadrature so node density is uniform by construction, and read the local
slope p in band ~ L^(-p). Positive p means the tail is dying.

## Result: both models converge, far above the deployed cutoff

Local slope of the octave-band contribution to zeta_TTT at gamma = 0.5':

| band | tree, 396.6 Mpc | tree, 1200 | tree, 2326.6 | BiHalofit, 396.6 | BiHalofit, 1200 | BiHalofit, 2326.6 |
|---|---|---|---|---|---|---|
| (250, 500] | -0.195 | -1.214 | -2.996 | -2.135 | -1.691 | -2.996 |
| (500, 1000] | +0.146 | -0.626 | -2.221 | -1.766 | -1.824 | -2.223 |
| (1000, 2000] | +0.388 | -0.167 | -1.418 | -1.187 | -2.025 | -1.430 |
| (2000, 4000] | +0.567 | +0.186 | -0.761 | -0.640 | -1.737 | -0.820 |
| (4000, 8000] | +0.755 | +0.498 | -0.185 | -0.127 | -1.064 | -0.433 |
| (8000, 16000] | +1.098 | +0.959 | +0.525 | +0.527 | -0.131 | -0.098 |

The slopes rise monotonically toward positive values in every column. The
tree vertex turns over between ell ~ 1000 (near shell) and ell ~ 8000 (far
shell). The BiHalofit vertex turns over near ell ~ 8000 at the near shell
and is only just reaching its turnover at ell ~ 16000 for the middle and
far shells, where the slope is -0.13 and -0.10.

**So ell_high_max is a numerical convergence parameter, not a physical
smoothing scale.** The vertex integral converges; the deployed value of
1000 simply sits on the RISING side of the integrand, below its peak.

The mechanism is the effective spectral slope. The collapsed vertex
integrates the squeezed bispectrum over the hard pair's transverse
wavenumber, so a band contributes as INT dl l P(l/chi) ~ L^(2 + n_eff)
with n_eff = dlnP/dlnk. Measured on the fiducial spectra:

| k [h/Mpc] | n_eff linear | 2 + n_eff | n_eff nonlinear (z=1) | 2 + n_eff |
|---|---|---|---|---|
| 0.1 | -1.65 | +0.35 | -1.55 | +0.45 |
| 1.0 | -2.29 | -0.29 | -1.09 | +0.91 |
| 3.0 | -2.48 | -0.48 | -1.72 | +0.28 |
| 10.0 | -2.61 | -0.61 | -2.10 | -0.09 |

The linear spectrum crosses n_eff = -2 near k ~ 0.5 h/Mpc, so tree turns
over early. Halofit's one-halo term FLATTENS the spectrum, pushing n_eff
back up to about -1.1 at k ~ 1 h/Mpc, so the nonlinear turnover is delayed
to k ~ 5 to 10 h/Mpc. That is exactly the ordering the band measurement
shows, and it is why the nonlinear case looked divergent over the first
decade.

## How far from converged the deployed value is

Summing the measured bands, comparing the total up to ell = 16000 with the
part up to ell = 1000:

| gamma, shell | tree | BiHalofit |
|---|---|---|
| 0.5', 396.6 Mpc | x1.56 | x10.1 |
| 0.5', 1200 Mpc | x2.64 | x57.3 |
| 0.5', 2326.6 Mpc | x13.8 | x17.7 |
| 5.45', 1200 Mpc | x1.56 | x3.87 |
| 42.2', 2326.6 Mpc | x20.5 | x25.1 |

These are lower bounds on the converged value: the BiHalofit columns at
the middle and far shells are still rising at ell = 16000, and the bands
below 125 are not included. The tree vertex is undercounted by factors of
1.5 to 20 depending on shell and angle, the nonlinear one by 2 to 57.

## What this means

1. The deployed FK amplitude is not "conditioned on an arbitrary cutoff".
   It is UNDERCONVERGED, and the converged value is substantially larger.
   The factor 1.97 measured earlier for a cutoff doubling is simply the
   partial sum climbing toward its peak.
2. The fix is mechanical rather than conceptual: raise ell_high_max until
   the band contribution has turned over and decayed, which the slopes put
   at ell of order 30000 or more for the nonlinear model at the far shell.
   Cost scales as the number of quadrature nodes at constant density, so
   this is expensive but not difficult.
3. There is a genuine physical caveat, and it is sharper than the cutoff
   question. The vertex is dominated by hard modes at k of several h/Mpc,
   which is where the matter bispectrum model is least trustworthy:
   BiHalofit is calibrated to k < 3 to 10 h/Mpc, and baryonic effects are
   significant above k ~ 1 h/Mpc. So the converged FK amplitude will rest
   on the part of the bispectrum with the largest model error, and the
   baryonic boost stops being a small correction there. That is the honest
   statement to attach to any converged FK number.
4. A practical constraint surfaced during the measurement: the fiducial
   CAMB table stops at k = 206 h/Mpc, and the third leg reaches 2 x
   ell_max, so at the nearest shell the tree branch runs out of table above
   ell ~ 27000. A converged production build needs a longer P(k) table or
   an explicit high-k extrapolation.

## Correction to an earlier note

`finding_ell_max_uv_sensitivity.md` reasoned from the coincident-pair
moment that the vertex might be genuinely UV-divergent, and this note's
first draft repeated that for the nonlinear case on the strength of the
n_eff estimate at k ~ 1 h/Mpc. The band measurement shows otherwise: the
radial collapse leaves a two-dimensional transverse integral, and that
integral converges once n_eff drops below -2, which it does for both
spectra, only later for the nonlinear one. The divergence argument was
right about the mechanism, the coincident hard pair with no angular
damping, and wrong about the conclusion.

## Correction and sharpening (2026-08-25, later): why the scale is so large

The estimate of "ell of order 30000" above is too small, and the reason is
worth stating precisely, because the number is counter-intuitive until one
sees that it is not a resolution scale at all.

After the radial collapse the hard pair contributes as a two-dimensional
transverse integral, INT dk k P(k). Its integrand peaks where
n_eff = dlnP/dlnk equals -1, and it converges only once n_eff < -2. So the
convergence scale is set by where the spectrum steepens past -2, and the
approach is a power law, not an exponential cut.

Measured cumulative fraction of INT dk k P(k) at z = 1:

| | 50% below | 90% below | 99% below |
|---|---|---|---|
| linear | k = 0.54 | 8.0 | 82 h/Mpc |
| halofit | k = 9.5 | 80 | 180 h/Mpc |

In multipoles, ell = k * chi_h at the three shells
(chi_h = 293 / 4445 / 5582 Mpc/h):

* linear, 90%: ell = 2300 (near shell) to 45000 (far shell)
* halofit, 90%: ell = 23000 to 450000

The cause is the one-halo term. Halofit flattens the spectrum so that
n_eff only reaches -1.09 at k = 1 h/Mpc and -2.10 at k = 10 h/Mpc; beyond
that the tail of INT dk k P falls as k^(1+n_eff) ~ k^(-0.1), so essentially
every decade keeps contributing. The linear spectrum crosses -2 near
k = 0.5 h/Mpc and therefore converges at multipoles of a few thousand.

### The consequence is physical, not numerical

Chasing ell ~ 10^5 would be the wrong response. What the measurement says
is that, as currently formulated, the FK vertex is dominated by modes at
k = 10 to 100 h/Mpc, which is inside single haloes. There:

* BiHalofit is used an order of magnitude beyond its calibration (k < 10)
* baryonic physics dominates the matter spectrum entirely
* the coincident-pair moment is UV-dominated by construction, the same
  reason the nonlinear density variance is set by the smallest scales
  retained

So the honest treatment is to report the FK amplitude as a function of the
small-scale cutoff, and to state that the tree-level case converges by
ell of a few thousand while the nonlinear case does not converge on any
scale the model controls. A converged nonlinear FK number would be a
statement about the one-halo regime, not about the bispectrum physics the
draft is describing.

## A structural floor on lambda that the P(k) table imposes

Discovered while building a table with shells extended toward the observer
(2026-08-25). The Limber branch evaluates the third leg at k = 2 ell_max /
chi_h, so as the shell approaches the observer the required wavenumber
diverges. With the fiducial table's k_max = 206 h/Mpc and ell_max = 1000,
the innermost usable shell is

    chi_h = 2 ell_max / k_max = 9.7 Mpc/h,  i.e. chi = 14.5 Mpc.

A build at lambda = 6 Mpc fails loudly with "k values above k_max", which
is the correct behaviour: canoes refuses to extrapolate P(k) silently.

Two consequences. First, the coupling between ell_max and the smallest
reachable shell is real: raising the cutoff also pushes the innermost
usable shell outward, unless the P(k) table is extended. Second, the
deployed floor at lambda = 396.6 Mpc is NOT set by this constraint, since
14.5 Mpc would be allowed; it is inherited from the shared L2 grid.
