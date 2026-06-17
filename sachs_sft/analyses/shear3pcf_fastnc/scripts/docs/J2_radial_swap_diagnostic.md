# Diagnostic A: decouple radial from angular in the J_2 shear-3PCF residual

Builder: `scripts/probe_J2_radial_swap.py`
Output:  `outputs/J2_radial_swap.npz`
Interpreter: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`

## TL;DR -- VERDICT

The surviving phi-dependent ~factor-2 J_2 residual `[0.68 @30, 0.43 @60,
2.03 @120]` (gamma=10') is **100% ANGULAR/CHANNEL, NOT radial.** Feeding the
canoes B_delta radial shape through the PROVEN-correct J_2 orientation kernel
reproduces NONE of the phi-trend: it gives a perfectly FLAT 0.453x rescale.

> The radial ell-shape difference between canoes B_delta and the B_kappa-phase
> reference is, AFTER the three committed canoes fixes (h^4, LOW group-zeroing,
> (1+z)^-4 jacobian), a **constant ~0.453** with NO ell-slope (the documented
> +0.44 slope is GONE -- it was a symptom of the now-fixed (1+z) radial measure).
> A constant radial rescale factors out of the linear J_2 assembly, so it CANNOT
> produce a phi-trend. The [0.68,0.43,2.03] swing therefore lives entirely in the
> canoes great-circle spin-2 cumulant route (the angular/channel projection),
> NOT in the B_delta radial integrand.

## Method (clean decoupling)

The B_kappa-phase reference `natural_components` (`stage1_bkappa_phase_reference.py`)
IS the proven-correct flat-sky spin-2 J_2 assembly: it takes a radial input
`B(l1,l2,l3)` and folds it through the centroid-projected 2D-Fourier orientation
integral `exp[2i s_j(phi_lj - alpha_j)]` (== `2 pi i^S J_S(|Q|) e^{iS arg Q}`
after the alpha-integral; `zetaD_HIGH_fix.md` Gate B verified it to <1%). I feed
this SAME assembly two radial inputs, holding the angular kernel BYTE-IDENTICAL:

  (a) `refB` = `B_kappa(l1,l2,l3)`  -- the trusted reference radial shape;
  (b) `canB` = `W_can(l1,l2,l3)`    -- the canoes HIGH convergence-3pt scalar
       weight (the deployed zeta_D radial integrand, PPP response), built as a
       general triangle function exactly as `probe_HIGH_B_ratio_ellscaling.W_can`.

Identity check: my assembly's `fi_ref` reproduces the stored FINAL `fi_ref` to all
digits (`1.975916e-07, 2.264624e-07, 1.280506e-07`), so the kernel is the same
object that built the trusted reference -- only the radial function is swapped.

## Numbers (gamma=10', phi=30/60/120)

### Radial swap: canB-through-J2 / refB-through-J2 (frame-invariant)

| phi  | fi_ref(refB,J2) | fi_can(canB,J2) | ratio canB/refB | per-comp [G0,G1,G2,G3] |
|------|-----------------|-----------------|-----------------|------------------------|
| 30   | 1.975916e-07    | 8.957289e-08    | **0.4533**      | [0.453,0.453,0.454,0.453] |
| 60   | 2.264624e-07    | 1.026765e-07    | **0.4534**      | [0.454,0.453,0.453,0.453] |
| 120  | 1.280506e-07    | 5.796493e-08    | **0.4527**      | [0.453,0.452,0.452,0.455] |

Normalized to median: `[1.000, 1.00015, 0.99856]` -- FLAT (swing 1.0016x).
Per-component flat across ALL components (0.452-0.455 everywhere): the radial
swap does NOT reshuffle the channel split either.

### Equilateral sanity: R = W_can/B_kappa vs ell

`R = [0.4446, 0.4489, 0.4535, 0.4534, 0.4535]` at `L = [60,100,200,500,1000]`;
log-log slope = **+0.0065** (FLAT). The stored `Wcan_vs_Bkappa_decomp.npz`
+0.44 slope and ~5-18x magnitude were from a PRE-FIX canoes (radial jacobian
`ones`); after `1ddc141` ((1+z)^-4) the slope is gone and R is a flat ~0.45.

### Decisive comparison

| quantity                              | phi=30 | phi=60 | phi=120 | swing |
|---------------------------------------|--------|--------|---------|-------|
| FULL great-circle route can/ref       | 0.4755 | 0.2966 | 1.4152  | --    |
| FULL residual (can/ref)/median = TARGET| 0.6808 | 0.4246 | 2.0262  | 4.77x |
| **RADIAL-only canB/refB / median**     | 1.0000 | 1.0002 | 0.9986  | **1.0016x** |
| leftover = TARGET / radial_norm       | 0.6808 | 0.4245 | 2.0291  | 4.78x |

**The leftover (after dividing the radial-only trend OUT of the target) is the
ENTIRE target, unchanged (swing 4.78x ~ target 4.77x).** Radial explains 0% of
the phi-trend.

## What this localizes

- **Angular integral / J_2 kernel: EXONERATED for good.** It reproduces the
  reference byte-for-byte and, fed the canoes radial shape, produces NO phi-trend.
- **Radial B_delta ell-shape: EXONERATED as the source of the phi-trend.** After
  the fixes it is a flat ~0.453 constant -- it sets the residual OVERALL SCALE
  (~0.45, consistent with the median ~0.70 once the full-route channel mixing is
  included) but carries ZERO phi structure.
- **The phi-dependent ~factor-2 residual is in the canoes GREAT-CIRCLE spin-2
  cumulant route** -- specifically the per-leg spin-2 great-circle -> centroid
  projection / channel-split (`_spin_aware_three_pt.py` canonical-frame harmonic
  sum + the per-leg conjugation/helicity assignment that maps Gamma^0..3 to
  zeta_PPP / zeta_D). That is where the full-route `fi_can/ref = [0.475, 0.297,
  1.415]` per-config (and the per-component reshuffle G1,G2 over / G3 under at
  phi=120) is generated, since the flat-sky J_2 assembly fed the same radial
  shape does NOT produce it.

This is consistent with, and sharpens, `shear3pcf_natural_map.md` caveat 1
(great-circle Wigner-3j/Limber projection != fastnc fixed-frame Fourier spin-2
projection) and the channel-label/leg-conjugation hint (the fnc/ref Gamma^3
reshuffle 2.40x): the residual is a PROJECTION/CHANNEL effect of the
great-circle spin-2 route, not a B_delta radial-shape effect.

## Next step (scoped to the angular/channel route, NOT the kernel or radial)

Derive (xAct) the per-leg great-circle -> centroid spin-2 projection +
channel-conjugation map for the canoes `_spin_aware_three_pt` canonical frame,
and diff it against the flat-sky centroid `exp[2i s_j(phi_lj-alpha_j)]` the
reference/fastnc use. The residual is the difference between these two spin-2
PROJECTIONS of the same B_kappa, evaluated per natural component.
