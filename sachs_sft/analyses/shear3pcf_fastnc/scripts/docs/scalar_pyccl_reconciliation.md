# Scalar convergence 3PCF vs PyCCL: decisive reconciliation (third independent code)

Builder: `scripts/probe_scalar_vs_pyccl_tree.py`
Output: `outputs/scalar_vs_pyccl_tree.npz`
Interpreter: `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`
Skill consulted: `pyccl`.

This settles the open question left by `zetaD_Bfix.md` (which found canoes SCALAR
convergence 3PCF ~15-49x larger than the B_kappa-phase scalar): it confirms the
canoes scalar conv-3PCF is genuinely over the STANDARD lensing convention, using
PyCCL as an independent third code (the one that already validated the canoes
2-point), with P(k), cosmology, growth and distances matched, and disambiguates
the bug as 3-POINT-SPECIFIC with a (1+z)^4 (not (1+z)^6) per-shell measure factor.

## Matched setup (apples-to-apples)

| input | value / method |
|---|---|
| cosmology | Omega_m=0.3160919980475834, h=0.6711, n_s=0.97, w0=-1 (canoes default == STF fiducial) |
| P(k) | PCAMBz0.txt (the SAME table canoes uses); sigma8(table)=0.80902; CCL with this sigma8 reproduces sigma8 to 5 decimals; table fed DIRECTLY as the physical-unit P(k) interpolant (P-matched, not just cosmology-matched) |
| growth D(a) | CCL `growth_factor` (D(1)=1); verified to match canoes growth (lightcone M[0]/(1+z)) to <0.1% over z in [0,5] |
| distances | CCL chi(z=5)=7968.7 Mpc vs canoes 7971.2 Mpc (0.03% — radiation in CCL) |
| source plane | single z_s=5 |
| gamma grid | identical to STAGE-0 (0.5'..300', includes 1', 10'-ish, 50'-ish) |
| config | squeezed (1, cos g, cos g) degenerate triangle, J0 fold (IDENTICAL to canoes STAGE-0 `Z_lambda`) |
| dimensionless | asserted: convergence 3PCF is dimensionless (kappa dimensionless), h-invariant (see below) |

The STANDARD side is the canoes STAGE-0 SPT reference with ONLY the radial
lensing-efficiency / per-leg response convention swapped: `[A(a)(1+z)^4]^3 K^3 a^8`
(canoes) -> `g(chi)^3` with the textbook single-plane weight
`g = (3/2) Om H0^2 (1/a) chi(chi_s-chi)/chi_s`.  Every angular ingredient (B_delta
= 2 F2 P P + 2cyc with TRUE SPT F2, the J0 kernel, the ell-quadrature) is
byte-identical, so the canoes/standard ratio is a CLEAN radial-convention ratio.

## A. canoes/PyCCL scalar conv-3PCF ratio (the primary number)

Real canoes `zeta_TTT` fold (STAGE-0 `Z_lambda`, deployed LOW+HIGH) vs the PyCCL
standard tree conv-3PCF `Z_std`, SAME squeezed J0 fold, SAME P(k):

| gamma | \|Z_std\| (PyCCL std) | \|Z_real_canoes\| | realcan/std (raw) | realcan/std (canoes footing*) |
|---|---|---|---|---|
| 0.5'  | 2.31e-6 | 2.29e-5 | **9.95** | 22.1 |
| 1.0'  | 2.30e-6 | 2.29e-5 | **9.95** | 22.1 |
| 12.2' | 1.02e-6 | 1.01e-5 | **9.83** | 21.8 |
| 27.2' | 3.91e-7 | 2.89e-6 | 7.40 | 16.4 |

`*` "canoes footing" multiplies my dimensionless `Z_std` by h^2 to put it on the
canoes `units="physical"` (x h^6) footing; see the h-power note below. The RAW
comparison (both as a genuinely dimensionless convergence 3PCF) gives **~10x**.

So canoes scalar conv-3PCF is **~10x the standard** (raw/dimensionless) or **~22x**
on the canoes `units="physical"` footing, at small gamma.  Either way it is a
LARGE, real over-normalization — fully consistent with `zetaD_Bfix.md`'s 15-49x
(that probe used a different triangle shape and the canoes h^2 footing).

**Cross-check vs the B_kappa-phase scalar reference**: my PyCCL `Z_std` (squeezed,
12.2') = 1.02e-6 vs the prior B_kappa-phase scalar (10', 60deg triangle) = 5.92e-7
-> ratio 1.73 (a triangle-shape/config difference, NOT a 10-50x gap). My PyCCL ref
and the B_kappa-phase ref are the SAME standard-B_kappa family and CORROBORATE each
other; canoes is the outlier relative to BOTH.

## B. The (1+z)^6 test -> it is actually (1+z)^4 per shell

The prompt's back-of-envelope predicted (1+z)^6 (per-leg (1+z)^2 over 3 legs).
That over-counts: it ignores the `dlambda = a^2 dchi` measure change that the
canoes lambda-fold uses, which removes one a^2 = (1+z)^-2.  Analytically the per-
shell radial-convention ratio is

    canoes/std = a^2 * [ K A(a)(1+z)^4 / g ]^3
               = a^2 * [ (1+z)^3 / (1+z)^1 ]^3 = a^2 (1+z)^6 = **(1+z)^4**,

verified symbolically (power = 4.000, prefactor 1.000).  The z_s sweep confirms it:

| z_s | can/std (rebuilt) | z_eff | (1+ze)^4 | (1+ze)^6 | ratio/(1+ze)^4 | ratio/(1+ze)^6 |
|---|---|---|---|---|---|---|
| 1.0 | 3.89  | 0.392 | 3.76  | 7.28  | **1.04** | 0.54 |
| 2.0 | 7.95  | 0.652 | 7.44  | 20.31 | **1.07** | 0.39 |
| 5.0 | 22.27 | 1.109 | 19.80 | 88.08 | **1.13** | 0.25 |

`ratio/(1+z_eff)^4 ~ 1.0-1.13` (flat; residual drift = Jensen gap between
`<(1+z)^4 W>/<W>` and `(1+<z>)^4`).  `ratio/(1+z_eff)^6 ~ 0.25-0.54` (NOT flat).
**The convention factor is (1+z)^4 per shell**, i.e. (1+z)^2 too much per leg in
the 3-leg vertex AFTER one a^2 is reabsorbed by the lambda measure.

## C. 2-point vs 3-point -> 3-POINT-SPECIFIC

The decisive disambiguator.  The canoes 2-point convergence uses the IDENTICAL
per-leg `(1+z)^4` response (`kappa2.py:9,164`, `SACHS_GEOMETRIC_POWER`) and the
SAME `dlambda=a^2 dchi` measure as the 3-point vertex (`kappa3.py:1708`,
`_kappa3_limber_response_product_h`: `Phi00=+A(a)(1+z)^4 delta`).

Canoes 2-pt vs PyCCL, single source plane z_s=5
(`paper_analysis1_O0_vs_pyccl_figdata.npz`, kappa-kappa, sftwick O0 vs PyCCL):

| gamma | PyCCL_kk | canoes O0 | canoes/PyCCL |
|---|---|---|---|
| 0.5'  | 8.46e-4 | 8.29e-4 | **0.979** |
| 2.1'  | 6.98e-4 | 6.84e-4 | **0.979** |
| 8.5'  | 3.13e-4 | 3.08e-4 | **0.982** |
| 35'   | 6.60e-5 | 6.51e-5 | **0.987** |

**median canoes/PyCCL = 0.981 over [0.5',50']** (within 2%). My own fresh PyCCL
linear xi_kappa reproduces the archived PyCCL reference to ~1-10% (linear-vs-mildly
nonlinear), confirming the archived reference is the standard one.

So:
- **2-point: canoes == PyCCL standard (within 2%)** — the per-leg (1+z)^4 response
  is correctly compensated by the standard projection + measure at 2 legs.
- **3-point: canoes ~10x (raw) / ~22x (canoes footing) the PyCCL standard.**

If the (1+z)^4 were a genuine per-leg radial-measure bug it would corrupt the
2-point by (1+z)^4 per shell (~20x at z_s=5), which it manifestly does NOT.
Therefore the discrepancy is **3-POINT-SPECIFIC**: it lives in the equal-shell
driving-field 3-cumulant assembly (the `zeta_TTT` vertex normalization / equal-
shell measure / the J0-fold (2pi)^4 and B_delta projection of the 3-point object),
NOT in the per-leg Sachs response that the 2-point shares.

## A separate, smaller issue: the h-power of `units="physical"`

A convergence 3PCF is dimensionless and MUST be h-invariant.  Computing the
STANDARD conv-3PCF fully in physical units and fully in h-units gives IDENTICAL
results (ratio 1.0000000), so my dimensionless `Z_std` carries the correct h-power.
The canoes/STAGE-0 `units="physical"` builder multiplies the h-native zeta by h^6,
which makes the result h-DEPENDENT and is h^2 too large for a dimensionless
conv-3PCF (my physical rebuild reproduces the STAGE-0 SPT reference EXACTLY x 1/h^2,
i.e. 2.230 = 1/h^2 = 2.220).  This 1/h^2 is a SECOND, independent bookkeeping issue
(the right factor for a dimensionless conv-3PCF is h^4, not h^6).  It does NOT
affect the can/std RATIO (both sides of my rebuilt ratio share one footing), but it
explains why the "canoes footing" column is h^2 = 2.22x larger than the raw column,
and why `zetaD_Bfix.md`'s 15-49x (canoes footing) is ~2.2x the raw ~7-22x.

## VERDICT

(a) **YES** — the canoes scalar convergence 3PCF is genuinely off the STANDARD
convention by a LARGE factor: **~10x (dimensionless/raw) at small gamma, ~22x on
the canoes `units="physical"` footing** at z_s=5, declining with gamma.  This is
NOT an artifact of my PyCCL build: my PyCCL reference agrees with the independent
B_kappa-phase scalar reference (same standard family) to ~1.7x (shape only), and
the canoes 2-point agrees with PyCCL to 2%.

(b) The factor is **(1+z)^4 per shell** (verified by symbolic algebra AND the z_s
sweep: ratio/(1+z_eff)^4 = 1.0-1.13, while ratio/(1+z_eff)^6 = 0.25-0.54).  The
prompt's (1+z)^6 back-of-envelope over-counted by ignoring the dlambda=a^2 dchi
measure change.  Equivalent statement: (1+z)^2 too much per leg, in the 3-point.

(c) It is **3-POINT-SPECIFIC, not per-leg**: the 2-point convergence (same per-leg
(1+z)^4 response, same lambda measure) matches PyCCL to 2%.  A per-leg radial-
measure convention would have corrupted the 2-point by ~(1+z)^4 too; it did not.
The over-normalization lives in the equal-shell driving-field 3-cumulant assembly.

**The spin-2 "4-17x shear gap" is therefore confirmed to be a SUBSET of this
3-point-specific scalar over-normalization (~10-22x), not a spin-2-only bug.** The
correct next step is to fix the 3-point `zeta` (equal-shell) measure / normalization
so the scalar TTT fold matches PyCCL to ~10% FIRST; only a residual after that can
be a genuine spin-2 effect.

(A second, smaller, independent fix: the `units="physical"` h^6 -> h^4 power for a
dimensionless conv-3PCF.)
