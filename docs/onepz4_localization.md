# Localization of the canoes 3-pt "(1+z)^4 over-normalization"

Status: RESOLVED (2026-06-09). Synthesis of CHECK A (reference trustworthiness)
and CHECK B (zeta vs fold). Both probes are READ-ONLY; outputs in
`scripts/sachs_sft/analyses/shear3pcf_fastnc/outputs/checkA_reference_2pt_vs_ccl.npz`
and `.../checkB_zeta_vs_fold.npz`.

## TL;DR

The presupposed "(1+z)^4 over-normalization" **does not exist**. There is no
redshift-dependent error anywhere in the 3-pt chain. The only real discrepancy
is a **constant, redshift-independent factor h^-2 = 2.2204** (h = 0.6711), and
its single locus is a **physical-units conversion power in the REFERENCE**
(`stage0_spt_reference.py`, line 208: `* h**6` should be `* h**4`), plus the
stale `stage0_ours_COMMITTED_h6.npz` that shared the same h^6 convention.

> Single most-supported locus: **stage0_spt_reference.py line 208 (reference-side
> h-power), NOT canoes zeta_TTT and NOT the stage0_ours analysis fold.**

This is a clean resolution of the decision logic, but with an important twist:
the decision tree's three candidate loci were
{PyCCL-3pt-reference, canoes zeta_TTT, stage0_ours fold}, and the true locus is a
**fourth thing** -- the reference's *unit-conversion power*, which is a property
of the reference but is NOT a radial-measure / W-kernel / Jacobian bug (the
thing CHECK A was built to exonerate). CHECK A correctly clears the reference's
*radial machinery*; the h-power bug is orthogonal to that and is caught by
CHECK B's dimensional analysis.

## What each check establishes

### CHECK A -- the reference's radial machinery is CCL-correct (TRUSTWORTHY)

The reference's own 2-point analogue (same g(chi), same CCL distances, same CCL
growth D(a), same PCAMBz0 P(k) via the identical `pk_phys_z0` closure, same
single source plane z_s=5), reduced to the Limber projection
`C_l^kk = INT dchi g^2/chi^2 D^2 P_delta((l+1/2)/chi)`, reproduces CCL-native
`WeakLensingTracer.angular_cl` to <1% across the lensing band.

Verified ratios (densest reference grid n_lambda=600), C_l^kk ref/CCL:

| l    | 10    | 18    | 48    | 125   | 261   | 675   | 1271  | 2391  | 5000  |
|------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| ratio| 1.031 | 1.013 | 1.002 | 1.002 | 1.000 | 1.000 | 1.004 | 1.006 | 1.009 |

- Band l in [100,2000]: median ref/CCL = **1.0007**, mean 1.0013, min..max 0.9991..1.0048.
- Reference's own n_lambda=40 vs 600 radial drift in-band: median 0.18%.
- Real-space xi_kappa(theta) ref/CCL: median over 1'..100' = **1.0022**
  (the 0.82 at theta=300' is the xi zero-crossing, ill-conditioned, not an error).

**No (1+z)^n trend anywhere** -- the reference 2pt is flat-unity, not scaled by
any power of (1+z). Therefore the reference's W kernel, chi-measure, and Limber
reduction carry NO measure error, and the project lead's Jacobian-matching
warning is satisfied on the reference side (it integrates g(chi) in CCL-native
chi coords; no lambda-vs-chi mismatch).

### CHECK B -- the only discrepancy is a constant h^-2 in the reference's unit conversion

The decisive B1 probe compared canoes `zeta_TTT` against the independent
first-principles SPT tree+Limber reference **in h-units** (no physical
conversion, isolating the physics), backed by dimensional analysis. Findings,
all re-verified against the npz and source on 2026-06-09:

- **Apples-to-apples (h-units) canoes zeta_TTT / reference = 0.810** at
  gamma=1', mid shell z=2.625 (~19% -- the expected SPT internal-consistency band
  for two independent quadratures of the same tree+Limber family; magnitude and
  sign agree). This is the physics check, and it PASSES.
- **canoes physical-conversion factor = physical/h-unit = 0.202838 = h^4 EXACTLY**
  (h=0.6711). Confirmed: zeta_h = -5.316e-13; zeta_h * h^4 = -1.078e-13 = the
  canoes "physical" value. This is dimensionally correct: the equal-shell scalar
  3-cumulant is natively (h/Mpc)^4 (= A(a)^3 P_delta P_delta / chi^4 after the
  dimensionless Poisson factor A/k^2 and dimensionless ell-integral), so the
  correct conversion is *h^4. canoes `units="physical"` does this correctly
  (kappa3.py line 4216).
- **The reference applies *h^6 instead** (stage0_spt_reference.py line 208),
  making it h^2 too small. The stale `stage0_ours_COMMITTED_h6.npz` shared this
  h^6 convention.
- **Inter-version constant: current(canoes h^4) / committed(h^6) = 2.2203708769849917
  = h^-2 EXACTLY** (ln(ratio)/ln(1/h) = 2.0000000 to machine precision,
  **z-INDEPENDENT across all 40 shells** -> definitively excludes any
  (1+z)-power explanation).

Folded Z_kappa(gamma=1', physical):

| quantity            | value        | note                                  |
|---------------------|--------------|---------------------------------------|
| current (canoes h^4)| -5.074e-05   | correct                               |
| committed (h^6)     | -2.285e-05   | h^2 too small (shares ref's bug)      |
| reference (h^6)     | -2.701e-05   | h^2 too small (the bug, line 208)     |
| current / ref       | **1.878**    | apparent "mismatch" = ref not fixed   |
| committed / ref     | **0.846**    | agree because both carry h^6          |

The 1.878 "over-normalization" that motivated this investigation is therefore an
**artifact of comparing the correctly-h^4 current `stage0_ours` against the
NOT-YET-FIXED h^6 reference**. Once the reference is corrected to h^4, ours/ref
returns to ~0.85, i.e. the genuine ~15-19% SPT internal-consistency level.

## Reconciling A and B (no inconsistency)

A and B are fully consistent; they exonerate different facets:

- CHECK A clears the reference's **radial measure / W-kernel / Limber reduction**
  (the *shape* and *chi-integration* machinery). True, and it is what A tests.
- CHECK B clears canoes **zeta_TTT** (physics, 0.81 in h-units) AND the
  **stage0_ours K^3/lambda fold** (matches the validated 2-pt convention;
  radial_measure='lambda' is a verified NO-OP, per-leg (1+z)^4 applied
  identically on both sides, lambda-route vs chi-route agree to <1e-2).
- The bug CHECK B finds (h^6 vs h^4) is a scalar overall normalization that
  **does not touch the radial integrand shape** -- so CHECK A's 2pt-vs-CCL test,
  which uses the reference's radial machinery but its OWN 2pt normalization path
  (not the 3pt h-power line), is blind to it by construction and correctly
  reports ~unity. No contradiction.

The decision-tree branch that fires is therefore **NOT** the clean "A trustworthy
+ B zeta correct => bug is in the stage0_ours fold" branch. Both zeta AND the
fold are correct. The (1+z)^4 framing was a false premise: the real defect is a
constant h-power in the reference, which the tree did not enumerate as a distinct
locus.

## Principled fix (single locus, reference-side)

The fix is **reference-side / analysis-side bookkeeping only -- NOT canoes-side,
NOT a change to zeta_TTT, NOT a change to the K^3/lambda fold geometry.**

1. `scripts/sachs_sft/analyses/shear3pcf_fastnc/stage0_spt_reference.py`
   - line 208: `zeta_TTT = zeta_h * (h ** 6)` -> `zeta_h * (h ** 4)`
   - docstring line 56: "convert ... by h^6 (zeta carries h^6)" -> "by h^4 (zeta carries h^4)".
   - Rationale: the equal-shell scalar 3-cumulant is natively (h/Mpc)^4; the
     correct h-unit -> physical conversion is *h^4, matching canoes
     `units="physical"` (which is already correct).
2. `scripts/sachs_sft/analyses/shear3pcf_fastnc/stage0_ours.py`
   - docstring line 38: "units='physical': zeta carries h^6" -> "h^4".
   - (Code is already correct -- it consumes canoes `units='physical'` directly
     with no extra h factor. The 13:16->13:21 change from h^6 to h^4 today was a
     correct fix. Only the docstring is stale.)
3. Regenerate/retire the stale `stage0_ours_COMMITTED_h6.npz` (it carries the
   wrong h^6 convention and should be archived, not consumed).

After this fix, the ours/ref fold ratio returns to ~0.85, which is the expected
SPT internal-consistency level (different ell-grids: canoes LOW Wigner-3j + HIGH
Limber vs reference 220-node log-ell + 96 phi), NOT a normalization error.

No change to canoes. No change to the per-leg (1+z)^4 Sachs prefactor
(`born_limber.py` T_phi00=(1+z)^4 L2/2chi^2). No change to the lambda-measure
NO-OP (radial-measure memory 2026-05-17 stands).

## Confidence

**HIGH (~0.9).** The h-power is pinned to machine precision (h^-2 exact,
power = 2.0000000, z-independent across all 40 shells), the provenance arithmetic
reproduces exactly (zeta_h*h^4 = physical value; zeta_h*h^6 = committed value),
and the two checks are mutually consistent with the validated 2-pt anchor
(CCL to 0.1-0.2% in-band, A; canoes/CCL 2pt to ~2%, memory). The dimensional
argument (native (h/Mpc)^4 => *h^4) is independent confirmation of the numeric
factor.

## Remaining ambiguity

1. **The ~19% h-unit residual (canoes zeta_TTT / ref = 0.81)** is attributed to
   SPT internal-consistency (different quadratures of the same family), not a
   bug. This is plausible and within the documented ~30% expectation band, but it
   is NOT pinned to a specific cause -- it could hide a small (few-%) shape error
   in either quadrature. It is not a (1+z)-power (CHECK B confirms z-independence
   of the version-to-version factor), so it does not affect the locus verdict,
   but the absolute 3pt normalization is only validated to ~20%, not to ~1%.
2. **CHECK A validates only the scalar/radial+P(k)+Limber machinery.** The 3pt
   reference additionally uses an F2-SPT B_delta and a J0 angular fold not
   exercised by the 2pt analogue. A pure bispectrum-shape or J0-quadrature error
   would not be caught by A. CHECK B mitigates this: the ell-grid/J0 machinery is
   byte-identical between Z_std and Z_can in the probe, so the can/std RATIO is
   insensitive to it -- but an absolute B_delta/J0 error common to both would
   survive. This again does not look like (1+z)^4 and does not move the locus.
3. The fix is a **recommendation** -- both probes were strictly READ-ONLY; the
   h^6->h^4 edits to `stage0_spt_reference.py` (line 208, docstring 56) and the
   docstring edit to `stage0_ours.py` (line 38) have NOT been applied, and the
   stale `stage0_ours_COMMITTED_h6.npz` has NOT been retired.

## Inputs / provenance

- CHECK A output: `scripts/sachs_sft/analyses/shear3pcf_fastnc/outputs/checkA_reference_2pt_vs_ccl.npz`
- CHECK B output: `scripts/sachs_sft/analyses/shear3pcf_fastnc/outputs/checkB_zeta_vs_fold.npz`
- Reference builder: `scripts/sachs_sft/analyses/shear3pcf_fastnc/stage0_spt_reference.py` (line 208 = locus)
- Ours fold: `scripts/sachs_sft/analyses/shear3pcf_fastnc/stage0_ours.py` (correct; docstring stale)
- canoes physical conversion: `canoes/sachs/kappa3.py` line 4216 (already correct, *h^4)
- Stale artifact to retire: `.../outputs/stage0_ours_COMMITTED_h6.npz`
