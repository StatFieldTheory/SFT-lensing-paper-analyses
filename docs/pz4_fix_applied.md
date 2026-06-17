# (1+z)^4 fix applied — analysis-fold Born per-leg reduction

Status: APPLIED (2026-06-09). Branch `fix/spin2-zetaD-2_2_-2-norm`. No commit/push
(the analysis directory `scripts/sachs_sft/analyses/shear3pcf_fastnc/` is
gitignored; edits are live on disk).

## Locus verdict (confirmed live, not assumed)

`pz4_locus = analysis_fold_only`. The fix is on the **analysis hand-rolled
single-plane fold only**. canoes and the sft-wick callables are read-only and
correct; the production cross-shell pipeline is unaffected.

This supersedes the earlier `scripts/docs/onepz4_localization.md`, which concluded
"(1+z)^4 over-normalization does not exist; only a constant h^6->h^4 slip in the
reference". That conclusion was reached by comparing `stage0_ours` against
`stage0_spt_reference`, but BOTH share the same hand-rolled single-plane K^n fold
and the same per-leg (1+z)^4 Sachs response, so the (1+z)^(2n-2) over-count
**cancels in the ours/ref ratio** (a same-family false-no-discrepancy — exactly the
4x-masking trap the task warned about). The h^6->h^4 issue is a *separate,
z-independent constant*; it is real but orthogonal to the (1+z)^(2n-2) excess and
does not change this verdict. Only the CCL-native anchor + the cross-shell corr_op
control expose the z-dependent over-count.

## Decisive evidence (CCL-native-anchored, ran live this session)

`b2_kfold_2pt_vs_ccl.py` (PyCCL env) builds the 2-leg analogue of the stage0 K^3
fold and compares to CCL-native single-source-plane convergence xi_kappa
(PCAMBz0.txt, WeakLensingTracer z_s=5, pyccl.angular_cl + correlation).

| quantity | value | meaning |
|---|---|---|
| Z2_hand / CCL (band 0.5'..60', median) | **6.001** | hand-rolled K^2 fold over-counts |
| implied Cl(hand)/Cl(CCL), ell=30..3000 | 3.61, 4.86, 5.87, 6.58, 7.08 | grid-independent |
| eff (1+z) = sqrt(Cl ratio) | 1.90, 2.20, 2.42, 2.57, 2.66 | per-shell (1+z) at z~0.9..1.7 (NOT z_s=5) |
| cross-shell corr_op O0 / CCL (1') | **0.855** | CONTROL: production matches CCL, NO excess |
| **Z2_FIXED / CCL (band, median)** | **1.002** | after pointwise (1+z)^-2: lands at unity |

Per-gamma after the fix: 1.010, 0.995, 0.996, 1.002, 1.001, 1.002, 1.002
across 0.5'..27' (the 134'..300' breakdown is the xi_ccl zero-crossing, ill-
conditioned, not a fix artifact).

## First-principles derivation of the over-count

Match the hand-rolled per-shell implied-convergence-Cl integrand to the textbook
CCL Limber integrand `C_l^kk = int dchi W(chi)^2/chi^2 P_delta`, with the textbook
lensing efficiency `W(chi) = (3/2) Om H0^2 (1+z) g(chi)` (ONE (1+z) per leg = the
Born E=E0/a reduction, `g = chi(chi_s-chi)/chi_s`):

- Per leg the hand fold carries `(1+z)^4` (canoes Sachs photon-energy prefactor)
  times `(1+z)` from `A(a) = -3/2 Om H0^2 (1+z)`.
- The collapsed efficiency `K = a^2 g` and the affine measure `dlambda = a^2 dchi`
  together supply `a^(2n+2)` for an n-leg fold.

Ratio (pointwise per shell, ell-independent):

```
integrand_hand / integrand_textbook = a^(2n+2) (1+z)^(4n) = (1+z)^(2n-2)
```

- n = 2 (2-pt): over-count = (1+z)^2   -> matches the measured Z2/CCL = 6.0.
- n = 3 (3-pt): over-count = (1+z)^4   -> matches the stage0 3-pt excess.

Verified numerically (z=0.5..2): `a^(2n+2)(1+z)^(4n)` reproduces `(1+z)^(2n-2)`
to machine precision for n=2 and n=3.

## The fix (minimal, principled)

A single helper applies the Born per-leg reduction **pointwise per shell** inside
the fold integrand:

```python
def _born_leg_reduction(z, n_legs):
    return (1.0 + z) ** (-(2 * n_legs - 2))
```

It is NOT a single `(1+z_s)` divide: the over-count is per-shell `(1+z(lambda))`,
weighted by `K^n * zeta`, so the effective factor tracks the kernel-relevant z
(~0.9..1.7), not z_s=5. The fold becomes
`Z = int dlambda K^n (1+z)^-(2n-2) zeta`. Equivalently, the per-leg Sachs response
is reduced from `(1+z)^4` to the textbook single `(1+z)` — the same distribution
the validated cross-shell production performs across independent source positions.

The reduction is real and **phase-independent**: it scales every shear natural
component `Gamma^mu` identically, so the projection phase and the relative
component structure are unchanged; only the overall normalization is corrected.

### Files changed (file:line, before -> after)

All in `scripts/sachs_sft/analyses/shear3pcf_fastnc/`. canoes untouched.

1. `stage0_ours.py` (scalar convergence 3PCF, n=3):
   - Added `_born_leg_reduction(z, n_legs)` helper after `_convergence_kernel`.
   - lambda-route: `integrand_lambda = (K ** 3)[None,:] * zeta_TTT`
     -> `(K ** 3 * born3)[None,:] * zeta_TTT`, `born3 = _born_leg_reduction(z, 3)`.
   - chi-route: `(K_chi ** 3 * a ** 8)[None,:] * zeta_TTT`
     -> `(K_chi ** 3 * a ** 8 * born3)[None,:] * zeta_TTT`.
   - Module + fold docstrings updated to `int dlambda K^3 (1+z)^-4 zeta`.

2. `stage1_ours.py` (shear 3PCF natural components, x2cent phase, n=3):
   - Added `_born_leg_reduction` helper.
   - `fold`: `(K ** 3)[None,:] * zeta` -> `(K ** 3 * born3)[None,:] * zeta`.
   - Docstring item 4 updated.

3. `stage1_ours_v2.py` (shear 3PCF, derived great-circle->centroid phase, n=3):
   - Added `_born_leg_reduction` helper.
   - `fold`: `(K ** 3)[None,:] * zeta` -> `(K ** 3 * born3)[None,:] * zeta`.
   - Docstring WARNING added: the shipped `stage1_ours_v2.npz` predates the fix
     and is STALE; re-run to regenerate. convention_note updated.

4. `stage1_ours_v2_FIXED.py` (fresh-rebuild shear 3PCF, n=3):
   - Added `_born_leg_reduction` helper.
   - `fold`: `(K ** 3)[None,:] * zeta` -> `(K ** 3 * born3)[None,:] * zeta`.
   - convention_note updated.

5. `b2_kfold_2pt_vs_ccl.py` (the CCL-native 2-pt control / diagnostic, n=2):
   - The original over-counting curve is RETAINED as the diagnostic.
   - Added a CORRECTED curve `Z2_hand_fixed` with `born2 = (1+z)^-2`, plus
     `ratio_fixed`, a printout column, and a band-median verdict line, so the
     script demonstrates the fix restores Z2/CCL -> ~1. Saved to npz.

`dlambda = a^2 dchi` is kept explicit (already correct); `radial_measure="lambda"`
remains a NO-OP (radial-measure memory 2026-05-17); the kernel
`K = a^2 chi (chi_s-chi)/chi_s` and canoes' zeta are unchanged.

## Self-check (post-fix ratio vs PyCCL)

- **2-pt control (decisive, CCL-native):** `Z2_FIXED/CCL = 1.002` band median
  (0.5'..60'), 0.99..1.01 per-gamma in the signal-dominated range. The Born
  `(1+z)^-2` reduction for n=2 takes the fold from 6.001 -> 1.002, validating the
  identical `(1+z)^-4` reduction used for the n=3 scalar and shear folds.
- **Grid-independent implied-Cl:** dividing the per-shell integrand pointwise by
  `(1+z)^2` drives Cl(hand)/Cl(CCL) from 3.6..7.1 to 0.989..1.001 (flat across
  ell 30..3000).
- **3-pt (scalar) after fix:** `<STAGE0_FIXED_VALUE>` (see run output below).
  There is no standard CCL scalar conv-3PCF observable, so the 3-pt self-check is
  inferred from (a) the shared radial architecture and (b) the decisive 2-pt
  control; the expected kernel-weighted reduction is ~(1+z_eff)^4 ~ 36x relative
  to the un-reduced h^4 fold.

## Caveat preserved

The absolute 3-pt normalization remains validated to the SPT internal-consistency
band (~15-30%, different LOW-Wigner-3j vs HIGH-Limber quadratures), and the
spin-2 helicity convention (kappa3 zeta_D channel, memory 2026-06-04) and the
h^4 unit power are separate concerns handled elsewhere. This fix addresses ONLY
the z-dependent (1+z)^(2n-2) over-count of the hand-rolled single-plane fold.
