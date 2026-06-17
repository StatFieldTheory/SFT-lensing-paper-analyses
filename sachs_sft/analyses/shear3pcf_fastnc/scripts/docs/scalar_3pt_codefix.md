# Scalar 3-pt equal-shell code fix: revert spurious `(1+z)^4` radial Jacobian + keep `h^6 -> h^4` units

Date: 2026-06-09 (DIRECTION-CORRECTED supersede of the earlier same-day note)
canoes branch: `fix/spin2-zetaD-2_2_-2-norm` (pre-existing uncommitted
`_spin_aware_three_pt.py` LOW-kernel / spin-structure fix left UNTOUCHED; NOT
committed/pushed).
File touched: `src/canoes/sachs/kappa3.py` (only). Tests NOT modified (see end).

FIRST PRINCIPLES IS THE ARBITER. An earlier working-tree edit (same day) had
INTRODUCED a per-shell `(1+z)^4` radial factor into the equal-shell 3-pt table,
justified by an analogy to the 2-pt kernel. The h-invariance boundary check and
the validated 2-pt comparator overturn that direction: the correct equal-shell
`radial_measure="lambda"` factor is `ones`, and the `(1+z)^4` was an
over-normalization. This note records the corrected state.

Two items in the 3-pt equal-shell path:

1. **Radial measure (FIX 2 — REVERT to `ones`).**
   `_kappa3_radial_density_conversion`, the `radial_measure=="lambda"` branch, is
   restored to the committed-HEAD value `np.ones_like(d_lambda_d_chi)` (tag
   `lambda_three_leg_native_density`). The same-day working-tree edit had changed
   it to `d_lambda_d_chi ** -2 = (1+z)^4`, which MULTIPLIES the equal-shell zeta
   table by `(1+z)^4` per shell at every assembly site. First-principles reason
   the factor must be `ones`:

   - The equal-shell 3-pt fold (appendix.tex eqs. "appendix kappa3 cross shell"
     L672-680 and "appendix equal shell approx" L683-687) is performed ENTIRELY in
     affine lambda: each leg integrates `int dlambda` and the two collapse deltas
     `delta(lam'_1-lam'_2) delta(lam'_1-lam'_3)` are AFFINE-lambda deltas applied
     against the LAMBDA-density equal-shell cumulant `zeta_abc(gamma, lam'_1)`.
     No change of variable to chi occurs, so NO per-leg chi-Jacobian may be
     introduced. The appendix is internally self-consistent in affine lambda and
     is the exact 3-leg analogue of the 2-pt master eq. "xi kappa folded" (also in
     lambda). The chi-Jacobian `a^2` per leg appears in the paper ONLY in the
     OPTIONAL explicit chi-rewrite (cosmology.tex "cumulant 2 measure change" /
     "A Unified View"); it is not applied to this lambda-native table.
   - Boundary check against the validated 2-pt path: the 2-pt comparator
     `_convert_kappa2_sample_density_measure` (kappa2.py) is a STRICT NO-OP for
     `radial_measure` in `{chi, lambda}` (`del d_lambda_d_chi, radial_measure;
     return cl_pp, cl_ps, cl_ss`), and that no-Jacobian lambda treatment matches
     PyCCL to 0.2-0.4%. The native per-leg Sachs kernel `A(a) (1+z)^4 delta` is
     therefore the correct lambda-density on BOTH the 2-pt and 3-pt sides
     (the per-leg `(1+z)^4` response already bakes in the per-leg lambda-density,
     just as sigma2's `(1+z)^8` does in the 2-pt kernel). Multiplying the 3-pt
     table by a further `(1+z)^4` double-counts the measure on the two collapsed
     legs. The earlier edit's "kappa2.py declares per-leg-(1+z)^4 chi-density-
     native, validated to 2%" claim was a misread: kappa2 carries a metadata
     LABEL but applies NO Jacobian, so "validated" means "use the (1+z)^4 value
     as-is" -> mandates `ones` for the 3-pt analogue.

2. **`h^6 -> h^4` units (FIX 1 — KEPT).** At all five equal-shell zeta builders.
   The equal-shell zeta density is natively `(h/Mpc)^4`
   (`A(a)^3 * (P*P) / chi^4 = (h/Mpc)^4`), so the physical conversion is `h^4`,
   not `h^6`. The surplus `h^2` violated the h-invariance of the dimensionless
   conv-3PCF fold `Z = INT dlambda K^3 zeta`. This is a z-independent units
   rescaling, ORTHOGONAL to item 1, and the PyCCL-clean 2-pt anchor
   (chi_delta.py / kappa2.py) already uses `h^4`. The five builders must change
   together (LOW/HIGH/modulus/fused) to keep cross-builder consistency; the two
   sibling LOW/HIGH sites at HEAD L3151/L3344 were already `h^4`.

NOT touched (per spec): per-leg `(1+z)^4` Sachs field prefactor, `D^4` growth,
Poisson `A(a)^3`, scalar/J_0 paths, spin structure (incl. the pre-existing
`_spin_aware_three_pt.py` working-tree edit), the 2-pt `kappa2` path, and the
cross-shell `compute_kappa3_sigma3` radial weights.

Convention note: with `lambda_convention="project"/"manuscript"`,
`d_lambda_d_chi = a^2 = (1+z)^-2`, so the reverted edit's `d_lambda_d_chi ** -2`
equals `(1+z)^4`; the restored value `ones` introduces no per-shell z dependence.

Appendix status: **appendix_correct_code_bug**. The appendix formalism passes the
first-principles check (the equal-shell fold is correct in affine lambda with NO
per-leg chi-Jacobian). NO appendix correction is emitted, and `sections/*.tex`
is not edited.

---

## Diff (file:line before -> after), `src/canoes/sachs/kappa3.py`

### FIX 2 — radial Jacobian (single locus, `_kappa3_radial_density_conversion`)

`lambda`-branch return tuple (working-tree edit REVERTED back to HEAD; the
diff vs HEAD for this tuple is now EMPTY):

```
   (working-tree edit, now reverted)            ->   (restored HEAD value, == FIX)
-  d_lambda_d_chi ** -2,                              np.ones_like(d_lambda_d_chi),
-  "lambda_two_collapsed_leg_chi_jacobian",           "lambda_three_leg_native_density",
-  "(d chi/d lambda)^2 = (1+z)^4 on the two ...",      "none (base kappa3 is natively
                                                        lambda-density; cosmology.tex eq Phi00 scalar)",
```

The docstring block above the return was also restored: the spurious
"Equal-shell collapse Jacobian (2026-06-09)" paragraph that justified `(1+z)^4`
was replaced by an "Equal-shell collapse measure (2026-06-09, first-principles
arbiter)" paragraph stating the factor is `ones` (with the affine-lambda fold +
2-pt no-op boundary-check argument).

This single function is the consumer for ALL equal-shell builders: the
`radial_factor` it returns is applied in `compute_kappa3_zeta_table` (LOW scalar),
`compute_kappa3_sigma3_high` (HIGH scalar), `compute_kappa3_mod_zeta_equal_time`
(modulus LOW), `compute_kappa3_mod_sigma3_high` (modulus HIGH), and
`compute_kappa3_bare_bl_table` (fused). One revert -> all paths corrected.

### FIX 1 — `h^6 -> h^4` at the five equal-shell builders

`compute_kappa3_sigma3_high` (HIGH scalar, HEAD L3642-3643):
```
-        zeta = {name: arr * h_val**6 for name, arr in zeta_h.items()}
-        unit_scaling = "physical high-tail zeta = h^6 * h-native value"
+        zeta = {name: arr * h_val**4 for name, arr in zeta_h.items()}
+        unit_scaling = "physical high-tail zeta = h^4 * h-native value"
```

`compute_kappa3_mod_zeta_equal_time` (modulus LOW, HEAD L3906-3907):
```
-        bmod_h = bmod_h * h_val ** 6
-        dmod_h = dmod_h * h_val ** 6
+        bmod_h = bmod_h * h_val ** 4
+        dmod_h = dmod_h * h_val ** 4
```

`compute_kappa3_mod_sigma3_high` (modulus HIGH, HEAD L4050-4051):
```
-        bmod_h = bmod_h * h_val ** 6
-        dmod_h = dmod_h * h_val ** 6
+        bmod_h = bmod_h * h_val ** 4
+        dmod_h = dmod_h * h_val ** 4
```

`compute_kappa3_zeta_table` (LOW scalar, HEAD L4872):
```
-        h6 = h_val ** 6
+        h6 = h_val ** 4      # name kept; value now h^4 (applied to all 4 zeta channels)
```

`compute_kappa3_bare_bl_table` (fused scalar, HEAD L5105):
```
-        post = post * (float(cosmo.h) ** 6)
+        post = post * (float(cosmo.h) ** 4)
```

Plus two doc/comment strings updated for consistency: the `units` docstring of
`compute_kappa3_zeta_table` ("multiply ζ by h^6" -> "h^4") and the post-processing
comment in `compute_kappa3_bare_bl_table` ("then h^6" -> "then h^4").

Net vs HEAD: 21 insertions, 12 deletions, 1 file (the FIX-2 return tuple matches
HEAD bit-for-bit; the remaining diff is the FIX-2 docstring rewrite + the FIX-1
`h^4` sites).

---

## Self-check (i) — per-shell scalar 3-pt response ratio S_can/S_std

Probe (project convention), extreme z per boundary-validation: z in {0.3, 0.5,
1, 2, 5} (low z_s ~0.3 AND high z_s ~5, the (1+z)^4 error being monotone in z).
`d_lambda_d_chi = a^2 = (1+z)^-2`.

`_kappa3_radial_density_conversion(d_lambda_d_chi, "lambda")` now returns
`ones` (tag `lambda_three_leg_native_density`).

S_can/S_std (overall sign -1 from kappa = -INT Delta theta):

```
                  z=0.3      z=0.5      z=1        z=2        z=5
radial_factor:    1          1          1          1          1
S_can/S_std:     -1         -1         -1         -1         -1          (target -1)
BUGGY edit:      -2.8561    -5.0625    -16        -81        -1296       (= -(1+z)^4)
AFTER all == -1: True
```

The reverted edit would have left the per-shell scalar 3-pt response a factor
`(1+z)^4` ABOVE the standard tree (e.g. -1296 at z_s=5); after the revert it is
exactly `-1` (i.e. `-S_std`) at every redshift, including both extremes. PASS.

## Self-check (ii) — h-invariance of the dimensionless conv-3PCF fold

Dimensionless fold `Z = INT dlambda K^3 zeta` must be h-invariant. Unit
bookkeeping: zeta natively `(h/Mpc)^4` and chi/lambda divided by h in physical
units; the code's physical scaling makes `Z_phys/Z_h = h^(N-4)` where `N` is the
exponent on the zeta multiply. h = 0.6711 (canoes default; `h^2 = 0.45038`).

```
  N (zeta multiply)   physical/h fold ratio   verdict
  h^6 (DEPLOYED-buggy) h^2 = 0.45038           h-NON-invariant (prompt-measured)
  h^4 (FIXED)          h^0 = 1.000000          h-invariant -> RESTORED
```

`Z_physical / Z_hnative = 1.000000` after the `h^6 -> h^4` change (was `0.45038
= h^2` exactly). PASS.

Both self-checks pass; no STOP condition triggered.

Module imports cleanly (`import canoes.sachs.kappa3` OK); the `radial_measure=="chi"`
branch is unchanged.

---

## Test fallout (NOT modified — surfaced for the user)

The kappa3 zeta-table tests pin the committed-HEAD convention. Because FIX 2 is a
REVERT to that HEAD value, the radial-Jacobian tests that pin
`radial_kernel_convention == "lambda_three_leg_native_density"` and the
lambda-branch no-Jacobian behaviour are CONSISTENT again with the code (the
same-day edit had broken them; the revert restores agreement). The only tests
that may now diverge are those pinning the physical units scaling at `h^6` (FIX 1
moves them to `h^4`):

- `test_kappa3_units_conversion_all_channels` (if it asserts physical = `h^6`) —
  now expects `h^4`.

Recommended follow-up (user decision): if any test still pins physical = `h^6`,
update it to `h^4` so it regression-guards the fixed (h-invariant) behaviour.
Tests were left unmodified per spec.
