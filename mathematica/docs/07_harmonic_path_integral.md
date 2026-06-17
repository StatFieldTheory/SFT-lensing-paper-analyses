# 07 — `harmonic_path_integral_verification.wl`

> Numerical verification of the three identities underlying the draft
> subsection
> [`sections/harmonic_path_integral.tex`](../../sections/harmonic_path_integral.tex)
> ("Harmonic-Space Representation of the Sachs path integral").

## 1. Purpose and paper anchor

Script: [`scripts/harmonic_path_integral_verification.wl`](../harmonic_path_integral_verification.wl).

Paper anchors (labels the draft subsection introduces):

- `eq: harmonic fields` — scalar-harmonic transform of $\Dsachs$ and $\RespField$ at fixed $\lambda$.
- `eq: harmonic resp op` — $(LM, \lambda)$ response propagator (block-diagonal).
- `eq: harmonic corr op`, `eq: harmonic corr op reduced` — correlation propagator diagonalised in $(LM)$ with reduced spectrum $\tilde{\CorrOp}^{ij}_{L}(\lambda_{1}, \lambda_{2})$.
- `eq: src cumulant harmonic` — $Y_{LM}$ form of the isotropic two-cumulant.
- `eq: harmonic vertex projection`, `eq: harmonic vertex gaunt` — local-vertex reduction to a triple (or higher) $Y^{*}$ integral.

This script checks the three algebraic identities the draft relies on,
against direct spherical-sphere numerical integration.

## 2. First principles

- **Addition theorem** for Legendre polynomials:
  $P_{L}(\direc_1 \cdot \direc_2) = \frac{4\pi}{2L+1}\sum_{M} Y^{*}_{L M}(\direc_2)\,Y_{L M}(\direc_1)$.
- **Orthogonality** of spherical harmonics:
  $\int Y^{*}_{L M}\,Y_{L' M'}\,\mathrm d\Omega = \delta_{L L'}\,\delta_{M M'}$.
- **Wigner-3j form** of the conjugated triple-$Y$ integral:
  $\int Y^{*}_{L_1 M_1} Y^{*}_{L_2 M_2} Y^{*}_{L_3 M_3}\,\mathrm d\Omega = (-1)^{M_1+M_2+M_3}\sqrt{\frac{(2L_1+1)(2L_2+1)(2L_3+1)}{4\pi}}\,\binom{L_1\;L_2\;L_3}{0\;0\;0}\binom{L_1\;L_2\;L_3}{-M_1\;-M_2\;-M_3}$,
  nonzero only when $M_1+M_2+M_3 = 0$.

## 3. Tests

| Test | Claim | Tolerance | Cases |
|---|---|---|---|
| V1 | For $\kappa(\direc_1,\direc_2) = \sum_{L}\frac{2L+1}{4\pi}c_{L}P_{L}$, $\int\int Y^{*}_{L_1 M_1}\,Y_{L_2 M_2}\,\kappa\,\mathrm d\Omega_{1}\,\mathrm d\Omega_{2} = \delta_{L_1 L_2}\delta_{M_1 M_2}\,c_{L_1}$ | $10^{-4}$ | 11 tuples, $L\in[0,3]$, diagonal and off-diagonal |
| V2 | Triple-$Y^{*}$ integral = Wigner-3j decomposition; selection rule $M_1+M_2+M_3 = 0$ | $10^{-6}$ | 6 tuples, triangle-allowed |
| V3 | Separable $C_{L}(\chi,\chi') = f(\chi)f(\chi')g(L)$ reduces at Born to $\aps_{L}(\lambda_1,\lambda_2)$ by pure radial relabelling $\chi \leftrightarrow \chi(\lambda)$ | $10^{-12}$ | 5 tuples, $L\in\{0,1,2,3,5\}$ |

V1 drives the `eq: src cumulant harmonic` and `eq: harmonic corr op`
results. V2 underwrites `eq: harmonic vertex gaunt` (and reuses the
same kernel as the separate
[`jacobian_remap_verification.wl`](../jacobian_remap_verification.wl),
with which it shares the `Gaunt3j` helper). V3 spot-checks the
bridge sentence in the draft's final paragraph (Born-level
$(\ell m, \chi) \leftrightarrow (\ell m, \lambda)$) by construction:
the two parameterisations of a separable spectrum differ only by a
monotonic radial substitution with no Jacobian.

## 4. Reproduction

```bash
wolframscript -file scripts/harmonic_path_integral_verification.wl
```

Run time: ~60 s on a laptop. V1 dominates because of the 4D sphere
integral; V2 and V3 are near-instant.

## 5. Notes on V1's numerical integration

The 4D `NIntegrate` occasionally emits `NIntegrate::slwcon` warnings
when the exact result is zero (off-diagonal $L_1 \neq L_2$ cases),
because the integrator cannot confirm a zero to $10^{-15}$ precision
against a double-sphere background. The returned value is still
$\sim 10^{-10}$, well inside the tolerance, and the test passes. If
tighter tolerances are needed, increase `WorkingPrecision` and
`AccuracyGoal` in the `NIntegrate` call.

## 6. Inputs / outputs

**Inputs:** hard-coded Legendre coefficients $c_{L}$, test tuples, and
a toy Born cosmology ($a(\tau) = \tau$, $\tau_{0} = 1$, $E_{0} = 1$,
$\chi(\lambda) = \tau_{0} - \sqrt{\tau_{0}^{2} - 2 E_{0}\lambda}$).

**TeX markers:** none. This is a pure-numerical verification feeding
no displayed equation; its role is to confirm the algebraic identities
independently of the draft's derivation.

## 7. Upstream / downstream

- **Upstream:** none beyond Mathematica's built-ins.
- **Downstream:** supports the draft subsection
  [`sections/harmonic_path_integral.tex`](../../sections/harmonic_path_integral.tex).
  When/if that draft is merged into
  [`sections/path_int.tex`](../../sections/path_int.tex), this script
  stays in place and remains the algebraic check for the merged
  subsection.

## 8. Related scripts

- [`06_jacobian_remap_verification.md`](06_jacobian_remap_verification.md)
  — verifies the gradient-gradient identity used in Appendix
  `append: jacobian remap`; the helper `Gaunt3j` is shared (duplicated
  in-script for self-containment).
