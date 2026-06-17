# 06 — `jacobian_remap_verification.wl`

> Numerical verification of the three harmonic-space identities
> underlying Appendix "First-Order Jacobian-Kernel Remapping of
> Angular Statistics".

## 1. Purpose and paper anchor

Script: [`scripts/jacobian_remap_verification.wl`](../jacobian_remap_verification.wl).
Paper anchor: `sections/appendix.tex` `\label{append: jacobian remap}` and
the equations `\label{eq: appendix dCl rad}`,
`\label{eq: appendix dCl ang}`.

This script verifies by direct numerical integration that the key
harmonic-space identity
$$
F[l,m,L,M,l',m'] \;=\; \tfrac{1}{2}\bigl[L(L+1) + l'(l'+1) - l(l+1)\bigr]\,G[l,m,L,M,l',m']
$$
holds between the **gradient-gradient kernel**
$$
F[l,m,L,M,l',m'] \;\equiv\; \int_{S^2} Y^{*}_{l,m}\,(\nabla^A Y_{L,M})\,(\nabla_A Y_{l',m'})\,\mathrm d\Omega
$$
and the **conjugated Gaunt kernel**
$$
G[l,m,L,M,l',m'] \;\equiv\; \int_{S^2} Y^{*}_{l,m}\,Y_{L,M}\,Y_{l',m'}\,\mathrm d\Omega .
$$
This identity is the central engine of `eq: appendix dCl ang`:
the transverse-Jacobian convolution kernel in the first-order power-spectrum
correction reduces to the Gaunt kernel times an algebraic factor of the
Laplacian eigenvalues, matching the standard CMB-lensing literature.

## 2. First principles

- **Orthogonality of spherical harmonics** on the unit 2-sphere with
  metric $\mathrm ds^2 = \mathrm d\theta^2 + \sin^2\theta\,\mathrm d\phi^2$.
- **Wigner-3j form of the Gaunt integral** (Condon–Shortley
  convention): $G = (-1)^m \sqrt{(2l+1)(2L+1)(2l'+1)/(4\pi)}\,
    \bigl(\substack{l\;L\;l' \\ 0\;0\;0}\bigr)\,
    \bigl(\substack{l\;L\;l' \\ -m\;M\;m'}\bigr)$.
- **Gradient-gradient identity** follows from
  $\nabla^2(fg) = f\,\nabla^2 g + g\,\nabla^2 f + 2\,(\nabla^A f)(\nabla_A g)$
  combined with $\nabla^2 Y_{l,m} = -l(l+1)\,Y_{l,m}$; integrating by
  parts on the (boundary-less) sphere converts the kernel to the Gaunt
  kernel times the Casimir difference
  $\tfrac{1}{2}[L(L+1)+l'(l'+1)-l(l+1)]$.

## 3. Tests performed

| Test | Claim | Tolerance | Cases |
|---|---|---|---|
| 1 | $\int Y^*_{l,m}\,Y_{l',m'}\,\mathrm d\Omega = \delta_{ll'}\delta_{mm'}$ | $10^{-8}$ | 6 tuples |
| 2 | $G$ numerical = Wigner-3j decomposition | $10^{-6}$ | 6 tuples (with $m=M+m'$) |
| 3 | $F = \tfrac{1}{2}[L(L+1)+l'(l'+1)-l(l+1)]\,G$ | $10^{-4}$ | 7 tuples |
| 4 | $L=0$ or $l'=0 \Rightarrow F = 0$ (since $\nabla Y_{0,0}=0$) | $10^{-6}$ | 6 tuples |

Selection-rule warning: Mathematica's `NIntegrate` emits `slwcon` /
`eincr` diagnostics when the exact integral is zero (selection-rule
forbidden cases); the returned value is still $\sim 10^{-17}$, well
inside the tolerance.

## 4. Reproduction

```bash
wolframscript -file scripts/jacobian_remap_verification.wl
```

Runs in a few tens of seconds on a laptop (dominated by the
gradient-gradient integrals, which have more complex integrands).

## 5. Output

All 25 checks print `PASS` with residuals at or below $10^{-10}$
for Tests 1 and 3, and machine precision for Test 2 (where both
sides go through the same `ThreeJSymbol` pipeline up to the direct
numerical integration).

## 6. Inputs / outputs

**Inputs:** hard-coded $(l, m, L, M, l', m')$ test tuples; no external
files.

**TeX markers:** none. This is a pure-numerical verification that
feeds no displayed equation — its role is to confirm, independently
of the paper's algebraic derivation, that the identity underlying
Appendix D's power-spectrum formulas is correct.

## 7. Upstream / downstream

- **Upstream:** none. The script uses only Mathematica's built-in
  `SphericalHarmonicY` and `ThreeJSymbol`.
- **Downstream:** supports the appendix `\label{append: jacobian remap}`.
  If the Jacobian remapping kernels are ever re-derived or generalised
  (e.g., to include vector/tensor-sector $B$-mode deflections), this
  script provides the numerical baseline against which the new
  kernels should be checked.

## 8. Extensions (future work)

- **Spin-s generalisation.** The scalar-sector deflection is a pure
  gradient; vector and tensor sectors add a curl component and a
  spin-2 deflection. Analogous identities involve spin-weighted
  harmonics ${}_s Y_{\ell m}$ and can be checked by adapting
  `GradGradNum` to use spin-weighted covariant derivatives.
- **Full $\Delta C_\ell^{\rm ang}$ check.** Squaring
  `eq: appendix alm remap` and summing $m$ gives the standard
  Lewis–Challinor lensing-correction formula. A Monte-Carlo
  cross-check (drawing Gaussian $a_{\ell m}$, $\psi_{LM}$, $\delta\chi_{LM}$
  and computing the empirical power spectrum) would verify the full
  `eq: appendix dCl ang` end-to-end — not implemented here.
