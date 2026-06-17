# Shear-3PCF natural-component map: driving-field cumulants -> Gamma^0..3

Derivation: `scripts/derive_shear3pcf_natural_components.wl` (plain symbolic
spin-weighted algebra; runs clean, all geometric/limit checks below).
Consumer: `stage1_ours_v2.py`. Comparison: `stage1_compare_v2.py` ->
`outputs/stage1_compare_v2.npz`.

This supersedes the v1 build (`stage1_ours_v1_archived_2026-06-09.py`), whose
projection phase used fastnc's `x2cent` (the x-projection -> centroid bridge),
which is NOT the correct bridge from canoes' canonical great-circle frame.

---

## 1. The physical chain (sources pinned)

1. **Shear <- spin-2 driving field.** `cosmology.tex` (`eq: gamma from Dsachs23`):
   `gamma = gamma_+ + i gamma_x = -int dl (Dsachs2 + i Dsachs3) = -int dl Psi0`,
   where `Psi0 = Psi_+ + i Psi_x` is the spin-(+2) screen-Hessian driving field
   (`appendix.tex` `eq: appendix screen hessian`, multiplier
   `sqrt(L^2(L^2-2))/chi^2`). The spin-2 eth^2 multiplier is baked into the
   canoes modulus channels; not re-applied.

2. **Natural components (Schneider-Lombardi; Sugiyama+2024 = fastnc).**
   With `g|_a = gamma e^{-2 i a}` the shear projected onto reference angle `a`:
   - `Gamma^0 = <g|_a1  g|_a2  g|_a3 >`  (all unconjugated)
   - `Gamma^1 = <g*|_a1 g|_a2  g|_a3 >`  (leg 1 conjugated)
   - `Gamma^2 = <g|_a1  g*|_a2 g|_a3 >`  (leg 2 conjugated)
   - `Gamma^3 = <g|_a1  g|_a2  g*|_a3>`  (leg 3 conjugated)
   `a_j` = CENTROID projection angle at vertex j.

3. **Helicity bookkeeping confirmed against fastnc.** fastnc's Bessel-order map
   (`fastnc.py:365`) at multipole M=0 gives `(m,n)` per mu =
   `{(-3,-3),(-1,-1),(1,-3),(-3,1)}`, i.e. total leg-spin `{-6,-2,-2,-2}`:
   **Gamma^0 = all-same-helicity (total spin 6), Gamma^{1,2,3} = one helicity
   flipped (total spin 2).** This fixes the channel assignment:
   - `Gamma^0 <- zeta_PPP` (`<P1 P2 P3>`, all +2; gamma^4-suppressed, `d^L_{2,-2}`),
   - `Gamma^{1,2,3} <- zeta_D` (`<...>` with leg mu conjugated; modulus, `d^L_{2,2}~O(1)`).
   This is the SAME channel logic the appendix uses at 2-point: `xi_+ <- zeta_B`
   (modulus, finite), `xi_- <-` un-conjugated (suppressed) [appendix 606-622].

4. **canoes canonical great-circle frame.** `_spin_aware_three_pt.py` reduces
   `zeta^{s1s2s3}` with `n1` at the pole, `n2` in the xz-plane, `n3` at azimuth
   `phi3`, spin reference = great-circle tangent, `{}_sY_{lm} = sqrt[(2l+1)/4pi]
   d^l_{m,-s}(theta) e^{i m phi}`. canoes returns the REAL cumulant (parity kills
   the imaginary part). The great-circle projection is BAKED into the harmonic
   sum -- it is not an external phase.

---

## 2. The derived map (centroid projection, SAS isoceles family)

SAS triple `(cos t3, cos g, cos g)`: legs 1,2 = base endpoints (separated by
`t3 = 2 g sin(phi/2)`), leg 3 = apex.

**Centroid reference angles** `a_j` (vertex -> centroid):
```
a1 = -arccot(3 tan(phi/2))
a2 = -pi + arccot(3 tan(phi/2))
a3 = pi/2
```
**Canonical great-circle reference angles** `b_j` (geodesic to leg 1):
```
b1 = 0
b2 = pi
b3 = pi - arctan(cot(phi/2))
```
**Per-leg great-circle -> centroid spin-2 rotation** `(a_j - b_j)` (pure phi,
g-independent -- verified `d/dg = 0`):
```
a1 - b1 = -arccot(3 tan(phi/2))
a2 - b2 = -2pi + arccot(3 tan(phi/2))
a3 - b3 = -pi/2 + arctan(cot(phi/2))
```
**The map** (`s_j^mu = +1` unconjugated, `-1` conjugated leg):
```
Gamma^mu = Z_channel(mu) * exp(-2 i sum_j s_j^mu (a_j - b_j))
Gamma^0 = zeta_PPP * exp(-2i[ (a1-b1)+(a2-b2)+(a3-b3) ])   phase = e^{+i phi}
Gamma^1 = zeta_D1  * exp(-2i[-(a1-b1)+(a2-b2)+(a3-b3) ])
Gamma^2 = zeta_D2  * exp(-2i[ (a1-b1)-(a2-b2)+(a3-b3) ])
Gamma^3 = zeta_D3  * exp(-2i[ (a1-b1)+(a2-b2)-(a3-b3) ])   phase = e^{-i phi}
```
with the overall radial fold `Z = -int_0^{l_s} dl K^3 zeta`,
`K = a^2 chi (chi_s-chi)/chi_s`, `Dbar = a chi`, `z_s = 5`, `radial_measure=lambda`.

---

## 3. Limit checks (from the .wl)

- **(i) Magnitudes are phase-invariant.** `|exp(-2i sum ...)| = 1` exactly, so
  `|Gamma^mu| = |Z_channel(mu)|`. The projection CANNOT redistribute magnitude
  between components -- a frame rotation is ruled out as the fix (lead's
  guardrail). PASS.
- **(ii) Small-gamma channel structure.** `a_j - b_j` are pure functions of phi
  (scale-free; `d/dg = 0`). So at small gamma the magnitude structure is carried
  entirely by the cumulants: `zeta_PPP ~ gamma^4` (`d^L_{2,-2}`) suppresses
  `|Gamma^0|`, modulus `zeta_D ~ O(1)` (`d^L_{2,2}`) keeps `|Gamma^{1,2,3}|`
  finite. Matches appendix 615-619. PASS.
- **(iii) 2-point xi_+/xi_- analogue.** One-conjugation (xi_+) picks the finite
  modulus `zeta_B`; no-conjugation (xi_-) picks the gamma^4-suppressed
  un-conjugated channel. Reproduces appendix 606-622 ("lifts xi_+ from zero to a
  comparable level"). PASS.
- **(ii') Isoceles symmetry.** Apex-conjugated (leg 3) phase is `e^{-i phi}`
  (real combination on the mirror), consistent with fastnc's `Gamma`-real /
  conjugate-pair structure for `t1=t2`. The phase-only mirror identity has a
  residual `(e^{-2i phi}-1)` term: see the CAVEAT below.

---

## 4. CAVEATS (reported honestly, no fitting)

**VERDICT (stage1_compare_v2.npz): NOT CLOSED.** PRIMARY point (gamma=10',
phi=60deg): ours_v2 Gamma^0=-1.72e-8-2.98e-8j, Gamma^1=-3.54e-7+6.12e-7j,
Gamma^2=-4.28e-7, Gamma^3=+2.15e-7-3.72e-7j vs fnc Gamma^0=+6.81e-8,
Gamma^1=+1.38e-7, Gamma^2=Gamma^3*=+1.24e-7. |ratio| = {0.50, 5.14, 3.45, 3.46};
frame-invariant ratio 3.99. Grid: frame-invariant ratio 3.06-17.36, median 8.04,
scatter std/mean 0.43; per-component |ratio| medians 4.7-6.2 with range 0.01-28x.
The residual is STRUCTURED (gamma- and phi-dependent), NOT a constant the
projection or a normalization could absorb.

1. **Magnitude does NOT close.** This is the genuine, dominant finding. The
   per-component `|Gamma^mu|` (phase-invariant, set by the cumulants) disagree
   with fastnc tree-level by structured factors (frame-invariant ratio ~3-17x,
   median ~8x; per-component 0.1-20x, growing with gamma). See
   `outputs/stage1_compare_v2.npz`. A pure projection phase cannot fix this.
   **Diagnosis:** fastnc builds ALL FOUR natural components from the SAME scalar
   convergence bispectrum `b_L` (`bispectrum.py:746` `kappa_bispectrum_multipole`)
   via the spin-2-2-2 mode-coupling `G_LM = 4pi int dx P_L(cos x) cos[2bbar+Mx]`
   (`coupling.py:250` `MCF222LegendreFourier`) -- a fixed-frame Fourier-space
   spin-2 projection. Our driving-field cumulants are the great-circle-frame
   `<Psi0 Psi0 Psi0>`-type objects with the `sqrt(L^2(L^2-2))/chi^2` per-leg
   multiplier and the exact-Wigner-3j (ell<=60) + Limber (ell<=1000) split. These
   are genuinely different projections of the same `B_kappa`; they agree in
   ORDER (the spin-2 channels are the right family -- a naive scalar `zeta_TTT`
   fold is ~40x larger than fastnc's frame-invariant, vs ~3-17x for the spin-2
   channels) but not in the configuration-dependent details. The gamma-growing
   residual points to the discrete-Wigner-3j (low-ell) vs flat-sky-FFTLog
   (continuous, high-ell) branch difference dominating at larger separation.

   **Decisive cross-check (rules out a constant-factor fix).** The ratio
   `fnc_frame_invariant / |fold(zeta_TTT)|` (scalar convergence 3PCF) over the
   grid is NOT constant: range 0.02-25.6, median 0.48, scatter std/mean = 3.35.
   So fastnc's shear 3PCF is neither (a) a constant spin-2 factor times our
   scalar convergence 3PCF, nor (b) a constant times our great-circle spin-2
   channels (frame-invariant ratio 3-17x, also non-constant). Both candidate
   "single normalization" fixes are EXCLUDED by the data. The shear natural
   components are a genuinely different, configuration-dependent, fixed-frame
   Fourier-space spin-2 projection (fastnc's 2DFFTLog `MCF222` mode coupling),
   not reducible to our real-space great-circle cumulants by any scalar factor
   or frame phase.

2. **Phase has a residual great-circle-frame-convention ambiguity.** The derived
   `Gamma^0` phase is `e^{+i phi}`, but fastnc's `Gamma^0` is ~real on the
   isoceles family. The `b_j` (canonical great-circle reference) carries a
   per-leg azimuthal-origin convention (vertex 1 sits at the pole; its spin
   reference is the canonical azimuth-0 axis) that is not fixed by the flat-sky
   geometry alone. Because the magnitude is the dominant, phase-independent
   residual, this phase offset is sub-leading to the closure question; it is
   flagged rather than tuned (tuning a phase to match fastnc is forbidden).

3. **Not a reimplementation of fastnc.** Building fastnc's full 2DFFTLog spin-2
   mode-coupling on our side would be an independent reimplementation of
   Sugiyama+2024 and is out of scope. The transparent route here -- the
   driving-field spin-2 channels with the derived projection -- captures the
   correct family and helicity structure but leaves the structured magnitude
   residual documented above.

---

## 5. Source references

- `sections/appendix.tex` 528-622 (`eq: appendix screen hessian`,
  `eq: appendix zeta squeezed`, `eq: appendix modulus reconstruction`).
- `sections/cosmology.tex` 940-974 (`eq: gamma from Dsachs23`, xi_pm defs).
- fastnc: `fastnc/fastnc.py` 350-409 (GammaM, mu->(m,n) map), 700-784
  (x2ortho/x2cent projection phases), `fastnc/coupling.py` 250-268
  (MCF222LegendreFourier `G_LM`), `fastnc/bispectrum.py` 697-782
  (decompose, kappa_bispectrum_multipole -- the SCALAR `b_L` source).
- canoes: `src/canoes/nuell/correlation/_spin_aware_three_pt.py` 14-130
  (canonical frame + spin-weighted SH), `three_pt.py:163`
  (`_canonical_frame_angles`), `src/canoes/sachs/kappa3.py` 3737-4006
  (`compute_kappa3_mod_zeta_equal_time`, `_spin_weights_dmod`).
- Project memory: `project_kappa3_spin2_helicity_bug_2026-06-04.md`
  (the 2-point helicity fix this generalizes).
