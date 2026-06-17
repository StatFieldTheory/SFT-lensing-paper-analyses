# equal_time_limber kappa3 vertex — complete math flow (matter P(k) → ζ)

This document traces the **entire production path** that turns the bare matter
power spectrum P_δ(k, z=0) into the equal-time third-cumulant table ζ_XYZ that
the sft-wick non-local K-vertex consumes. Every computational step is paired
with (a) the equation it implements, (b) the exact `file:line` in code, and
(c) the authoritative paper equation it corresponds to.

Math is written in Unicode/ASCII (chi=χ comoving distance, lambda=λ affine
projection param, ell=ℓ multipole, gamma=γ angular separation, Phi00/Psi0 the
two driving fields, H_0 Hubble, Omega_m, D growth, A Poisson amplitude). Powers
are `^`, subscripts `_`.

Code root: `CAN = /Users/zzhang/projects/canoes/src/canoes`.
Paper root: `PAP = /Users/zzhang/Documents/MyDrafts/STF_lensing/sections`.

---

## 0. Pipeline at a glance

```
                        bare matter  P_δ(k, z=0)            cosmology (Ω_m,h,n_s)
                        PCAMBz0.txt, h-units                 + λ-grid (16 shells)
                              │                                + cosine-triples (1671)
                              ▼
        ┌─────────────────────────────────────────────────────────────┐
        │  build_equal_time_limber_table.py  (conductor; this folder)  │
        └─────────────────────────────────────────────────────────────┘
            │                                          │
            │ ℓ ≤ 60                                    │ 60 < ℓ ≤ 1000
            ▼                                          ▼
   ┌──────────────────────┐                  ┌──────────────────────────┐
   │ LOW  (Stage 2)       │                  │ HIGH (Stage 3)           │
   │ compute_kappa3_      │                  │ compute_kappa3_          │
   │   zeta_table         │                  │   sigma3_high            │
   │ exact Wigner-3j      │                  │ flat-sky Born-Limber     │
   │ FFTlog→I_ℓ→Sachs6×6×6│                  │ k=ℓ/χ, SPT F2, J_S Bessel│
   └──────────┬───────────┘                  └────────────┬─────────────┘
              │ ζ_LOW^XYZ(triple, λ)                        │ ζ_HIGH^XYZ(triple, λ)
              └───────────────┬────────────────────────────┘
                              ▼  Stage 4
                   kappa3_combine_low_high
                   ζ^XYZ = ζ_LOW^XYZ + ζ_HIGH^XYZ   (hard ell_cut split, no taper)
                              ▼  Stage 5
                   Kappa3Output.save_npz  →  equal_time_limber_kappa3_*.npz
                   (zeta_TTT/TTP/TPP/PPP over cosine_triples × λ_shells,
                    physical units = ζ_h · h^6, λ-density measure)
                              ▼  Stage 6
                   equal_time_limber_kappa3_callable.py
                   coupling_fn(n_list, t_list) → (3,3,3) tensor
                   (cKDTree kNN in cos-space + linear-in-λ)
```

The build script itself is only the conductor — it passes arguments and sums
the two branches:
`build_equal_time_limber_table.py:133-151` (LOW call, HIGH call, combine).

---

## 1. Stage 0 — Inputs

### 1.1 Matter power spectrum (the single field-theory primitive)
```
P_δ(k, z=0)   loaded from PCAMBz0.txt   [k in h/Mpc,  P in (Mpc/h)^3]
low-k:   P(k) = P(k_min) (k/k_min)^{n_s}          (primordial tilt, n_s=0.97)
high-k:  raises (no extrapolation)
in-range: log-log linear interpolation
```
- code: `equal_time_limber/_local_cosmo_pk.py:42-53` (`load_pk_delta`, feeds the
  *bare matter* P_δ, not P_Φ); class `_CambTablePk` at `CAN/cosmo/pk.py:85-154`
  (eval `:152-153`, low-k `∝k^{n_s}` `:145-151`, high-k raises `:139-144`).
- paper: P_Φ(k; z=0) is named "the single time-independent primitive"
  (`PAP/appendix.tex:496-500`); the δ↔Φ link is `P_δ(k)=[k²/A(a=1)]² P_Φ(k)`
  in `eq: appendix bphi tree` (`PAP/appendix.tex:608-614`).

### 1.2 Cosmology
```
Ω_m = 0.3160919980475834 = (Ω_c h² + Ω_b h²)/h² = (0.12029+0.02207)/0.6711²
h   = 0.6711,   n_s = 0.97
```
- code: `equal_time_limber/_local_cosmo_pk.py:33-39`; defaults
  `CAN/cosmo/background.py:42-51`.

### 1.3 The radial grid (λ-shells) and the angular grid (cosine-triples)
```
λ-shells:  16 nodes 396.6 … 2326.6 Mpc  (z = 0.1 … 5.7), read from the L2 NPZ
triples:   (cos γ12, cos γ23, cos γ31), 15-node cos grid (cos=k/7), triangle-
           filtered → 1671 realizable triples
```
- code: λ-grid `build_equal_time_limber_table.py:61-69`; triple filter
  `:72-81` (`g12+g23≥g31 & … & g12+g23+g31 ≤ 2π`).

---

## 2. Stage 2 — LOW branch (exact Wigner-3j, ℓ ≤ ell_cut=60)

Entry: `compute_kappa3_zeta_table` `CAN/sachs/kappa3.py:3777`.
The vertex factorizes into a **radial** stage (build the reduced bispectrum
b_{ℓ1ℓ2ℓ3}(λ)) and an **angular** stage (Wigner-3j → ζ).

### 2.1 FFTlog decomposition of P_δ (Nmax = 128 Mellin modes)
```
P_δ(k) ≈ Σ_n c_n k^{ν_n},   ν_n = bias + 2πi n /(Nmax·Δκ),  bias = -1.6,
         Nmax = 128 log-k samples (Hermitian half = 64 modes used)
```
`Nmax` is the number of FFTlog k-samples (NOT a Bessel order, NOT the
bandlimit ell_max).
- code: `CAN/nuell/transfer/fftlog.py:43-131`; kappa3 half-spectrum build
  `CAN/sachs/kappa3.py:458-466`.

### 2.2 Radial double-Bessel transform (the I_ℓ kernel)
```
I_ℓ(ν, t) = 4π ∫_0^∞ dv v^{ν-1} j_ℓ(v) j_ℓ(v t)          (ASZ-2017 Mellin integral)
P_ℓ(r; χ) = Σ_n c_n χ^{-ν_n} I_ℓ(ν_n + shift, r/χ) / (4π)   (power-law "PL" leg)
S_spec_{ℓ}(r; χ)                                            (closure δ² leg, Γ-regularized)
```
- code: I_ℓ def `CAN/nuell/transfer/kernels/il.py:1-6`, entry `il.py:1288`;
  PL leg `CAN/nuell/transfer/_bl_jax.py:679-698`; spec leg `_bl_jax.py:700-717`.
- paper: per-leg windowed transfer Δ_ℓ^X(k;λ) is `eq: appendix delta ell`
  (`PAP/appendix.tex:562-568`); the windowed per-leg kernel M'_ℓ is
  `eq: appendix Mprime` (`PAP/appendix.tex:653-659`).

### 2.3 Tree-level SPT layout + the Sachs 6×6×6 basis contraction
```
b_{ℓ1ℓ2ℓ3}(λ) = (2/π)^3
              · Σ_{k,l,m ∈ basis(6)} (C_a·M)_k (C_b·M)_l (C_c·M)_m
              · Σ_{3 cyclic spec} Σ_{6 F2 monomials} coeff
              · ∫_0^∞ dr r^2  S_spec_{ℓc}(r;χ) P_{ℓa}(r;χ) P_{ℓb}(r;χ)
              · Δ_tri(ℓ1,ℓ2,ℓ3)                              (triangle mask)
```
The 18 = (3 cyclic spectator) × (6 F2 monomials) terms encode the SPT
tree matter bispectrum `B_δ = 2 F2 P P + cyc`. `(C_a, C_b, C_c)` are the per-leg
Sachs CoefVecs (Stage 2.6); `M` is the lightcone derivative basis (Stage 2.5).
- code: r-integral leaf `_bl_jax.py:914-917`; 6×6×6 contraction (fused multi
  channel) `_bl_jax.py:1081-1087`; SPT monomials `CAN/nuell/transfer/
  _phi_accumulator.py:49-69`; `(2/π)^3` `_bl_jax.py:966-971`.
- paper: **this is exactly `eq: appendix b prime`** (`PAP/appendix.tex:667-675`):
  `B'^{XYZ}_{LLℓ3} = (2/π)^3 ∫dr r^2 Σ_{t=0}^{17} c_t M'_spec M'_PL M'_PL`.

### 2.4 The cosmology prefactor stack (applied per λ-shell)
The dimensionless [angular × SPT-shape] integral is multiplied by, in order:

| # | factor | equation | code |
|---|--------|----------|------|
| i | Poisson³ | `A(a)^3 = [−(3/2) H_0² Ω_m (1+z)]^3`, H_0 = 100/c [h/Mpc] (h-independent), units (h/Mpc)^6 | `kappa3.py:85, 4222-4224`, cubed `:4268-4269` |
| ii | Sachs geometric | `prefactor = ((1+z)^P)^3 = (1+z)^12`, `P = SACHS_GEOMETRIC_POWER = 4` | `kappa3.py:499`; `CAN/nuell/transfer/sachs_driver.py:97` |
| iii | per-leg M-factors | scale vectors `sX = C_X · M`, `M[0]=D(1+z)=D/a`, `M[1]=−d_χ[D(1+z)]`, `M[2]=d_χ²[D(1+z)]` (each linear in D) | `kappa3.py:482-484, 512`; `CAN/cosmo/lightcone.py:342, 346-350` |
| iv | spectator growth | `× spec_growth = M[0]/(1+z) = D(z)` (promotes D³ → D⁴) | `kappa3.py:4495-4506` |
| v | radial measure | `radial_measure="lambda" → × 1` (native λ-density, **no (dχ/dλ) Jacobian**) | `kappa3.py:1761-1805, 4508-4520` |
| vi | units | `units="physical" → ζ = h^6 · ζ_h`, χ→χ/h, λ→λ/h | `kappa3.py:4522-4530` |

- paper: per-leg Poisson `A(a)=−(3/2)H_0² Ω_m/a` and per-leg `(1+z)^4` Sachs
  weight in `PAP/appendix.tex:603-604, 536-543`; λ-density measure
  `dλ = a²(χ) dχ` in `eq: cumulant 2 measure change`
  (`PAP/cosmology.tex:979-988`); `D̄=a·χ` in `eq: D equals a chi`
  (`PAP/appendix.tex:286-289`).

> ⚠ Redshift scaling is NOT a clean power. See §7. The growth is clean
> (**D⁴**), but the (1+z) content is heterogeneous because COEFS_PHI00 (2.6)
> contracts M[0], M[1], M[2] together with the ℋ²/ℋ′/L²-χ⁻² operator terms.

### 2.5 The lightcone derivative basis M[0..2]
```
M[0] = D(z)·(1+z) = D/a            (one growth amplitude + (1+z))
M[1] = −d/dχ [ D(1+z) ]            (linear in D; χ-derivative ⇒ ℋ-content)
M[2] = +d²/dχ² [ D(1+z) ]
ℋ ≡ calH = E(z)/((1+z)·(c/100)),   E(z)=sqrt(Ω_m(1+z)^3 + 1−Ω_m)
```
- code: `CAN/cosmo/lightcone.py:328-329, 342-350`.

### 2.6 The per-leg Sachs CoefVecs (what makes T vs P)
```
COEFS_PHI00 = ( L²/χ² + 2ℋ² − 2ℋ′,  −2/χ + 2ℋ,  −1,  0,  −1,  2 )   # T = Φ00, spin 0
COEFS_PSI0  = ( −sqrt(L²(L²−2))/χ²,   0, 0, 0, 0, 0 )                 # P = Ψ0, spin 2
   with  L² = ℓ(ℓ+1)
```
A leg's operator value is `g_leg(ℓ,χ) = Σ_{k=0}^{5} COEFS[k] · M[a_k]`
(the basis index → τ-derivative order a_k = [0,0,0,1,2,1] gathers M).
- code: `CAN/nuell/transfer/sachs_driver.py:74-90`; gather `kappa3.py:482-484`;
  `BASIS_TO_AB`, `L²=ℓ(ℓ+1)` `CAN/nuell/transfer/basis.py:45-52, 111-120`.

### 2.7 Angular stage — exact Wigner-3j → ζ
For each cosine-triple (γ12, γ23, γ31), with channel spins (s1,s2,s3):
```
ζ_LOW^XYZ(γ12,γ23,γ31; λ) = D(z) · Σ_{ℓ1ℓ2ℓ3 ≤ 60}  b_{ℓ1ℓ2ℓ3}(λ)
    · h · w3j(ℓ1,ℓ2,ℓ3; 0,0,0) · sqrt[(2ℓ1+1)/4π] · (−1)^{s1}
    · Σ_{m2}  w3j(ℓ1,ℓ2,ℓ3; −s1, m2, s1−m2)
            · {}_{s2}Y_{ℓ2 m2}(γ12, 0) · {}_{s3}Y_{ℓ3, s1−m2}(γ31, φ3)

h = sqrt[(2ℓ1+1)(2ℓ2+1)(2ℓ3+1)/4π]
```
Two Wigner-3j symbols appear: the parity/triangle anchor `(ℓ1ℓ2ℓ3;0,0,0)`
(kept for ALL four channels) and the m-resolved `(ℓ1ℓ2ℓ3;−s1,m2,s1−m2)`.
The angular separations enter through the spin-weighted spherical harmonics
with θ = arccos(cos γ). The leading `D(z)` is the spectator growth (2.4-iv).
- code: scalar (TTT) `CAN/nuell/transfer/three_pt.py:516-645`; spin-weighted
  `CAN/nuell/transfer/_spin_aware_three_pt.py:15-23, 200-289`; GEMV apply
  `CAN/nuell/transfer/_jax_kernels.py:1012-1032`.
- paper: **this is `eq: appendix zeta squeezed`** (`PAP/appendix.tex:684-704`),
  `ζ_XYZ = Σ_{L,ℓ3} G(L,L,ℓ3) P_{ℓ3}(cos γ) B'^{XYZ}_{LLℓ3}`, with the Gaunt
  weight `G = (2L+1)²(2ℓ3+1)/(4π)² [w3j(L,L,ℓ3;0,0,0)]²`. **The paper writes the
  squeezed (L,L,ℓ3) specialization** (Legendre P_{ℓ3}(cos γ), two rays
  coincident); the code's general m-sum reduces to it for that config.

---

## 3. Stage 3 — HIGH branch (flat-sky Born-Limber, 60 < ℓ ≤ 1000)

Entry: `compute_kappa3_sigma3_high` `CAN/sachs/kappa3.py:3406-3684`.

### 3.1 Flat-sky Limber collapse (k = ℓ/χ, single shell)
```
two flat-sky vectors: u=|ℓ2|, v=|ℓ3|, φ = angle(ℓ2,ℓ3);  closure |ℓ1| = w
w = sqrt(u² + v² + 2uv cos φ)
k_i = ℓ_i / χ(λ)                     (equal-χ Limber, NO factor of a)
radial Limber kernel = 1/χ^4         (three 1/χ² collapsed by two equal-χ δ's)
```
- code: quadrature `kappa3.py:1608-1649` (Gauss-Legendre n_ell=96 in ℓ,
  n_phi=64 in φ); k=ℓ/χ `:3570-3572`; 1/χ⁴ `:3585`; closure w `:3513`.

### 3.2 Tree-level SPT matter bispectrum (F2 kernel)
```
B_δ(k1,k2,k3) = 2 [ F2(k1,k2)P1P2 + F2(k1,k3)P1P3 + F2(k2,k3)P2P3 ]
F2(k_a,k_b) = 5/7 + (1/2) μ (k_a/k_b + k_b/k_a) + (2/7) μ²,   μ = cos(k_a,k_b)
P_i = P_δ(k_i, z=0)                  (bare matter, h-units)
```
- code: F2 `kappa3.py:1690-1702`; B assembly `:3573-3579`; internal cosines
  from the closed ℓ-triangle `:3530-3540`.
- paper: same SPT tree `B^Φ_tree` (with F2) in `eq: appendix bphi tree`
  (`PAP/appendix.tex:608-614`).

### 3.3 Cosmology factors (HIGH mirror of 2.4)
```
growth:   × D(z)^4         (B_δ ∝ P·P ∝ D^4 at tree level; growth=M[0]/(1+z)=D)
Poisson + Sachs per leg:  R_channel = Π_legs (±) A(a)(1+z)^4,
          A(a) = −(3/2) Ω_m H_0² (1+z);  sign +: Φ00 (T), −: Ψ0 (P)
prefactor: 1/(2π)^4        (the d²ℓ_a d²ℓ_b normalization)
```
- code: growth⁴ `kappa3.py:3580-3583`; A(a) `:3498-3500`; response product
  `_kappa3_limber_response_product_h` `:1669-1687`, applied `:3589-3594`;
  1/(2π)⁴ `:3541`.

### 3.4 Spin via the analytic orientation (Bessel) integral
```
position triangle: θ_ij = arccos(cos γ_ij);  r2 = θ12, r3 = x3 + i y3 (law of cos)
Q = u·r2 + v·e^{−iφ}·r3,   S = s1+s2+s3
α_channel = Re[ 2π i^S e^{iS argQ} J_S(|Q|) · e^{i(s1 arg ℓ1 + s3 φ)} ]
  (S=0 ⇒ α = 2π J_0(|Q|), real)
```
spins: TTT (0,0,0)→S=0; TTP (0,0,2)→S=2; TPP (0,2,2)→S=4; PPP (2,2,2)→S=6.
All four channels are non-zero in HIGH.
- code: triangle vertices `kappa3.py:1705-1717`; α-kernel
  `_kappa3_limber_alpha_kernel` `:1720-1758`; spins
  `_kappa3_limber_channel_spins` `:1664-1666`.

### 3.5 HIGH assembled
```
ζ_HIGH^XYZ(triple,λ) = (1/(2π)^4)(D^4/χ^4) R_channel
   · ∫_{ell_cut < max(ℓ) ≤ ell_high_max} du dv dφ  (u v · w_ℓ w_ℓ w_φ)  B_δ  α_channel
```
- code: radial_base `:3585`, weighted×α sum `:3595-3596`; HIGH window
  `(max_leg>ell_cut)&(max_leg≤ell_high_max)` `:3515-3518`.

---

## 4. Stage 4 — Combine LOW + HIGH

```
ζ^XYZ(triple, λ) = ζ_LOW^XYZ(triple, λ) + ζ_HIGH^XYZ(triple, λ)
```
A plain elementwise array sum on the IDENTICAL (cosine-triple, λ) grid.
**Hard split at ell_cut, no tapering / no overlap blend**: LOW does exact 3j for
ℓ ≤ ell_cut, HIGH does flat-sky for max(ℓ) ∈ (ell_cut, ell_high_max]; the
boundary max(ℓ)=ell_cut is in LOW only. Grid/units/measure compatibility is
asserted before the add; cosmo_meta is merged.
- code: `kappa3_combine_low_high` `kappa3.py:3710-3774` (sum `:3766-3771`,
  grid checks `:3741-3742`, meta merge `:3750-3764`).

---

## 5. Stage 5 — Output table

`Kappa3Output` holds `cosine_triples (M,3)`, `lambda_shells_Mpc`,
`chi_shells_Mpc`, `z_shells`, `zeta_TTT/TTP/TPP/PPP (M, N_λ)`, `ell_max`,
`cosmo_meta`. `save_npz` writes them in:
```
units  = physical:  λ,χ in Mpc; ζ = h^6 · ζ_h
measure = lambda :  ζ is a per-leg λ-density three-leg cumulant (no Jacobian)
```
- code: dataclass `kappa3.py:164-237`; `save_npz` `:187-206`.
- This is the NPZ the callable loads. The build then overwrites the metadata
  block (`build_equal_time_limber_table.py:152-170`) recording
  `already_R_contracted=False`, `equal_time=True`, the FK baseline string, etc.

---

## 6. Stage 6 — Callable consumption (what sft-wick sees)

`coupling_fn(n_list, t_list)` turns three sky directions + three λ positions
into the (3,3,3) driving-field-cumulant tensor:
```
cos γ_ij = n_i · n_j  →  sort descending  →  cKDTree query (k=4 inverse-distance)
λ_mean = mean(t_list) →  linear interp between the two bracketing λ-shells
ζ_channel = Σ_k w_k ζ_channel[node_k, λ] / Σ_k w_k        (exact node if dist<1e-10)
tensor:  out[0,0,0]=ζ_TTT; out[{0,0,1}perm]=ζ_TTP; out[{0,1,1}perm]=ζ_TPP;
         out[1,1,1]=ζ_PPP; all index-2 (B-mode) entries = 0
```
- code: `equal_time_limber_kappa3_callable.py:129-170` (scalar),
  `:173-228` (batch). `already_R_contracted=False` ⇒ sft-wick applies the
  response R itself (this table is the bare vertex K, not R-folded).
- See `slices/` for how ζ behaves over (γ, z): a contact spike at the collapsed
  apex cos=(1,1,1) sitting ~100× above a flat wide-separation plateau, growing
  toward high z.

---

## 7. Redshift scaling of ζ (the honest statement)

ζ_TTT(apex, z) does **NOT** factor as a clean `(1+z)^P · D(z)^Q × z-independent
shape`. Breakdown:

- **Growth is clean: Q = 4.** Three per-leg M-factors (each linear in D) ×
  one spectator D = D⁴. No extra D hides elsewhere (the χ-derivatives in
  M[1], M[2] differentiate the *shape* D(χ)(1+z), not the amplitude count).
  Metadata `growth_scaling="D^4"` is correct on the D-count.
- **(1+z) is NOT clean.** Explicit factors give a clean `(1+z)^15` =
  `(1+z)^12` (Sachs prefactor 2.4-ii) × `(1+z)^3` (Poisson cube 2.4-i). But the
  per-leg operator `Σ COEFS_PHI00[k]·M[a_k]` mixes M[0]=D(1+z) with M[1], M[2]
  (χ-derivatives, carrying ℋ-content) and adds the `2ℋ²−2ℋ′` and `L²/χ²` terms
  inside one slot. ℋ = E(z)/((1+z)·c/100) is not a power of (1+z). So the exact
  (1+z) exponent is not well defined.
- **Leading-term scaling** (keep only the ID·M[0] `L²/χ²` piece, dominant at
  ℓ≳ a few):
  ```
  ζ_TTT(apex,z) ≈ [−(3/2)H_0² Ω_m]^3 · (1+z)^15 · D(z)^4 · χ^{−6}
                  · [Π_i ℓ_i(ℓ_i+1)] · [Wigner-3j] · [SPT ∫r² dr shape]
  ```
  Use this as the *leading* scaling, not an exact factorization.

This explains the §6 slices: ζ grows steeply toward high z (the (1+z)^15·D⁴
front factor), and is sharply peaked at the collapsed apex (the Wigner-3j /
SPT-shape angular structure).

---

## 8. REVIEW — production step ↔ equation completeness

Every step in the production path, with its equation and a verdict on whether
an equation is attached.

| step | production action | equation | code | paper | covered? |
|------|-------------------|----------|------|-------|----------|
| S0.1 | load matter P_δ | §1.1 (P_δ, n_s tilt) | `_local_cosmo_pk.py:42-53`; `pk.py:85-154` | `appendix.tex:496-500, 608-614` | ✅ |
| S0.2 | cosmology / λ-grid / triples | §1.2-1.3 | `build…:61-81` | — | ✅ |
| S2.1 | FFTlog of P_δ | §2.1 `P_δ≈Σc_n k^{ν_n}` | `fftlog.py:43-131` | (numeric, no paper eq) | ✅ |
| S2.2 | radial double-Bessel I_ℓ | §2.2 `I_ℓ`, `P_ℓ(r;χ)` | `il.py:1-6`; `_bl_jax.py:679-717` | `eq: appendix delta ell`, `Mprime` | ✅ |
| S2.3 | SPT 18-term + 6×6×6 contraction | §2.3 `b_{ℓ1ℓ2ℓ3}` | `_bl_jax.py:914-1087`; `_phi_accumulator.py:49-69` | **`eq: appendix b prime`** | ✅ |
| S2.4i | Poisson A(a)³ | §2.4-i | `kappa3.py:4222-4224, 4269` | `cosmology.tex:405-407`; `appendix.tex:603-604` | ✅ (sign caveat R1) |
| S2.4ii | (1+z)^12 Sachs prefactor | §2.4-ii | `kappa3.py:499`; `sachs_driver.py:97` | `appendix.tex:536-543` | ✅ |
| S2.4iii | per-leg M-factors (D/a etc.) | §2.5, 2.6 | `kappa3.py:482-512`; `lightcone.py:342-350` | `eq: appendix Mprime` | ✅ |
| S2.4iv | spectator D (→D⁴) | §2.4-iv | `kappa3.py:4495-4506` | — (code convention, R3) | ✅ (caveat R3) |
| S2.4v | λ-density measure (×1) | §2.4-v | `kappa3.py:1761-1805` | `eq: cumulant 2 measure change`; `D=aχ` | ✅ |
| S2.4vi | h⁶ unit conversion | §2.4-vi | `kappa3.py:4522-4530` | (h-units bookkeeping) | ✅ |
| S2.7 | exact Wigner-3j angular | §2.7 | `three_pt.py:516-645`; `_spin_aware_three_pt.py:15-289` | **`eq: appendix zeta squeezed`** + `G` | ✅ (spin caveat R2) |
| S3.1 | flat-sky Limber k=ℓ/χ, 1/χ⁴ | §3.1 | `kappa3.py:3570-3585, 1608-1649` | (Limber, standard) | ✅ |
| S3.2 | SPT F2 matter bispectrum | §3.2 | `kappa3.py:1690-1702, 3573-3579` | `eq: appendix bphi tree` | ✅ |
| S3.3 | HIGH cosmology factors | §3.3 | `kappa3.py:3498-3500, 3580-3594, 1669-1687` | `appendix.tex:603-604` | ✅ |
| S3.4 | spin Bessel orientation J_S | §3.4 | `kappa3.py:1705-1758` | (flat-sky 3PCF, no closed paper eq) | ✅ |
| S4 | LOW+HIGH sum | §4 | `kappa3.py:3710-3774` | (numerical ℓ-split) | ✅ |
| S5 | NPZ save (units, measure) | §5 | `kappa3.py:164-206` | — | ✅ |
| S6 | callable → (3,3,3) | §6 | `…callable.py:129-228` | sft-wick contract | ✅ |

**Verdict:** every production step has a corresponding equation. Three points
are code↔paper *conventions to keep on the radar* (not gaps in the flow):

- **R1 — Poisson sign.** `cosmology.tex:405-407` writes Φ00=+Aδ with A>0;
  `appendix.tex:603-604` and the code use A(a)=−(3/2)H_0²Ω_m/a (A<0). Because
  A enters cubed (odd power), the sign propagates to ζ. The code is internally
  consistent with the appendix; the main-text caption sign differs. Worth a
  one-line reconciliation in the paper, not a code change.
- **R2 — spin-2 Wigner-3j.** Both code and paper keep the *scalar* parity
  anchor `(ℓ1ℓ2ℓ3;0,0,0)` for ALL four channels including the spin-2 (Ψ0=P)
  legs (only the m-sum 3j carries the −s1/s1−m2 shifts). This is the documented
  "pragmatic W3j_0 anchor" (project memory `project_canoes_patches_done`).
  Code and paper AGREE; flagged because both are the pragmatic convention, not
  the full (±2,…) triple-spin anchor.
- **R3 — explicit D⁴ growth.** The code multiplies an explicit linear growth
  D(z)⁴ (3 per-leg + spectator). The authoritative `\input`-ed paper sections
  do NOT write a `D⁴` law; they encode time-dependence via the per-leg
  A(a)(1+z)^4 windows on the z=0 primitive P_Φ. Physically D⁴ is the tree-level
  bispectrum scaling (B_δ ∝ P·P ∝ D⁴), so the two are consistent, but the
  explicit D⁴ is a code-side bookkeeping choice. (An explicit "D⁴/per-leg
  time-factor" statement exists only in the untracked `cosmology2.tex:966`.)

---

## 9. End-to-end composite (LOW, one channel)

```
P_δ(k,z=0)
  → FFTlog(Nmax=128):  P_δ(k) ≈ Σ_n c_n k^{ν_n}
  → P_ℓ(r;χ) = Σ_n c_n χ^{−ν_n} I_ℓ(ν_n+shift, r/χ)/(4π),  I_ℓ(ν,t)=4π∫dv v^{ν−1}j_ℓ(v)j_ℓ(vt)
  → b_{ℓ1ℓ2ℓ3}(λ) = (2/π)^3 · A(a)^3 · (1+z)^12
                    · Σ_{basis} (C_a·M)(C_b·M)(C_c·M)
                    · Σ_{3 cyc spec}Σ_{6 F2 mono} coeff ∫dr r² S_spec_{ℓc} P_{ℓa} P_{ℓb}
  → ζ_LOW^XYZ(γ12,γ23,γ31;λ) = D(z) · Σ_{ℓ1ℓ2ℓ3≤60} b_{ℓ1ℓ2ℓ3}
       · h·w3j(…;0,0,0)·sqrt[(2ℓ1+1)/4π]·(−1)^{s1}
       · Σ_{m2} w3j(…;−s1,m2,s1−m2)·{}_{s2}Y_{ℓ2 m2}(γ12)·{}_{s3}Y_{ℓ3,s1−m2}(γ31)
  → (× 1 for λ-measure;  × h^6 for physical units)

ζ_HIGH^XYZ = (1/(2π)^4)(D^4/χ^4) R_channel ∫_{60<ℓ≤1000} du dv dφ (uv·w_ℓ w_ℓ w_φ) B_δ α_channel,
       B_δ = 2[F2 P P + cyc],  α_channel = Re[2π i^S e^{iS argQ} J_S(|Q|) e^{i(s1 argℓ1 + s3 φ)}]

ζ^XYZ(triple, λ) = ζ_LOW^XYZ + ζ_HIGH^XYZ        (hard split at ℓ=60)
```

Channels: T=Φ00 (spin 0, COEFS_PHI00); P=Ψ0 (spin 2, COEFS_PSI0).
TTT/TTP/TPP/PPP = (PHI00,PHI00,PHI00)/(PHI00,PHI00,PSI0)/(PHI00,PSI0,PSI0)/(PSI0,PSI0,PSI0),
spins (0,0,0)/(0,0,2)/(0,2,2)/(2,2,2).

---

## Provenance
- Code traced at `CAN = /Users/zzhang/projects/canoes/src/canoes` (commit state
  of 2026-06-02).
- Authoritative paper equations: `eq: appendix b prime`, `eq: appendix zeta
  squeezed` (`PAP/appendix.tex:667-704`); Poisson `PAP/cosmology.tex:405-407,
  561-570`; Sachs (1+z)^4 `PAP/appendix.tex:536-543`; measure
  `PAP/cosmology.tex:979-988`; `D=aχ` `PAP/appendix.tex:286-289`.
- Reproduce the table: `python build_equal_time_limber_table.py`.
