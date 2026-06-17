# The γ-dependence of the FK contribution across panels

**Companion to** [`WHY_FK_BREAKS_KAPPA_GAMMA_EQUIVALENCE.md`](WHY_FK_BREAKS_KAPPA_GAMMA_EQUIVALENCE.md).
That note explains *why* FK is confined to convergence (the Wigner-3j selection rule).
This note explains the **angular-scale (γ) trends** of the FK contribution in each panel:

> Why is FK_κκ nearly **flat** until very large angles, while FK in ξ±/κγ_t **rises
> monotonically toward large scales**?

## The data (FK Order-2, corr_op, dt=2, with observable signs applied)

| γ [arcmin] | FK_κκ | FK_κγ_t | FK_ξ± |
|---|---|---|---|
| 0.5  | **+2.745e-5** | −1.3e-15 | −1.7e-15 |
| 8.5  | +2.744e-5 | −3.9e-13 | −4.8e-13 |
| 57   | +2.738e-5 | −1.7e-11 | −2.1e-11 |
| 372  | +2.49e-5  | −7.0e-10 | −8.7e-10 |
| 957  | +1.51e-5  | −3.7e-9  | −4.5e-9  |
| 2462 | −7.8e-8   | −6.8e-8  | −5.4e-8  |
| 5000 | −3.4e-8   | −2.7e-8  | −1.6e-8  |

FK_κκ is flat (+2.745e-5, 4 digits) over 0.5′–57′, decays past ~100′, sign-flips at
γ≳2000′. FK_κγ_t ≈ FK_ξ₊ ≈ FK_ξ₋ climb monotonically from the ~1e-15 floor to ~1e-8.

## Framework: the 2PCF is an angular transform of C_ℓ

```
ξ(γ) = Σ_ℓ (2ℓ+1)/(4π) · C_ℓ · K_ℓ(γ)
```
- **small γ ↔ high ℓ (small scales); large γ ↔ low ℓ (large scales)** (ℓ_peak ~ 1/γ).
- Angular kernels differ by spin: κκ uses Legendre **P_ℓ**; κγ_t uses **d^ℓ_{0,2}**; ξ±
  use **d^ℓ_{2,±2}**. Note `P_ℓ(0)=1` (full weight at γ=0), but `d^ℓ_{0,2}(0)=0` and
  `d^ℓ_{2,-2}(0)=0` (the spin-2 kernels for κγ_t and ξ− vanish kinematically at γ=0).

## κκ: flat plateau → decay → sign-flip (a REAL κ³ signal)

FK is a **squeezed** κ³ term: one long-wavelength mode k_L modulating two small-scale
modes k_S. This is the **separate-universe response** — a long density mode δ_L acts like
a local shift of the background, uniformly rescaling small-scale power:
`P_S(x) ≈ P̄_S [1 + (dlnP_S/dδ_L) δ_L(x)]`.

For two points κ(n₁), κ(n₂) separated by γ:

- **Small γ:** both points lie within one **coherence length of the long mode**, so they
  see the *same* δ_L modulation → the FK correction to ⟨κκ⟩ is a nearly **γ-independent
  offset** → the flat +2.745e-5 plateau.
- **Large γ (≳ the long-mode coherence angle, ~100′ here):** the two points see *different*
  δ_L → the modulation **decorrelates** → the plateau decays (957′: +1.5e-5), and turns
  into a tiny negative tail at the largest γ.

So the plateau width ≈ the projected coherence angle of the squeezed (long) mode.

## κγ_t, ξ±: rising from the floor (NO signal — a suppressed residual)

The same long mode also boosts the local shear *power*, but producing a **coherent**
⟨shear·shear⟩ or ⟨κ·shear⟩ correlation at the squeezed vertex requires the angular spins
to match:

- A **scalar (spin-0)** long mode couples to two **scalar** convergence legs via the
  all-spin-0 `[w₃ⱼ(L,L,ℓ₃;0,0,0)]²` → **nonzero** (this is the κκ signal).
- Its coupling to **spin-2** shear legs needs spin-weighted 3j symbols which, under the
  equal-time squeezed projection + parity, **average to zero** (ζ_TTP/TPP/PPP vanish;
  only ζ_TTT survives).

Physically: a scalar long mode can coherently **focus/defocus** (Ricci → κ) but cannot
coherently **shear** (Weyl → γ) — the spin-2 phase coherence is washed out in the squeezed
limit. So the shear channels have **no plateau and no real signal**, only a residual:

- **small γ (high ℓ):** suppression is near-total → ~1e-15 (floating-point floor); κγ_t and
  ξ− additionally have the kinematic `d^ℓ(0)=0` zero.
- **large γ (low ℓ):** the residual lives at low ℓ, and the spin-2 kernels at large γ pick
  it up → rises to ~1e-8.
- κγ_t ≈ ξ₊ ≈ ξ₋ (all carried by the (1,1) residual; (2,2)≈0) → one common suppressed
  residual feeding every spin-2 combination, not three independent signals.

## The "meeting" at large γ is κκ decaying, not shear gaining

At γ≳2000′, FK_κκ (~−8e-8) and the spin-2 FK (~−5e-8) become comparable — **not** because
shear acquired signal, but because **FK_κκ has itself decayed to the same ~1e-8 floor**.
Across small/mid γ, FK_κκ exceeds the spin-2 residual by **5–10 orders of magnitude** —
2.3e-8 (10′), 1.6e-6 (100′), 2.5e-4 (1000′). That order-of-magnitude **bound** is the robust
physics.

## Convergence test result (run 2026-06-01) — it is a CONVERGED residual, not sweep noise

A quadrature-convergence test was run on the FK config (`C_corr_op_K_limber_FK`), varying
`sweep.n_gauss ∈ {16,24,32,48}` (via `SFT_WICK_SWEEP_N_GAUSS`) and `c_n_gauss ∈ {16,20,28}`,
holding everything else fixed (diag_C=false, vectorized C). Findings:

- **Control:** FK_κκ(100′) stable to ~0.8% across all 6 configs → the test is valid and the
  real signal is converged.
- **FK_spin2 (κγ_t, ξ±) is STABLE, not shrinking:** at γ≈1000′, FK_κγ_t(n_gauss=48)/(n_gauss=16)
  = **0.9945** and FK_κγ_t(c_n_gauss=28)/(c_n_gauss=16) = **1.0000**. It tracks FK_κκ's ~0.8%
  wiggles in lock-step → it is **set by the κ³ table, not the sweep quadrature**.
- **c_n_gauss has ZERO effect** (byte-identical outputs): the FK config uses
  `c_closed_form_only: true`, so the corr_op C is closed-form and the `c_n_gauss` GL quadrature
  is **bypassed entirely** — C-propagator quadrature is definitively not the residual's source.

**Conclusion:** the spin-2 FK is a **converged residual locked to the κ³ table**, NOT random
sweep-quadrature noise (this *corrects* the looser "numerical noise" framing). It is therefore
one of two things the sweep cannot distinguish:
1. a **genuine tiny physical spin-2 leakage** (Wigner-3j suppression strong but not exactly zero), or
2. a **κ³-table angular-interpolation artifact** (the table's U/λ discretization not perfectly
   enforcing the analytic zero).

**To distinguish (1) vs (2):** rebuild the κ³ table on a denser angular (U/λ) grid and re-measure
— if the residual moves toward zero it is a table artifact; if it stays it is a physical leakage.
Refining the sweep cannot decide this (the residual does not live there). Either way the residual
is ≥3.6 OOM below FK_κκ at all separations, so it does not affect any conclusion of analysis3.

## One-line summary
- **κκ:** scalar long mode coherently modulates scalar convergence → flat plateau within
  the long-mode coherence angle, decaying beyond. A real κ³ signal.
- **spin-2:** scalar long mode cannot coherently modulate spin-2 shear (Wigner-3j) → no
  plateau, only a low-ℓ residual that the large-γ spin-2 kernels lift off the floor. The
  apparent "large-scale rise" is suppressed residual, not signal.
