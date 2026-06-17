# Why the FK correction breaks the κ–γ "same information" narrative

**Context.** analysis3 shows the Order-2 NLO decomposition O0 / FF / FK across the four
2-point observables. The striking feature: the **FK** (one F + one κ³ K vertex) channel
contributes a few percent to ⟨κκ⟩ but is suppressed by **5–6 orders of magnitude** in the
shear-shear (ξ±) and convergence-shear (κγ_t) channels. At first glance this seems to
contradict the textbook lensing statement that *κ and γ carry the same information*. This
note explains why the asymmetry is **expected physics**, not a pipeline artifact.

Equation references are to the STF_lensing draft (`sections/*.tex`).

---

## 1. The traditional "κ = γ" narrative and its assumptions

In standard weak lensing, both convergence and shear are second derivatives of a **single**
lensing potential ψ:

```
κ   = ½ ∇²ψ
γ₊  = ½ (∂₁² − ∂₂²) ψ ,   γ× = ∂₁∂₂ ψ
```

In Fourier space this gives `κ̃(ℓ) = −½ ℓ² ψ̃(ℓ)` and `γ̃(ℓ) = −½ ℓ² ψ̃(ℓ) e^{2iφ_ℓ}`, so

```
|κ̃(ℓ)| = |γ̃(ℓ)|   and   ⟨κκ⟩ = ⟨γγ⟩  (E-mode),
```

differing only by the spin-2 phase `e^{2iφ_ℓ}`. **This equivalence holds because κ and γ
are two projections of the *same* (Gaussian) scalar field ψ.** Its hidden assumptions are:
**(i)** linear theory, **(ii)** the Born approximation, **(iii)** a single scalar potential,
and (for the canonical Limber form) **(iv)** Limber whitening (`insights.tex:150–167`).

---

## 2. In the Sachs formalism, κ and γ have *different* geometric sources

The Sachs picture does **not** start from one potential. The two driving fields are
geometrically distinct (`sachs_dynamics.tex:55–85`):

- **Φ00 = −½ R_μν kᵘkᵛ = −4πG T_μν kᵘkᵛ** — the Ricci focusing scalar (spin-0), sourced
  **directly and locally by the stress-energy (matter density) along the line of sight**.
  It drives the expansion: `ρ̇ = ρ² + σσ̄ + Φ00` → convergence **κ**.
- **Ψ0 = −C_αβγδ kᵅZᵝkᵞZᵟ** — the Weyl shear scalar (spin-2, complex), sourced by the
  **integrated tidal / Weyl field** (a transverse Hessian). It drives the shear:
  `σ̇ = (ρ+ρ̄)σ + Ψ0` → shear **γ**.

The observables are signed line-of-sight integrals of the Sachs scalars
(`cosmology.tex:833, 849`):

```
κ        = −∫₀^λ X₁ dλ′                          (X₁ from Φ00 / Ricci)
γ₊ ± iγ× = −∫₀^λ (X₂ ± iX₃) dλ′                   (X₂,₃ from Ψ0 / Weyl)
γ_t      = −γ₊  (tangential rotation into the great-circle frame)
```

**At linear order in the scalar (Φ=Ψ) sector**, Φ00 ∝ ∇²(Φ+Ψ) and Ψ0 ∝ ẐⁱẐʲ∂ᵢ∂ⱼ(Φ+Ψ) are
exactly the spin-0 and spin-2 projections of the *same* Hessian of (Φ+Ψ) — so they reduce
to the textbook single-potential picture and the κ=γ E-mode equivalence is **recovered**.
The draft makes this explicit and shows ⟨κγ_t⟩ is the parity-even TP=TE cross
(`cosmology.tex:735–740, 872–873`): nonzero at Born (E-mode), while ⟨κγ×⟩ vanishes by parity.

> This is why ⟨κγ_t⟩ **must be nonzero at Born** — a zero would be the anomaly. (Indeed, a
> `propagators.diag_C` default-True bug was silently zeroing the off-diagonal C and hence
> κγ_t; fixed by setting `diag_C: false`. See the run READMEs.)

---

## 3. Why FK (κ³) breaks the equivalence: the Wigner-3j selection rule

The equivalence in §1–§2 is a **linear-theory + Born** statement. The FK correction is
neither: it injects a **three-point cumulant** (source non-Gaussianity) via the κ³ K vertex,
at post-Born Order-2. At this order Φ00 (local matter density) and Ψ0 (integrated tidal) are
**no longer two projections of one Gaussian potential** — the non-Gaussian source feeds them
differently.

Concretely, the windowed-squeezed κ³ angular projection is a parity-even, triangle-restricted
**Gaunt sum over Wigner-3j symbols** `[w₃ⱼ(L,L,ℓ₃; 0,0,0)]²` (`appendix.tex:682–704`). For
channels carrying a spin-2 (Ψ0/shear) leg, the spin-weighted phase factors **average toward
zero** at the equal-time K-vertex projection; only the **all-scalar TTT** (Φ00³, spin-0
Ricci) channel survives. In the current snapshot **only ζ_TTT is nonzero**; ζ_TTP, ζ_TPP,
ζ_PPP vanish (`appendix.tex:703–704`).

Physically: a **squeezed matter-density bispectrum projects straight onto the Ricci /
convergence (spin-0)** channel, because Φ00 *is* the local matter density. Building a nonzero
equal-time **spin-2** three-cumulant requires phase coherence that the Wigner-3j geometry
washes out in the squeezed limit. So the κ³ leakage is a Ricci (density-focusing) effect
confined to convergence.

The draft states and numerically confirms the resulting asymmetry
(`insights.tex:296–312`, `discussion.tex:115–118, 266–275`):

- FK on **ξ_κ (κκ)**: plateaus at a few ×10⁻⁵, comparable to the FF post-Born amplitude.
- FK on **ξ±, ξ_{κγ_t}**: several orders below FF and sub-percent — `|ΔS₈^FK|/|ΔS₈^FF|` =
  1.0×10⁻⁶ on ξ₊ and 4.5×10⁻⁵ on ξ₋, "confirming the selection rule to five and six orders
  of magnitude" (residual = numerical noise of the 𝒦⁽³⁾ table).

---

## 4. Numerical evidence from this analysis (sachs_sft)

At the canonical cell (γ = 0.5′, z_s = 5), corr_op C propagator, dt=2:

| channel | O0 (Born) | FK (Order-2) | FK / O0 |
|---|---|---|---|
| κκ          | +8.44e-4 | +2.74e-5  | **3.25 %** |
| κγ_t        | +5.67e-6 (E-mode) | ~1e-15 (Wigner-3j suppressed) | ~0 |
| γγ (ξ₊)     | ~8.4e-4  | ~1e-15 (Wigner-3j suppressed) | ~0 |

The measured FK κγ_t ≈ 1e-15 (vs FK κκ = 2.74e-5) is a direct, end-to-end confirmation of
the Wigner-3j suppression of FK in the spin-2-bearing channels.

For the **angular-scale (γ) shape** of these contributions — why FK_κκ is a flat plateau
while the spin-2 FK rises monotonically toward large γ — see the companion note
[`FK_GAMMA_TRENDS_ACROSS_PANELS.md`](FK_GAMMA_TRENDS_ACROSS_PANELS.md).

---

## 5. Bottom line

"κ and γ carry the same information" is true **only at linear order**, where both are
projections of one Gaussian potential. The FK correction is a **non-Gaussian, post-Born**
effect — precisely the regime in which that equivalence is *designed* to break. The κ³ source
non-Gaussianity is a **Ricci (density-focusing)** effect that the Wigner-3j geometry confines
to the **spin-0 convergence** channel; the **spin-2 shear**, being a transverse-tidal (Weyl)
object, does not receive the squeezed-κ³ leakage.

**Consequence:** at this order κ and γ carry *different* information — **convergence (κ)
becomes a more sensitive probe of source non-Gaussianity**, while **shear-shear (γγ) remains
a clean post-Born probe** (only the FF channel contributes at percent level). This is the
draft's headline physics result, not a contradiction of weak-lensing theory.

---

### A note on signs (Sachs-scalar → observable)
The sft-wick output legs are the **raw** Sachs-scalar LoS integrals ∫Xᵢ. Converting to the
observables (κ=−∫X₁, γ₊=−∫X₂, γ_t=−γ₊) flips the sign of **odd-γ_t** 2-point combinations
only: ⟨κγ_t⟩ = −(raw (0,1)); ⟨κκ⟩, ξ± are unchanged (sign² = +1). This is why κκ/γγ looked
correct throughout while κγ_t needed both the diag_C fix (to be nonzero) and the sign flip
(to be the physical, positive E-mode cross).
