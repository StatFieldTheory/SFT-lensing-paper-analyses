# ℓ-band decomposition of the collapsed-3PCF cumulant ζ

**Date:** 2026-06-04 (numbers refreshed 2026-06-13 after the h⁶→h⁴ and
radial-measure (1+z)⁻⁴ fixes; values dropped by ≈(1+z)⁴·h⁻² ≈ 9.1×10² but the
ℓ-shape and crossover are unchanged). **What:** how much each ℓ region contributes to the
equal-shell driving-field cumulant ζ at different angular scales γ, on the
squeezed/collapsed geometry `(cosγ12,cosγ23,cosγ31)=(1,cosγ,cosγ)` that the
lensing 2PCF samples.

**Headline:** at the arcminute scales relevant to the κκ 2PCF, **small-scale
(high-ℓ) power dominates ζ by ~2 orders of magnitude** and sets it negative;
there is a clean sign/dominance crossover near γ≈42′, beyond which the γ-flat
low-ℓ (ℓ≲30) contribution takes over (positive). This is the ζ-level companion
to the FK result (the post-fix κκ(0.5′)=+4.29e-6 baseline is high-ℓ-dominated).

## Method

ζ is `LOW (exact Wigner-3j, ℓ≤60) + HIGH (flat-sky Born-Limber, 60<ℓ≤1000)`. The
two branches are computed DIRECTLY (no callable/interpolation) on the squeezed
family at 28 log-spaced γ (0.5′–5000′) × the 16 λ-shells, then split into bands:
- LOW bands [2,30], (30,60] by differencing cumulative `compute_kappa3_zeta_table`
  (ell_max=30, 60);
- HIGH bands (60,125], (125,250], (250,500], (500,1000] DIRECTLY from
  `compute_kappa3_sigma3_high(ell_cut=lo, ell_high_max=hi)` (disjoint windows).

Script: `ell_band_decomp.py`. Figure: `figures/zeta_ell_band_decomposition_draft.{pdf,png}`.

**Self-consistency gate:** `Σ bands ≈ full(LOW≤60 + HIGH(60,1000])`, median rel
diff **7e-4** (TTT; 1e-4–1e-3 all channels). The large *max* rel diff is the
divide-by-near-zero artifact at the sign-crossing γ where the total ζ≈0 — not a
decomposition error (the per-window GL quadrature is in fact more converged than
the single full-range build).

## Per-band ζ_TTT at the source shell (z = 5.7), signed [Mpc-physical, h⁴ units]

| γ ['] | LOW [2,30] | LOW (30,60] | HIGH (60,125] | HIGH (125,250] | HIGH (250,500] | HIGH (500,1000] | total | high-ℓ share |
|---|---|---|---|---|---|---|---|---|
| 0.50 | +1.35e-16 | −1.92e-18 | −3.56e-17 | −5.87e-16 | −4.68e-15 | **−2.18e-14** | −2.70e-14 | 100% |
| 5.45 | +1.35e-16 | −1.91e-18 | −3.55e-17 | −5.77e-16 | −4.40e-15 | **−1.73e-14** | −2.22e-14 | 99% |
| 21.31 | +1.35e-16 | −1.89e-18 | −3.32e-17 | −4.52e-16 | −1.58e-15 | +2.67e-16 | −1.66e-15 | 93% |
| 42.16 | +1.35e-16 | −1.82e-18 | −2.68e-17 | −1.80e-16 | +1.87e-16 | −6.85e-17 | +4.61e-17 | 40% |
| 117.31 | +1.35e-16 | −1.26e-18 | +1.41e-18 | −6.16e-18 | +8.28e-17 | −3.08e-17 | +1.82e-16 | 26% |
| 326.43 | +1.35e-16 | −3.02e-19 | +1.28e-18 | +4.06e-18 | −7.18e-18 | +1.09e-17 | +1.44e-16 | 6% |

"high-ℓ share" = `|Σ_{ℓ>60} ζ| / (|Σ_{ℓ>60} ζ| + |Σ_{ℓ≤60} ζ|)`.

## Physical interpretation

This is the squeezed-bispectrum coupling: at the collapsed geometry the two
coincident legs carry multipole L and the transfer across γ rides on ℓ₃; a
long-wavelength γ-separation mode couples to two short-wavelength (high-L) modes,
so high-ℓ small-scale power feeds ζ at all γ. Concretely:
- The **highest band (ℓ∈(500,1000]) is the single largest contributor at γ≲20′**,
  ~160× the entire ℓ≤60 contribution; the bands stack monotonically with ℓ.
- Each high-ℓ band's contribution **decays with γ** (and the higher-ℓ band decays
  faster, its angular transfer P_{ℓ₃}(cosγ) suppressing at γ∼1/ℓ), so the high-ℓ
  dominance is strongest at small γ.
- The **low-ℓ [2,30] band is γ-flat and positive** (≈+1.2e-13); it only wins once
  the high-ℓ bands have decohered, beyond the **crossover γ≈42′**.
- The total ζ_TTT therefore **flips sign**: negative (high-ℓ-dominated) at arcmin
  scales, positive (low-ℓ-dominated) at degree scales.

So "small-scale power contributes significantly to ζ" holds — and at the
weak-lensing arcmin scales it is not just significant but **dominant**. The
caveat (honest): the dominance is scale-dependent, with the crossover at ~tens of
arcmin; the figure focuses on γ≲350′ (the regime relevant to the κκ 2PCF).
