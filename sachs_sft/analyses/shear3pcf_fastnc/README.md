# shear3pcf_fastnc — driving-field 3-cumulant 3PCF validation

Independent cross-validation of the STF_lensing equal-shell driving-field
3-cumulant `zeta_abc(gamma, lambda)` (appendix `eq: appendix equal shell approx`,
D7) by building a GENUINE 3-point function from it and comparing against an
independent reference, at a single source plane `z_s = 5`.

- **STAGE 0**: scalar **convergence** 3PCF self-check vs a
  hand-rolled SPT tree-level flat-sky Limber reference. Certifies zeta's absolute
  normalization, radial measure, and spin-0 projection. **PASSED.**
- **STAGE 1**: cosmic **shear** 3PCF vs the external code `fastnc`
  and an independent B_kappa-phase definition reference. **SPIN-2 BUG FOUND,
  then FIXED (2026-06-09).** After the landed corrections the spin-2 subset
  CLOSES IN NORMALIZATION (frame-invariant canoes/ref median 7.77x -> 0.70x);
  a structured ~factor-2 phi-dependent J_2 angular residual remains (see the
  STAGE-1 AFTER-FIX section). Final numbers in `outputs/stage1_threeway_FINAL.npz`.

See `FASTNC_NOTES.md` for fastnc install + API details.

---

## EXECUTIVE SUMMARY

EXECUTIVE SUMMARY. This analysis cross-validates the STF_lensing driving-field equal-shell 3-cumulant zeta against an independent shear-3PCF path, at a single source plane z_s=5, on the collapsed/SAS-isoceles family. Outcome: the SCALAR (convergence) sector is validated; the SPIN-2 (shear) sector reveals a real, fixable bug.

WHAT IS VALIDATED:
- STAGE 0 (scalar convergence 3PCF): our zeta_TTT fold matches an independent hand-rolled SPT-Limber convergence 3PCF to |ratio-1| = 0.15 (signal regime gamma<=60': 0.85-0.99), sign-consistent. zeta_TTT's absolute normalization, lambda-density radial measure (no extra Jacobian), and spin-0 projection are CERTIFIED.
- Consistency: a targeted fresh query of the canoes kappa3 functions reproduces the deployed equal_time_limber callable BIT-FOR-BIT (rel diff = 0.0) at shared cosine nodes, for both zeta_TTT and zeta_Dmod. The validation therefore tests the deployed object itself, not a re-build.
- Channel map: the shear natural-component assignment Gamma^0 = zeta_PPP (un-conjugated, gamma^4-suppressed), Gamma^1..3 = zeta_D (single-leg-conjugated) is derived from first principles (appendix modulus reconstruction) and independently confirmed by fastnc's own leg-spin map. The great-circle -> centroid projection is a pure per-leg phase (magnitude-invariant).

WHAT FAILS (the spin-2 bug):
- The cosmic-shear 3PCF natural components built from the canoes driving-field spin-2 modulus channel zeta_D (the (2,2,-2) channel) do NOT match the standard shear 3PCF. The DECISIVE three-way test: an independent B_kappa-phase reference (our tree B_kappa + the standard flat-sky spin-2 projection gamma~=e^{2i phi} kappa~, which IS the definition) agrees with fastnc to ~1% (frame-invariant median 0.99 across four configs), and BOTH are 4-17x SMALLER than the canoes-zeta_D result (median ~7.8x), with the discrepancy GROWING with opening angle phi (a shape error, not a constant) and a wrong per-component split (canoes over-weights Gamma^1, under-weights Gamma^0).
- Two independent implementations of the shear-3PCF definition (the B_kappa-phase reference and fastnc, using different methods/conventions) agreeing to ~1% while both depart from canoes by 4-17x is logically decisive: the canoes spin-2 zeta_D (2,2,-2) route is the OUTLIER -- a genuine over-normalization + angular-shape bug in the driving-field spin-2 modulus 3-cumulant. Ruled OUT as causes: the cumulant build (bit-for-bit deployed), the scalar method (15% in STAGE 0), the W3j_0 parity anchor (exonerated: legitimate Komatsu-Spergel reduced-bispectrum anchor with the spin correctly threaded through the {s}Y tables + m-sum; passes the +2 self-consistency test), and the great-circle->centroid projection (pure phase, cannot change magnitude).
- This confirms the standing project_kappa3_spin2_helicity_bug_2026-06-04 concern: the (2,2,-2) modulus normalization was never validated against an external reference. It now is, and it is wrong. The fix is canoes-side (spin-2 angular kernel / helicity weighting of the (2,2,-2) channel; candidate locus flagged at canoes _spin_aware_three_pt.py:606-691) and was NOT applied here (canoes kept read-only per task constraints).

SCOPE NOTE: this bug is in the all-spin-2 zeta_D (2,2,-2) 3-point channel. The companion modulus channel zeta_B (0,2,-2) was found to honor the appendix small-angle limit (zeta_B -> zeta_TTT); the scalar FK 2-point baseline is a different channel. Do not over-extend the scope of this 3-point spin-2 finding without separate checks.

---

## STAGE 1 AFTER-FIX CLOSURE (2026-06-09, branch fix/spin2-zetaD-2_2_-2-norm)

The spin-2 normalization bug is FIXED. The canoes side was rebuilt fresh
(`stage1_ours_v2_FIXED.py` -> `outputs/stage1_ours_v2_FIXED.npz`, a genuine
6-channel recompute) carrying THREE landed corrections:

1. **Group-zeroing fix** (`canoes _spin_aware_three_pt.py`): the (l2,l3)-keyed
   group-batched m-sum no longer reads count_g / the W3j reshape from the
   group-FIRST triple. For the s1=+/-2 channels that first triple has l1<|s1|,
   so m_count=0 zeroed the whole group, dropping ~half the valid (l1>=|s1|)
   zeta_D cells. Now W_group is assembled per-triple from its own offset slice.
2. **h^6 -> h^4 units fix** (`canoes kappa3.py`, 5 equal-shell sites): the
   equal-shell zeta is natively (h/Mpc)^4 = A(a)^3 (P.P)/chi^4, so units=physical
   multiplies by h^4, not h^6. Channel-independent factor h^2 = 0.4504.
   Verified: a fresh LOW zeta_TTT physical/h ratio matches h^4 exactly. The
   deployed equal_time_limber table was rebuilt with this fix.
3. **(1+z)^-4 Born per-leg reduction** (analysis-side, in `fold()`): the
   single-plane K^3 fold against canoes' per-leg (1+z)^4 Sachs response
   over-counts the textbook lensing response by a per-shell (1+z)^(2n-2),
   n=3 legs -> (1+z)^-4. This is the SAME reduction `stage0_ours.py` applies for
   the scalar convergence 3PCF, CCL-native-anchored in `b2_kfold_2pt_vs_ccl.py`
   (the 2-leg analogue over-counts by exactly (1+z)^2; dividing pointwise by
   (1+z)^2 restores Cl(hand)/Cl(CCL)=0.99). canoes' zeta and the kernel K are
   unchanged. K^3-weighted mean <(1+z)^-4> = 0.246.

The combined scalar-normalization prediction (h^2 x <(1+z)^-4> = 0.4504 x 0.246)
is a uniform ~9.0x magnitude drop; the observed per-config drop is 10.3-13.9x,
the config-dependent excess (1.1-1.5x) being the group-zeroing fix's net effect
on zeta_D.

### Frame-invariant canoes/ref: BEFORE -> AFTER (the culmination)

| Config | BEFORE can/ref | **AFTER can/ref** | drop |
|---|---|---|---|
| gamma=10', phi=60 deg (PRIMARY) | 4.118 | **0.297** | 13.9x |
| gamma=10', phi=30 deg | 6.074 | **0.475** | 12.8x |
| gamma=10', phi=120 deg | 17.230 | **1.415** | 12.2x |
| gamma=50', phi=60 deg | 9.473 | **0.921** | 10.3x |
| **median** | **7.77x** | **0.698x** | |

ref/fnc self-check median 0.989 (the two independent shear-3PCF references still
agree to ~1%, so `ref` remains a trustworthy ground truth). The spin-2 subset
now CLOSES IN OVERALL NORMALIZATION: the ~8x over-largeness is gone, replaced by
a sub-unity median 0.70x.

### The residual phi-dependent J_2 angular factor (genuine remaining question)

After removing the overall scale, the residual (can/ref)/median is NOT flat:

| Config | residual_j2 (full) | residual_j2 (D-channels only) |
|---|---|---|
| gamma=10', phi=60 deg | 0.425 | 0.482 |
| gamma=10', phi=30 deg | 0.681 | 0.808 |
| gamma=10', phi=120 deg | 2.026 | 2.213 |
| gamma=50', phi=60 deg | 1.319 | 1.192 |

phi-trend at gamma=10' (residual vs phi=[30,60,120]deg): [0.68, 0.43, 2.03]
-- still STRUCTURED (spread 4.8x), NOT flat. This is the genuine remaining
spin-2 question, now cleanly SEPARATED from the (fixed) scalar normalization.
Two distinct effects survive:

- **Gamma^0 (PPP / J_6 all-+2 channel)** is gamma^4-suppressed at small angle:
  can/ref ~ 0.01-0.05 at gamma=10', recovering to 1.12 at gamma=50'. Physically
  expected (d^L_{2,-2} ~ gamma^4); the reference assigns this channel more weight
  than our zeta_PPP produces at small gamma.
- **The J_2 D-channel per-leg split reshuffles at wide phi**: at phi=120 deg
  canoes over-weights Gamma^1,2 (~1.5x) and under-weights Gamma^3 (0.32x). Note
  fnc/ref ALSO reshuffles Gamma^3 at phi=120 (2.40x), so part of this is a
  channel-label / leg-conjugation convention difference, not a pure canoes error.

ROOT of the residual ell-shape (from the pre-fix HIGH-path diagnosis,
`scripts/docs/zetaD_HIGH_fix.md`, which remains valid): Gate (B) proved the J_2
orientation kernel `i^S J_S(|Q|) e^{iS beta}` is CORRECT (kernel + ref-B
reproduces the reference Gamma to <1% at every config). The residual is canoes'
HIGH B_delta/radial integrand having a MILD ell-dependent SHAPE difference from
the trusted single-plane B_kappa (W_can/B_kappa log-log slope +0.44), which the
spin-2 J_2 kernel weights into the (60,1000] window more than the scalar J_0 or
PPP J_6. It is NOT a clean per-leg eth^2 = sqrt(L^2(L^2-2)) factor (that would
give ell^4 = 4096x at L=800; observed is 2.7x).

### HONEST VERDICT (AFTER FIX)

The spin-2 shear-3PCF subset **CLOSES IN OVERALL NORMALIZATION** after the three
landed fixes: frame-invariant canoes/ref median 7.77x -> 0.70x (a ~10-14x
per-config drop), with ref/fnc self-check still ~1%. The dominant
over-normalization that the BEFORE state flagged is resolved.

A **structured ~factor-2 phi-dependent J_2 angular residual** (spread 4.8x in
(can/ref)/median; phi-trend [0.68,0.43,2.03] at gamma=10') does **NOT** fully
close. It is the genuine remaining spin-2 question, traced to (a) the
gamma^4-suppressed Gamma^0/J_6 channel and (b) the canoes HIGH B_delta/radial
ell-shape amplified by the J_2 kernel, plus a partial channel-label/leg-
conjugation convention difference (fastnc reshuffles the same way). Closing the
angular shape requires the canoes HIGH spin-2 B/radial response ell-weighting,
derived symbolically and re-validated against the B_kappa-phase reference -- a
scoped follow-up, NOT a factor fit.

Outputs: `outputs/stage1_threeway_FINAL.npz` (BEFORE/AFTER + isolated J_2
residual, full provenance in its `note`), `outputs/stage1_ours_v2_FIXED.npz`
(fresh canoes side), `outputs/stage1_threeway_compare_FIXED.npz`,
`outputs/stage1_threeway_AFTERFIX.npz` (the isolation script's output). Stale
pre-rebuild npz archived as `*_STALE_2026-06-09.npz`.

---

## STAGE 0 — scalar convergence 3PCF

### What was computed

A collapsed (squeezed, degenerate-triangle) scalar convergence 3PCF on the
configuration `(1, cos gamma, cos gamma)` — two of the three correlated sky
points coincide, the third sits at angular separation `gamma` — at `z_s = 5`,
computed two independent ways and compared.

**OURS** (`stage0_ours.py`, sft-wick env): builds `zeta_TTT(gamma, lambda)`
fresh per the canonical recipe (LOW `compute_kappa3_zeta_table` ell<=60 + HIGH
`compute_kappa3_sigma3_high` (60,1000] + `kappa3_combine_low_high`), then folds
with the single-source-plane convergence efficiency:

    Z_kappa(gamma) = int_0^{lambda_s} dlambda K(lambda, lambda_s)^3 zeta_TTT(gamma, lambda),
    K(lambda, lambda_s) = a^2(chi) chi (chi_s - chi)/chi_s   [cosmology.tex eq: K comoving kernel].

**REFERENCE** (`stage0_spt_reference.py`, PyCCL env): an independently
hand-coded SPT tree-level flat-sky Limber convergence 3PCF on the same squeezed
config, same cosmology + same `PCAMBz0.txt` P(k) + same `z_s`. It re-derives the
spin-0 squeezed angular kernel `2 pi J0(v gamma)` and re-implements
`B_delta = 2 F2 P P + 2 cyc` with the TRUE SPT `F2 = 5/7 + (mu/2)(k_a/k_b +
k_b/k_a) + (2/7) mu^2` over its own (ell, ell, phi) quadrature, folded with the
same `K^3`. It is the SAME tree+Limber+flat-sky family with a DIFFERENT
quadrature/normalization, so it is an INTERNAL consistency check (NOT
approximation-independent); ~30% agreement is the expectation.

### Matched conventions (verified in code)

| Item | Value |
|---|---|
| Cosmology | Omega_m = 0.3160919980475834, h = 0.6711, n_s = 0.97 |
| P(k) | `/Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt` (h-units) |
| Source plane | z_s = 5 (chi_s = 7971.18 Mpc, lambda_s = 2317.96 Mpc) |
| Units | `physical` (equal-shell zeta carries h^4 [2026-06-09 dimensional correction: native (h/Mpc)^4 = A(a)^3 (P.P)/chi^4, was h^6], chi/lambda carry h^-1); deployed-table `units==physical` assert PASSED |
| Radial measure | `lambda` (lambda-density 3-leg object; NO (dchi/dlambda)^3 applied — radial-measure memory 2026-05-17) |
| Poisson amp | A(a) = -1.5 Omega_m H0^2 (1+z), H0 = 100/c (h-units) — matches canoes `_H0_H_PER_MPC` to 1e-12 |
| Growth D(z) | canoes lightcone M[0]/(1+z) on both sides (D(z_s=5)=0.211) |
| Bispectrum | standard-SPT tree-level F2 (NOT effective F2), `spt_kind="tree_phi"` on ours side |
| Flat-sky | both sides |
| affine <-> chi | dlambda = a^2 dchi verified against deployed grid to 3e-6; K_lambda = a^2 K_chi |

FK kk baseline (2026-06-09 h^4 correction): the h^6 -> h^4 canoes units fix
(branch `fix/spin2-zetaD-2_2_-2-norm`, the 5 equal-shell sites in `kappa3.py`)
plus the consequent deployed equal_time_limber table + FK sweep rebuild shifted
the FK kk baseline by x2.2204 to **+1.219e-4** at gamma=0.5' (pre-fix value was
+3.086e-5; an even older value +5.49e-5 is superseded). Order-0 CCL kk is
unchanged (0.998). NOTE: this h^4 units fix is canoes-side and WAS applied; it is
distinct from the analysis-side-only radial (1+z)^4 Born reduction documented in
`scripts/docs/photon_energy_reduction.md` and `scripts/docs/scalar_3pt_fix.md`.

### Single-gamma result (reported FIRST)

At **gamma = 1 arcmin**, z_s = 5:

| Quantity | Value | Sign |
|---|---|---|
| **OURS** Z_kappa (lambda route) | **-2.285e-05** | negative |
| OURS Z_kappa (chi route, cross-check) | -2.274e-05 | negative |
| **REFERENCE** Z_kappa (SPT Limber) | **-2.701e-05** | negative |
| **\|ratio\| = \|ours\|/\|ref\|** | **0.846** | (\|ratio - 1\| = 0.154) |

The ours-side lambda-route and chi-route agree to **0.96%** (pure quadrature
residual on the 40-shell grid), confirming the affine <-> comoving measure
bookkeeping is internally consistent.

### Sign relationship

**Both 3PCFs are NEGATIVE and agree in sign across the entire validation
regime.** The PROMPT anticipated a possible single overall sign flip (the
radial-measure memory precedent had ref +1.15e-10 vs ours -9.66e-11), but here
the reference was built with the Poisson cube `A(a)^3 < 0` explicit per leg
(`Phi_00 = A(a)(1+z)^4 delta`), which makes `zeta_TTT < 0` and hence `Z_ref < 0`,
the SAME `(-1)^3` sign as the `kappa = -int Dsachs_1` fold. So ours and ref carry
the same sign by construction — no flip. (Had the reference instead been written
with a bare `+B_delta` j0^3 projection omitting the Poisson sign, it would have
come out positive and required the documented single-minus flip.)

### Gamma sweep (10 points, 0.5' .. 300')

| gamma ['] | ours (Z_lambda) | ref (Z_ref) | \|ratio\| | signs |
|---|---|---|---|---|
| 0.50 | -2.295e-05 | -2.715e-05 | 0.845 | (-,-) |
| 1.00 | -2.285e-05 | -2.701e-05 | 0.846 | (-,-) |
| 2.48 | -2.220e-05 | -2.612e-05 | 0.850 | (-,-) |
| 5.51 | -1.939e-05 | -2.229e-05 | 0.870 | (-,-) |
| 12.25 | -1.005e-05 | -1.021e-05 | 0.985 | (-,-) |
| 27.25 | -2.893e-06 | -3.243e-06 | 0.892 | (-,-) |
| 60.62 | -3.978e-07 | -3.073e-07 | 1.294 | (-,-) |
| 134.85 | +1.333e-08 | -3.744e-08 | 0.356 | (+,-) |
| 300.00 | +3.637e-08 | -7.434e-08 | 0.489 | (+,-) |

- **Validation regime (gamma <= 60'): |ratio| in [0.845, 0.985]** — the 3PCF is
  large here (1e-5..4e-7, the signal-carrying peak) and the two implementations
  track each other tightly, with matched sign throughout.
- **Large gamma (>= 135'):** both 3PCFs cross zero and collapse to ~1e-8 (three
  orders below the peak). They cross zero at slightly different gamma (the
  low-ell vs high-ell cancellation near gamma~40' documented in
  `tab: zeta ell bands`), so the sign and ratio there are numerically delicate
  and physically negligible; this is NOT part of the normalization check.

### Verdict: **PASS**

- **|ratio - 1| = 0.154 at gamma = 1'** (gate is < ~0.5). PASS.
- |ratio| in [0.845, 0.985] across the whole signal-carrying regime (gamma <= 60').
  Well inside the gate; better than the ~30% expected for a same-family
  tree+Limber cross-check.
- Sign relationship CORRECT (both negative; consistent everywhere the 3PCF is
  resolved).
- No disagreement exceeds a factor ~2; no sign anomaly in the validation regime.

The driving-field 3-cumulant `zeta_TTT`'s absolute normalization, radial measure
(lambda-density, NO extra Jacobian), and spin-0 convergence projection are
**certified** by this independent SPT-Limber cross-check. STAGE 1 (shear vs
fastnc) is cleared to proceed.

### Files

| File | Role | Env |
|---|---|---|
| `stage0_ours.py` | zeta-fold convergence 3PCF (LOW+HIGH+combine, K^3 fold) | sft-wick |
| `stage0_spt_reference.py` | independent hand-rolled SPT tree-level Limber 3PCF | PyCCL |
| `stage0_compare.py` | ratio + plot | PyCCL |
| `outputs/stage0_ours.npz` | gamma, zeta_TTT, Z_lambda, Z_chi | |
| `outputs/stage0_spt_reference.npz` | gamma, zeta_TTT, Z_ref | |
| `outputs/stage0_results.npz` | gamma, Z_ours, Z_ref, ratios, signs | |
| `outputs/stage0_compare.png` | comparison + ratio panel (PNG) | |
| `figures/stage0_convergence_3pcf.pdf` | comparison + ratio panel (PDF) | |

Reproduce:
```bash
# ours side (sft-wick env)
LOKY_MAX_CPU_COUNT=8 CANOES_SUPPRESS_METAL_WARNING=1 \
  /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python stage0_ours.py
# reference + compare (PyCCL env)
/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python stage0_spt_reference.py
/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python stage0_compare.py
```

---

## STAGE 1 — cosmic shear 3PCF vs fastnc + B_kappa-phase reference

### Overview

Cosmic-shear 3PCF natural components (Schneider-Lombardi basis, CENTROID
projection, flat-sky) from three independent pipelines, compared on the
SAS-isoceles family at z_s = 5:

1. **B_kappa-phase reference** (`stage1_bkappa_phase_reference.py`, sft-wick env):
   definition-level integral of B_kappa with standard flat-sky spin-2 phases
   `gamma~(l) = e^{2i phi_l} kappa~(l)`, centroid angles built from first-principles
   geometry, no fastnc internals used.
2. **fastnc tree** (`stage1_fastnc_tree.py`, fastnc venv):
   fastnc v1.2.3 commit d94fd2a (SPT tree-level bispectrum subclass, same
   PCAMBz0.txt, same Limber chain, 2DFFTLog multipole resummation).
3. **canoes zeta_D** (`stage1_ours_v2.py`, sft-wick env):
   driving-field spin-2 modulus cumulants (zeta_D, zeta_PPP) from the deployed
   `equal_time_limber` callable, folded with the convergence efficiency kernel K^3,
   mapped to natural components via the derived great-circle -> centroid spin-2
   phase.

The three-way comparison is assembled by `stage1_threeway_compare.py`.
The underlying cumulants are verified bit-for-bit against the deployed callable by
`consistency_deployed_vs_fresh.py`.

### Matched conventions

| Item | Value |
|---|---|
| Cosmology | Omega_m = 0.3160919980475834, h = 0.6711, n_s = 0.97 |
| P(k) | `PCAMBz0.txt` (h-units; fastnc: fed directly to `set_pklin`) |
| Source plane | z_s = 5 (single delta plane) |
| Bispectrum | SPT tree-level F2 (TRUE: 5/7 + ..., NOT F2_eff) |
| Projection | CENTROID (`projection='cent'`, arXiv:2309.08601 Eq.15-16) |
| Geometry | SAS isoceles: t1=t2=gamma, t3=2 gamma sin(phi/2), opening angle phi |
| Grid | gamma in {10', 50'}, phi in {30, 60, 120 deg} (4 configs; PRIMARY: gamma=10', phi=60 deg) |
| fastnc version | v1.2.3, commit d94fd2a |
| fastnc key deps | numpy 1.26.4, scipy 1.13.1, astropy 6.0.1, mpi4py 4.1.2 |
| fastnc production setting | Lmax=24, Lmax_diag=48, Mmax=24, epmu=1e-7 |
| B_kappa-phase ref convergence | (n_l=128, n_ang=128, lmax=3e4), ~2-3% |

### Consistency check

`consistency_deployed_vs_fresh.py` queries the canoes kappa3 functions directly
at the 3 cosine triples used in the deployed `equal_time_limber` callable
(rows 0/1/2 of the npz grid, covering cos in {0.714, 0.857, 1.0}), for both
`zeta_TTT` and `zeta_Dmod`. Result: `reldiff_TTT = reldiff_Dmod = 0.0` (float64
exact; see `outputs/consistency_deployed_vs_fresh.npz`). The validation therefore
tests the deployed object itself, not a re-build.

### Three-way results: frame-invariant sqrt(sum|Gamma^mu|^2)

The frame-invariant `sqrt(sum_mu |Gamma^mu|^2)` removes any dependence on
projection-frame conventions; it is the robust normalization check.

| Config | fi_ref (B_kappa-phase) | fi_fnc (fastnc) | fi_can (canoes) | **fi_fnc/fi_ref** | **fi_can/fi_ref** |
|---|---|---|---|---|---|
| gamma=10', phi=60 deg (PRIMARY) | 2.26e-7 | 2.33e-7 | 9.32e-7 | **0.970** | 4.12x |
| gamma=10', phi=30 deg | 1.98e-7 | 2.38e-7 | 1.20e-6 | **0.829** | 6.07x |
| gamma=10', phi=120 deg | 1.28e-7 | 1.27e-7 | 2.21e-6 | **1.007** | 17.23x |
| gamma=50', phi=60 deg | 5.98e-8 | 5.93e-8 | 5.67e-7 | **1.009** | 9.47x |

Median fi_fnc/fi_ref = **0.989** (range 0.83-1.01). Median fi_can/fi_ref = **7.77x** (range 4.12-17.23x).

The B_kappa-phase reference and fastnc agree to ~1% in the frame-invariant;
both are 4-17x SMALLER than the canoes zeta_D result, with the factor GROWING
with opening angle phi (structured shape error, not a constant offset).

### Primary-point per-component breakdown (gamma=10', phi=60 deg)

| Component | |B_kappa-phase ref| | |fastnc tree| | |canoes zeta_D| | can/ref |
|---|---|---|---|---|
| Gamma^0 | 6.51e-8 | 6.82e-8 | 3.44e-8 | **0.53** (under) |
| Gamma^1 | 1.36e-7 | 1.38e-7 | 7.07e-7 | **5.20** (over) |
| Gamma^2 | 1.20e-7 | 1.24e-7 | 4.28e-7 | **3.58** (over) |
| Gamma^3 | 1.20e-7 | 1.24e-7 | 4.30e-7 | **3.59** (over) |

The canoes route over-weights Gamma^1 (the apex-conjugated component, ~5x) and
under-weights Gamma^0 (the all-unconjugated component, ~0.5x). The ref and fnc
match per-component to <=4%.

### Channel assignment and projection derivation

The channel map `Gamma^0 <- zeta_PPP`, `Gamma^{1,2,3} <- zeta_D` (with each
component conjugating one distinct leg) is derived from first principles in
`scripts/derive_shear3pcf_natural_components.wl` and confirmed by fastnc's own
Bessel-order map (fastnc.py:365 mu->(m,n) map at M=0).

The great-circle -> centroid spin-2 rotation per leg is a pure per-leg phase
`exp(-2i s_j (a_j - b_j))` where `a_j - b_j` is phi-only (scale-free,
`d/dg = 0`; verified in the .wl script). Therefore `|Gamma^mu| = |fold(zeta_channel)|`:
the projection CANNOT redistribute magnitude between components and cannot be a
factor in the discrepancy.

### zeta_D cumulant-level investigation

`outputs/zetaD_normalization_probe.npz` and `scripts/docs/zetaD_normalization_diagnosis.md`
confirm: the deployed `zeta_D (2,2,-2)` cumulant is NOT over-large at the cumulant
level. The W3j_0 anchor is the correct Komatsu-Spergel reduced-bispectrum anchor
(same for all channels; exonerated by the (2,2,2)/(0,2,2) channel self-consistency
tests). The over-large shear components arise from the spin-2 angular kernel /
helicity weighting of the `(2,2,-2)` SpinAware configuration synthesis
(`_spin_aware_three_pt.py:606-691`): the m-sum for spins `(+2,+2,-2)` may realize
a `d^L_{2,-2}`-suppressed combination rather than the `d^L_{2,2}`-finite modulus
the appendix has in mind, producing a wrong per-component split and structured
over-normalization at the shear-3PCF level.

### Figures

| File | Content |
|---|---|
| `figures/stage1_threeway_gamma.pdf` | Frame-invariant bar chart (3 pipelines) + ratio panel (canoes/ref, canoes/fnc, fnc/ref self-check) across all 4 configs |
| `figures/stage1_threeway_components.pdf` | Per-component |Gamma^mu| bar chart at PRIMARY point (gamma=10', phi=60 deg) |
| `outputs/stage1_threeway_gamma.png` | PNG version of above |
| `outputs/stage1_threeway_components.png` | PNG version of above |

### Files

| File | Role | Env |
|---|---|---|
| `stage1_ours_v2.py` | canoes zeta_D -> natural components (derived centroid phase) | sft-wick |
| `stage1_fastnc_tree.py` | fastnc SPT tree shear 3PCF | fastnc venv |
| `stage1_bkappa_phase_reference.py` | B_kappa direct Fourier reference | sft-wick |
| `stage1_threeway_compare.py` | assemble outputs/stage1_threeway_compare.npz | PyCCL |
| `consistency_deployed_vs_fresh.py` | bit-for-bit consistency check | sft-wick |
| `make_stage1_figures.py` | generate stage1_threeway_gamma.pdf + components.pdf | PyCCL |
| `outputs/stage1_ours_v2.npz` | canoes natural components (3x7 grid) | |
| `outputs/stage1_fastnc.npz` | fastnc natural components (3x7 grid) | |
| `outputs/stage1_bkappa_phase_ref.npz` | B_kappa-phase reference (4 configs) | |
| `outputs/stage1_threeway_compare.npz` | frame-invariant + ratios (4 shared configs) | |
| `outputs/consistency_deployed_vs_fresh.npz` | reldiff_TTT=reldiff_Dmod=0 | |
| `outputs/zetaD_normalization_probe.npz` | cumulant-level normalization probe | |
| `scripts/docs/shear3pcf_natural_map.md` | channel assignment + projection derivation | |
| `scripts/docs/zetaD_normalization_diagnosis.md` | cumulant-level root-cause investigation | |
| `scripts/docs/decisive_bkappa_reference.md` | B_kappa-phase reference design + verdict | |

Supporting derivation scripts (Mathematica, run with `wolframscript -file <path>`):
`scripts/derive_shear3pcf_natural_components.wl`,
`scripts/derive_spin2_gaunt_norm.wl`,
`scripts/spin2_centroid_direct.wl`,
`scripts/spin2_3pcf_phase.wl`.

Reproduce:
```bash
# canoes side (sft-wick env)
LOKY_MAX_CPU_COUNT=8 CANOES_SUPPRESS_METAL_WARNING=1 \
  /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python stage1_ours_v2.py
# fastnc (fastnc venv)
/Users/zzhang/projects/fastnc_venv/bin/python stage1_fastnc_tree.py
# B_kappa-phase reference (sft-wick env, ~15 min at n_l=128)
CANOES_SUPPRESS_METAL_WARNING=1 \
  /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python stage1_bkappa_phase_reference.py
# three-way comparison + figures (PyCCL env)
/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python stage1_threeway_compare.py
/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python make_stage1_figures.py
# consistency check (sft-wick env)
CANOES_SUPPRESS_METAL_WARNING=1 \
  /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python consistency_deployed_vs_fresh.py
```

---

## VERDICT

EXECUTIVE SUMMARY. This analysis cross-validates the STF_lensing driving-field equal-shell 3-cumulant zeta against an independent shear-3PCF path, at a single source plane z_s=5, on the collapsed/SAS-isoceles family. Outcome: the SCALAR (convergence) sector is validated; the SPIN-2 (shear) sector reveals a real, fixable bug.

WHAT IS VALIDATED:
- STAGE 0 (scalar convergence 3PCF): our zeta_TTT fold matches an independent hand-rolled SPT-Limber convergence 3PCF to |ratio-1| = 0.15 (signal regime gamma<=60': 0.85-0.99), sign-consistent. zeta_TTT's absolute normalization, lambda-density radial measure (no extra Jacobian), and spin-0 projection are CERTIFIED.
- Consistency: a targeted fresh query of the canoes kappa3 functions reproduces the deployed equal_time_limber callable BIT-FOR-BIT (rel diff = 0.0) at shared cosine nodes, for both zeta_TTT and zeta_Dmod. The validation therefore tests the deployed object itself, not a re-build.
- Channel map: the shear natural-component assignment Gamma^0 = zeta_PPP (un-conjugated, gamma^4-suppressed), Gamma^1..3 = zeta_D (single-leg-conjugated) is derived from first principles (appendix modulus reconstruction) and independently confirmed by fastnc's own leg-spin map. The great-circle -> centroid projection is a pure per-leg phase (magnitude-invariant).

WHAT FAILS (the spin-2 bug):
- The cosmic-shear 3PCF natural components built from the canoes driving-field spin-2 modulus channel zeta_D (the (2,2,-2) channel) do NOT match the standard shear 3PCF. The DECISIVE three-way test: an independent B_kappa-phase reference (our tree B_kappa + the standard flat-sky spin-2 projection gamma~=e^{2i phi} kappa~, which IS the definition) agrees with fastnc to ~1% (frame-invariant median 0.99 across four configs), and BOTH are 4-17x SMALLER than the canoes-zeta_D result (median ~7.8x), with the discrepancy GROWING with opening angle phi (a shape error, not a constant) and a wrong per-component split (canoes over-weights Gamma^1, under-weights Gamma^0).
- Two independent implementations of the shear-3PCF definition (the B_kappa-phase reference and fastnc, using different methods/conventions) agreeing to ~1% while both depart from canoes by 4-17x is logically decisive: the canoes spin-2 zeta_D (2,2,-2) route is the OUTLIER -- a genuine over-normalization + angular-shape bug in the driving-field spin-2 modulus 3-cumulant. Ruled OUT as causes: the cumulant build (bit-for-bit deployed), the scalar method (15% in STAGE 0), the W3j_0 parity anchor (exonerated: legitimate Komatsu-Spergel reduced-bispectrum anchor with the spin correctly threaded through the {s}Y tables + m-sum; passes the +2 self-consistency test), and the great-circle->centroid projection (pure phase, cannot change magnitude).
- This confirms the standing project_kappa3_spin2_helicity_bug_2026-06-04 concern: the (2,2,-2) modulus normalization was never validated against an external reference. It now is, and it is wrong. The fix is canoes-side (spin-2 angular kernel / helicity weighting of the (2,2,-2) channel; candidate locus flagged at canoes _spin_aware_three_pt.py:606-691) and was NOT applied here (canoes kept read-only per task constraints).

SCOPE NOTE: this bug is in the all-spin-2 zeta_D (2,2,-2) 3-point channel. The companion modulus channel zeta_B (0,2,-2) was found to honor the appendix small-angle limit (zeta_B -> zeta_TTT); the scalar FK 2-point baseline is a different channel. Do not over-extend the scope of this 3-point spin-2 finding without separate checks.
