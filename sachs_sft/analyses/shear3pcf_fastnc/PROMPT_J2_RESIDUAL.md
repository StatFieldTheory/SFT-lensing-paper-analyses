# Handoff prompt — the J_2 angular residual in the cosmic-shear 3PCF

> **RESOLVED 2026-06-10.** Root cause = analysis-side natural-component
> slot+phase convention bug in `stage1_ours_v2_JACFIX.py` (apex-conjugated
> channel mis-slotted Gamma^3 -> Gamma^1; `_amb` great-circle phase instead of
> centroid phase). NOT canoes-core, NOT physics; paper observables unchanged
> (norm-preserving relabel). Decisive check: SAS isoceles mirror symmetry
> (Gamma^1 real, Gamma^2=conj(Gamma^3)) restored from |Im Gamma^1|/|Gamma^1|
> 0.3-0.97 -> 0 exactly. Boundary pytest
> `tests/test_natural_component_mirror_symmetry.py` (25 pass / 6 fail-without).
> See memory `project_j2_residual_closed_2026-06-10.md`. Optional follow-up only:
> the gamma-dependent apex SCALE + phi=90 outlier from the LOW covgrid16/ell_max=60
> truncation (OOM-level, not a shape residual). The investigation notes below are
> kept for provenance.

**Created:** 2026-06-09 (end of the shear3pcf_fastnc validation session)
**Branch:** `fix/spin2-zetaD-2_2_-2-norm` (canoes) — two fixes already committed
**Read first:** memory `project_shear3pcf_fastnc_validation_2026-06-09.md`
(the full cascade), and `project_kappa3_h6_overnorm_codebug_2026-06-09.md`.

---

## TL;DR for the next session

The shear-3PCF-vs-`fastnc` validation closed three confounds. Two were real
canoes bugs (now committed), one was a validation-script artifact. **One genuine
spin-2 question remains, cleanly isolated: a structured, phi-dependent, ~factor-2
angular residual in the natural-component cosmic-shear 3PCF that survives after
the overall normalization closes.** Your job is to find its root cause.

DO NOT re-litigate the three resolved items below. They are settled and committed
/ archived. Start from the open question.

---

## What is already RESOLVED (do not redo)

1. **h^6 -> h^4 units bug (committed: canoes `e19b3e4`).**
   The equal-shell 3-cumulant density is natively `(h/Mpc)^4 = A(a)^3 (P.P)/chi^4`,
   so h-native -> physical multiplies by `h^4`, not `h^6`. The surplus `h^2` broke
   h-invariance of the dimensionless 3PCF. Fixed at all 5 equal-shell sites in
   `canoes/src/canoes/sachs/kappa3.py`; cross-shell sigma3 sites correctly stay
   `h^4`. h-invariance now 1.000000.
   - Production impact: FK kappa-kappa baseline `+3.086e-5` -> `+1.219e-4`
     (x2.2204 = h^2). Order-0 (CCL-validated) kappa-kappa is byte-identical.
   - Deployed table rebuilt with h^4:
     `callables/kappa3_vertex/equal_time_limber/equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz`
     (old h^6 archived as `..._PREEH4_h6_2026-06-09.npz`). FK 2PCF sweep rebuilt
     under `sftwick_outputs/2PCF/C_corr_op_K_limber_FK/` (old `_PREEH4_h6`).

2. **LOW spin-2 group-batch cell-zeroing (committed: canoes `e9a34aa`).**
   In `_spin_aware_three_pt.py::_build_kernel_numpy`, count/m2_start/W_group were
   read from the group-FIRST triple; when `l_1 < |s_1|` masked it, whole groups
   were zeroed. Bit-identical for `s_1 == 0`; only the spin-2 path was affected.

3. **(1+z)^4 over-count = ANALYSIS-FOLD-ONLY, NOT production.**
   The production sft-wick pipeline (bare cross-shell R-propagator
   `R(t,l')=Theta(t-l')[D(l')/D(t)]^2`, D=a*chi, plus corr_op C) distributes the
   per-leg `(1+z)` correctly: O0/CCL = 0.855, the clean control. ONLY the
   hand-rolled single-plane K-fold (`stage0_ours.py`, `stage1_ours_v2_FIXED.py`)
   over-counted, so a per-shell Born reduction `(1+z)^-(2n-2) = /(1+z)^4` for
   n_legs=3 was applied **to the analysis scripts only**. canoes production was
   NOT changed for this. Scalar 3-pt then closes 22 -> 0.887 vs PyCCL.

The "(1+z)^(2n-2)" over-count law (per-leg Sachs response gives (1+z)^3 vs
textbook W~(1+z); single affine measure refunds (1+z)^2 once): 2-pt -> (1+z)^2,
3-pt -> (1+z)^4. See `scripts/docs/photon_energy_reduction.md`.

---

## THE OPEN QUESTION — the J_2 angular residual

After the h^4 fix + the analysis-side Born reduction, the frame-invariant
`canoes/ref` ratio of the natural-component shear 3PCF dropped:

| (gamma', phi deg) | BEFORE | AFTER |
|---|---|---|
| 10', 60  | 4.12  | 0.30 |
| 10', 120 | 17.23 | 1.42 |
| 50', 60  | 9.47  | 0.92 |
| **median** | **7.77x** | **0.70x** |

The OVERALL SCALE closes (matches the prediction `h^2 * <(1+z)^-4> = 0.45*0.246
= 9.0x`; ref/fnc self-check 0.989). **But the residual is NOT flat** — after
removing the median scale, `(can/ref)/median` is structured and phi-dependent.

**Symptom (the thing to explain):** at gamma=10' the phi-trend of the residual
is approximately `[0.68 @ 30deg, 0.43 @ 60deg, 2.03 @ 120deg]`. A ~factor-2
swing with opening angle phi that does not cancel against the scale.

### What is already known about the residual

- **It is radial, not angular-kernel.** The J_2 orientation kernel itself was
  PROVEN correct: `2*pi * i^S * J_S(|Q|) * exp(i*S*phi_Q)` with `S = s1+s2+s3`
  (`_kappa3_limber_alpha_kernel`, canoes `kappa3.py` ~line 1747). Do not suspect
  the angular integral.
- **It tracks the canoes HIGH B_delta / radial ell-shape.** The decomposition
  `outputs/Wcan_vs_Bkappa_decomp.npz` shows the canoes radial kernel `W_can`
  vs the trusted B_kappa phase reference has a log-log slope of **+0.44** in ell
  — i.e. canoes' tree-level radial/bispectrum kernel has a slightly different
  ell-dependence than the SPT B_kappa reference, and J_2 amplifies that
  differential into the phi-trend.
- **Part of it is a channel-label / leg-conjugation convention.** The natural
  components map as `Gamma^0 = zeta_PPP` (un-conjugated, gamma^4-suppressed) and
  `Gamma^{1,2,3} = zeta_D` (single-leg-conjugated modulus). The `fnc/ref`
  reshuffle of Gamma^3 is **2.40x** — the same kind of relabel — which is a
  strong hint the canoes-vs-fastnc channel/leg-conjugation assignment is part of
  the residual. See `scripts/docs/shear3pcf_natural_map.md`.

---

## Trusted references and key artifacts (all under this folder)

- `outputs/stage1_bkappa_phase_ref.npz` — the **decisive** B_kappa SPT phase
  reference (= fastnc to ~1%). This is ground truth for the radial shape.
- `outputs/stage1_fastnc.npz` — `fastnc` tree-level SAS run.
- `outputs/stage1_threeway_FINAL.npz` — BEFORE/AFTER frame-invariants + J_2
  residual (the table above lives here).
- `outputs/Wcan_vs_Bkappa_decomp.npz` — the W_can/B_kappa +0.44 ell-slope.
- `outputs/scalar_TTT_vs_Bkappa.npz`, `outputs/review_3pt_vs_pyccl_FINAL.npz`.
- Docs: `scripts/docs/decisive_bkappa_reference.md`,
  `scripts/docs/shear3pcf_natural_map.md`,
  `scripts/docs/scalar_pyccl_reconciliation.md`.
- Key scripts: `stage1_bkappa_phase_reference.py` (builds the reference),
  `stage1_threeway_compare_FIXED.py` / `stage1_threeway_AFTERFIX.py`,
  `stage1_fastnc_tree.py` (fastnc driver).

`fastnc` conventions (Sugiyama+2024, arXiv:2407.01798): h-unit native
(chi -> Mpc/h, k=ell/chi in h/Mpc); SAS isoceles t1=t2=gamma + opening angle
phi, t3=2*gamma*sin(phi/2); centroid projection (`x2cent`).

---

## Suggested first-principles plan (boundary-validation + xAct)

1. **Decouple radial from angular.** Feed the PROVEN-correct J_2 kernel two
   radial inputs in turn: (a) the B_kappa phase-reference radial shape, (b) the
   canoes HIGH B_delta radial shape. Confirm the phi-trend `[0.68, 0.43, 2.03]`
   is reproduced purely by swapping the radial shape. If yes, the residual is
   100% radial ell-shape; the angular side is exonerated for good.

2. **Characterize the +0.44 ell-slope as a first-principles object.** Compare
   the canoes tree-level matter bispectrum `B_delta(k1,k2,k3)` ell-shape against
   the SPT `B_kappa` phase reference directly (no J_2). Is the +0.44 a missing/
   extra factor of `ell` (or `k`, or a `(1+z)` power that didn't cancel) in the
   canoes HIGH radial kernel? Derive the expected SPT tree B_kappa ell-scaling
   symbolically (xAct / `.wl`) and diff against canoes.

3. **Pin the channel-label / leg-conjugation convention.** Derive the
   natural-component channel map (`Gamma^0..3` <-> `zeta_PPP` / `zeta_D`) from
   the Sachs spin-2 helicity structure via xAct, NOT by analogy. Then check
   whether relabeling / re-conjugating one leg in canoes reproduces the `2.40x`
   `fnc/ref` Gamma^3 reshuffle. The leg-conjugation assignment is the most likely
   single-line convention bug.

4. **Promote to a boundary test.** Once root-caused, add a parametrized pytest
   sweeping (gamma, phi) including the extreme corners (phi -> 0, phi -> 180,
   gamma very small + very large, ell=2 + ell=5000) per the boundary-validation
   rule, pinning known-good values from the B_kappa reference.

---

## Constraints (carry these — from CLAUDE.md / global rules)

- **First principles is the arbiter, NOT "paper wins."** The paper formalism may
  itself be wrong; verify against first principles (the user was explicit).
- Symbolic derivations: ALWAYS via Mathematica/xAct `.wl` scripts; never hand-
  derive in chat.
- Do NOT edit `sections/*.tex` during analysis (consume as ground truth). Never
  edit `biblio.bib` (surface references to the user).
- Keep local `E` primitive; `E = E_0/a` only in Born-labelled paragraphs.
- Conda interpreters (subagent Bash cannot `conda activate` — use full paths):
  - PyCCL:    `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`
  - sft-wick: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`
- Use the boundary-validation methodology (test BOTH methods at the threshold;
  always include extreme parameter corners). Use parallel Explore/subagents for
  broad sweeps; cross-check their claims.
- Archive (mv), don't delete, old artifacts. Don't commit/push without asking.
- End every reply with a 2-4 sentence Chinese summary titled "中文小结".
