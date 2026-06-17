# zeta_D spin-2 3-point B-weight: Phase-1 derivation + oracle (NO fix landed)

Branch (canoes): `fix/spin2-zetaD-2_2_-2-norm` (UNCHANGED by this work; the only
uncommitted diff is the prior LOW-kernel fix in `_spin_aware_three_pt.py`, left
in place per the guardrails). No edit to any scalar/J_0/COEFS/operator/response
path was made.

Interpreter: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`.
All artifacts under `analyses/shear3pcf_fastnc/`:
- xAct derivation: `scripts/derive_spin2_3pt_bispectrum_weight.wl`
- oracles: `scripts/probe_psi0B_vs_correct.py` -> `outputs/psi0B_oracle.npz`;
  `scripts/probe_Wcan_vs_Bkappa_decomp.py` -> `outputs/Wcan_vs_Bkappa_decomp.npz`;
  `scripts/probe_perleg_factor_pin.py`;
  `scripts/probe_low_high_split.py` -> `outputs/low_high_split.npz`;
  `scripts/probe_scalar_TTT_vs_Bkappa.py` -> `outputs/scalar_TTT_vs_Bkappa.npz`.

## VERDICT (one line)

**NOT cleanly pinned as a spin-2-specific bug -> STOPPED, no fix.** The xAct
derivation PROVES the per-leg spin-2 (eth^2) ell-weight is NOT the lever (it ->1
with an O(1/ell^2) correction), and a decisive SCALAR test shows the canoes
convergence normalisation already disagrees with the B_kappa-phase reference by a
LARGER factor (~15-49x) than the spin-2 "4-17x" -- i.e. the spin-2 discrepancy is
a SUBSET of a global scalar+spin convention/normalisation mismatch between the two
reference families, NOT a spin-2-specific over-normalisation. There is no clean
wrong expression to fix without first resolving which scalar normalisation is
correct.

## Phase 1a -- the correct per-leg spin-2 3-point ell-weighting (xAct, DERIVED)

`derive_spin2_3pt_bispectrum_weight.wl` (all PASS):
- [T1] `(ell-1)ell(ell+1)(ell+2) == L2(L2-2)`, `L2=ell(ell+1)`  -> 0 (PASS).
- [T2] eth^2 (s=0->2) eigenvalue `== sqrt(L2(L2-2))` (PASS).
- [T3] per-leg spin-2/scalar modulus weight `W2(ell)=sqrt(L2(L2-2))/L2 -> 1`
  (PASS); this is exactly canoes' `COEFS_PSI0` ID `-sqrt(L2(L2-2))/chi^2` over
  `COEFS_PHI00` ID lead `L2/chi^2`.
- [T4] `W2(ell)-1 = -(3-2ell+2ell^2)/(2 ell^4) = O(1/ell^2)`. Numerically
  `1-W2 = 2.7e-4 @ell=60`, `1e-4 @ell=100`, `1e-6 @ell=500`.
- [T5] equilateral 3-leg `W2(L)^3 -> 1`; `1-W2^3 = 8.2e-4 @L=60`, `3e-4 @L=100`.

**Correct per-leg spin-2 ell-weighting** = `sqrt((ell-1)ell(ell+1)(ell+2))/L2`,
which -> 1 at high ell (the `e^{2i phi}` spin phase belongs in the orientation
kernel, NOT the |b|). So the correct high-ell psi0-B/phi00-B magnitude ratio is 1
(slope 0), exactly as the prompt anticipated.

**Consequence:** dropping `W2` (the canoes HIGH path's ell-flat `Psi0=-Phi00`
response does drop it -- see Phase 1c) changes the per-leg weight by AT MOST 0.03%
at ell=60 and falls as 1/ell^2. It CANNOT produce a +0.44 log-log slope or a
4-17x at ell~100-800. **The per-leg spin-2 eth^2 weight is NOT the +0.44 lever.**

## Phase 1b -- where the deployed zeta_D signal lives (LOW vs HIGH)

`probe_low_high_split.py` (`low_high_split.npz`), deployed fold at PRIMARY
(gamma=10', phi=60deg), |INT dlambda K^3 zeta|:

| channel | \|fold LOW\| | \|fold HIGH\| | HIGH/LOW |
|---|---|---|---|
| TTT | 1.59e-09 | 9.23e-06 | 5.79e+03 |
| PPP | 1.41e-22 | 3.44e-08 | 2.45e+14 |
| D3  | 2.20e-16 | 7.07e-07 | 3.21e+09 |

**ALL channels (TTT, PPP, D3) are HIGH-dominated by 3-14 orders.** So the LOW
operator path (`compute_bispectra_deri_psi0` + `SpinAwareZetaEvaluator`) is
irrelevant to the deployed value for every channel; the deployed shear is set by
the HIGH Limber path `compute_kappa3_{,_mod}sigma3_high`. This also refutes the
hypothesis "TTT closes because its LOW path is correct": TTT is HIGH-dominated too.

## Phase 1c -- the HIGH path is spin-blind in its B/radial assembly

`compute_kappa3_mod_sigma3_high` (kappa3.py:3912) and `compute_kappa3_sigma3_high`
(kappa3.py:3433) share a BYTE-IDENTICAL radial assembly:
`measure * growth^4 * b_delta_today(l_i/chi) / chi^4 * prefactor * response * rf * h^6`,
where `b_delta_today = 2(F2 PP + cyc)` (scalar SPT) and
`response = _kappa3_limber_response_product_h(channel)` returns the ell-FLAT
per-leg `Phi00 = +A(a)(1+z)^4`, `Psi0 = -A(a)(1+z)^4` (kappa3.py:1696-1714).
The ONLY thing that differs between TTT (J_0), PPP (J_6), and D (J_2) is the
Bessel order `S=s1+s2+s3` in `_kappa3_limber_alpha_kernel` (kappa3.py:1747).

So the per-leg spin-2 weight that the LOW operator and the 2-point apply
(`-sqrt(L2(L2-2))/chi^2`) is INDEED absent from the HIGH path (it uses ell-flat
`Psi0=-Phi00`). By Phase 1a that omission is a <=0.08% high-ell error, not the
4-17x. (The deployed `stage1_ours_notes.md` claim that the screen-Hessian
`sqrt(L2(L2-2))/chi^2` is "already baked into the modulus channels" is true only
for the subdominant LOW path; the dominant HIGH path uses the ell-flat response.)

## Phase 1d -- the +0.44 slope is an L-INDEPENDENT radial convention factor

`probe_Wcan_vs_Bkappa_decomp.py` (`Wcan_vs_Bkappa_decomp.npz`): per-triangle
folded ratio `R(L)=W_can_fold/B_kappa` (equilateral, PPP=scalar response):

| L | R | R/R0 |
|---|---|---|
| 60 | 5.36 | 1.00 |
| 100 | 7.35 | 1.37 |
| 300 | 13.1 | 2.45 |
| 1000 | 18.4 | 3.43 |

log-log slope `d ln|R|/d ln L = +0.4395` (reproduces the prior +0.44). BUT the
PER-SHELL radial ratio `S_can/S_ref` is **L-INDEPENDENT** (max/min = 1.000000,
constant = 52.19 at chi-index 20 across all L). So the +0.44 slope is NOT a
per-leg ell factor and NOT a per-leg radial-normalisation error: it is the
`b_delta(l/chi)` being weighted across chi by two DIFFERENT LOS measures (canoes
`K_g^3 a^2 growth^4/chi^4 [A(a)(1+z)^4]^3 h^6` vs reference
`g^3 (1/chi)(1+z)^3 D^4`, `g=1.5 H0^2 Om(1-chi/chi_s)`). `probe_perleg_factor_pin.py`
confirms the net per-shell `S_can/S_ref = -15.77` factorises exactly into a huge
per-leg-efficiency ratio (`K_g A(a)(1+z)^4 / g`)^3 times the compensating
`a^2/chi^4` vs `1/chi` carriers -- a pure convention re-distribution, the SAME
for all spins.

## Phase 1e -- DECISIVE: the SCALAR already disagrees with B_kappa by MORE

`probe_scalar_TTT_vs_Bkappa.py` (`scalar_TTT_vs_Bkappa.npz`) folds the SAME
deployed canoes `zeta_TTT` and compares to the B_kappa-phase reference's OWN
scalar (spins=0) convergence 3PCF:

| config | \|canoes_TTT\| | \|Bk_scalar\| | can/Bk | signs |
|---|---|---|---|---|
| (10', 60deg) | 9.23e-06 | 5.92e-07 | **15.59** | +/+ match |
| (50', 60deg) | 1.04e-07 | 2.11e-09 | **49.12** | -/- match |

Side-by-side with the spin-2 verdict (`stage1_threeway_compare.npz`,
frame-inv can/ref):

| config | SCALAR can/Bk | SPIN-2 can/ref |
|---|---|---|
| (10', 60deg) | **15.59** | 4.12 |
| (50', 60deg) | **49.12** | 9.47 |

**The scalar canoes/B_kappa mismatch is LARGER than the spin-2 one at every
config.** If `zeta_D` had a spin-2-specific over-normalisation, the scalar
`zeta_TTT` would AGREE with B_kappa (ratio ~1) while only spin-2 disagreed.
Instead the scalar disagrees MORE. So the spin-2 "4-17x" is a SUBSET of a global
~15-49x scalar+spin convention/normalisation difference between the
canoes/STAGE-0 family and the B_kappa-phase/fastnc family; the J_2 kernel
partially COMPENSATES it relative to J_0.

### Why STAGE-0 "validated" the scalar to 0.85 but it is 15.6x vs B_kappa
STAGE-0 (`stage0_spt_reference.py`) compares canoes `zeta_TTT` to a hand-rolled
SPT reference built with the SAME canoes radial convention
(`(1/(2pi)^4)(D^4/chi^4)[A(a)(1+z)^4]^3`, folded with `K=a^2 chi(chi_s-chi)/chi_s`
and `a^8 K_chi^3`). That is an INTERNAL same-convention consistency check (the
STAGE-0 docstring says so: "~30% agreement is the expectation"), NOT a test
against the B_kappa normalisation. The B_kappa-phase reference uses the standard
convergence Limber kernel `g^3(1/chi)(1+z)^3 D^4`. These two scalar conventions
differ by ~15-49x (Phase 1d/1e). The scalar normalisation was therefore NEVER
validated against B_kappa; the spin-2 verdict silently inherits that gap.

## Phase 4 -- NOT REACHED (STOP rule invoked)

Per the prompt's explicit STOP rule ("If you cannot cleanly pin the exact wrong
expression + the correct form (and show it yields the +0.44 slope -> 4-17x AND is
3-point-specific so the 2-point stays right), STOP and report -- do NOT
guess-fix"): the discrepancy is NOT 3-point-specific (it is present and larger in
the scalar 2-point-like J_0 fold), NOT spin-2-specific (scalar disagrees more),
and NOT a per-leg eth^2 weight (xAct: ->1). No fix was applied. No
`stage1_threeway_compare_BFIX.npz` was written (no fix to certify). The 2-point
shear path is untouched (no code change at all).

## What a CORRECT next step looks like (for the human)

The open question is a NORMALISATION/CONVENTION reconciliation in the SCALAR
sector, prior to any spin-2 work:
1. Decide which convergence normalisation is physically correct for the paper's
   driving-field fold: the canoes/STAGE-0 `(D^4/chi^4)[A(a)(1+z)^4]^3 * K^3`
   family, or the standard `g^3(1/chi)(1+z)^3 D^4` convergence-Limber family.
   They differ by a structured ~15-49x with a +0.44 ell-slope -- a real LOS
   measure / lensing-efficiency bookkeeping difference, NOT a spin-2 bug.
2. Cross-check the SCALAR canoes `zeta_TTT` fold against the B_kappa-phase scalar
   3PCF (and/or pyccl's tree convergence 3PCF) and bring it to ~1 FIRST.
3. ONLY THEN re-run the spin-2 three-way; if a residual spin-2-specific factor
   survives after the scalar is reconciled, it can be pinned cleanly (and the
   xAct-derived `sqrt(L2(L2-2))/L2` per-leg weight added to the HIGH response if
   the residual matches its O(1/ell^2) signature, which the current data do not).

The boundary-validation + systematic-debugging rules forbid fitting a factor to
fastnc or forcing closure; the honest finding is that the "spin-2 bug" is not
cleanly separable from an unresolved scalar normalisation-convention difference.
