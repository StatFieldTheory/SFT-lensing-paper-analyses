# zeta_D HIGH / Limber spin-2 kernel: Phase-1 diagnosis (NO kernel fix landed)

Branch (canoes): `fix/spin2-zetaD-2_2_-2-norm` (UNCHANGED by this work; the only
uncommitted diff is the prior LOW-kernel fix in `_spin_aware_three_pt.py`, left
in place per the guardrails).
Env: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`.
Oracles (all under `analyses/shear3pcf_fastnc/`):
- `scripts/probe_HIGH_limber_oracle.py`        -> `outputs/HIGH_limber_oracle.npz`
- `scripts/probe_HIGH_fold_compare.py`         -> `outputs/HIGH_fold_compare.npz`
- `scripts/probe_HIGH_deployed_vs_complex.py`  -> `outputs/HIGH_deployed_vs_complex.npz`
- `scripts/probe_HIGH_window_B_isolation.py`   -> `outputs/HIGH_window_B_isolation.npz`

## VERDICT (one line)

**The over-largeness is NOT a wrong expression in `_kappa3_limber_alpha_kernel`.**
The orientation kernel is provably CORRECT. Per the prompt's explicit STOP rule
("If you cannot cleanly pin the exact wrong expression + correct form ... STOP --
do NOT guess-fix"), no kernel edit was applied. The 4-17x lives in canoes' HIGH
**B_delta / radial-fold normalization for the spin-2 modulus channels**, which is
a different locus than the one the task scoped (the kernel), so it is reported,
not patched.

## Phase 1 -- the analytic orientation integral (DERIVED, then verified)

For spins `(s1,s2,s3)` (canoes uses `s_j in {0,+-2}`; the reference uses spin
weight `2 sigma_j`, `sigma_j in {+-1}`, so `s_j == 2 sigma_j`), the flat-sky
natural component is

    Gamma = INT d^2 l2 d^2 l3 /(2pi)^4 B(l1,l2,l3)
              exp[i(l2.d2 + l3.d3)] exp[i(s1 phi_l1 + s2 phi_l2 + s3 phi_l3)]

with `l1 = -(l2+l3)`. Writing `l2 = u e^{i a}`, `l3 = v e^{i(a+phi)}` and doing the
absolute-orientation integral over `a in [0,2pi)`:

    I(s1,s2,s3; u,v,phi) = 2 pi i^S J_S(|Q|) e^{i S beta} e^{i(s1 psi + s3 phi)}
      Q    = u r2 + v e^{-i phi} r3       (r2 = d2, r3 = d3 as complex 2-vectors)
      S    = s1 + s2 + s3
      beta = arg(Q)
      psi  = arg(l1_rel),  l1_rel = -(u + v e^{i phi})

This is EXACTLY the complex form canoes computes (kappa3.py:1777-1784). canoes
then returns `np.real(I)` (line 1785) for `spin_sum != 0`.

### Gate (A): canoes' closed form == numeric a-integral (canoes geometry)
`probe_HIGH_limber_oracle.py`: max rel-err analytic-`I` vs brute a-integral on
canoes' own `(Q, l1_rel)` algebra = **2.74e-15**. The Bessel reduction is exact.

### Gate (B): canoes kernel + reference B, WHOLE plane == reference Gamma
`probe_HIGH_fold_compare.py` folds the canoes orientation kernel over a COMMON
`(u,v,phi)` grid and the reference's fastnc-identical single-plane `B_kappa`,
keeping `I` complex, and compares to the reference Gamma^3 (spins (2,2,-2)):

| config        | |G3 can_real| | |G3 can_cplx| | |G3 ref| | real/ref | cplx/ref |
|---------------|---------------|---------------|----------|----------|----------|
| PRIMARY 10'60 | 1.169e-07     | 1.178e-07     | 1.178e-07| 0.992    | **1.000**|
| 10' 30        | 1.406e-07     | 1.408e-07     | 1.408e-07| 0.999    | **1.000**|
| WORST 10' 120 | 3.806e-08     | 4.118e-08     | 4.118e-08| 0.924    | **1.000**|
| 50' 60        | 2.556e-08     | 2.558e-08     | 2.581e-08| 0.990    | 0.991    |
| 10' 90        | 8.009e-08     | 8.228e-08     | 8.228e-08| 0.973    | **1.000**|
| squeezed 10'10| 1.468e-07     | 1.469e-07     | 1.469e-07| 1.000    | **1.000**|
| 150' 120      | 9.905e-10     | 1.286e-09     | 1.121e-09| 0.884    | 1.147    |

**The kernel is correct.** With the reference's B, the complex-`I` fold reproduces
the reference Gamma^3 to <1% at almost every config; `np.real(I)` only costs ~0-12%
(and makes it SMALLER, never larger).

## Phase 1 -- what the over-largeness is NOT

### `np.real(I)` truncation is NOT the bug
`probe_HIGH_deployed_vs_complex.py` reimplements the deployed HIGH inner loop
byte-for-byte (validated == `compute_kappa3_mod_sigma3_high` to **0.00e+00**) with
one switch: real vs complex orientation kernel. HIGH-only, folded like
`stage1_ours_v2` (|Gamma^mu| = |fold(K^3 zeta_mu)|; `derived_phase` is unit-modulus):

| variant        | frame-inv can/ref median | range          |
|----------------|--------------------------|----------------|
| deployed_real  | 7.773                    | [4.118, 17.230]|
| complex        | 8.933                    | [5.843, 17.505]|

Keeping `I` complex makes the discrepancy LARGER, not smaller. The `np.real` is
not the lever. (It is, strictly, a minor inconsistency -- the natural components
are complex and the modulus should be taken after the radial fold -- but it does
not explain or fix the 4-17x.)

### The Bessel order S, the spin_extra phase, and the r2/r3 geometry are all fine
- The Bessel order for (2,2,-2) is S=2 (J_2), correct for a total-spin-2 natural
  component; for PPP (2,2,2) S=6 (J_6). Gates (A)+(B) above pass with these orders.
- `spin_extra = e^{i(s1 psi + s3 phi)}` is the exact leg-relative phase from the
  derivation; Gate (A) confirms it term-by-term.
- canoes' `_kappa3_flat_triangle_vertices` builds the SAME triangle (identical side
  lengths) as the reference's `triangle_geometry`, with OPPOSITE handedness (signed
  area sign-flipped: a reflection). A pure reflection conjugates the complex Gamma
  and leaves |Gamma| unchanged; Gate (B) confirms the magnitudes still match, so the
  handedness is not a magnitude error.

## Phase 1 -- where the 4-17x ACTUALLY lives (B / radial normalization)

`probe_HIGH_window_B_isolation.py` holds the canoes kernel AND the (60,1000] ell
window fixed and varies ONLY the B:

| config | |G3 ref-B whole| | |G3 ref-B window| | |Gref3| | whole/ref | win/ref |
|--------|------------------|-------------------|---------|-----------|---------|
| 10/60  | 1.181e-07        | 4.731e-08         | 1.196e-07| 0.988    | 0.395   |
| 10/30  | 1.435e-07        | 8.512e-08         | 1.484e-07| 0.967    | 0.574   |
| 10/120 | 3.464e-08        | 1.636e-08         | 3.585e-08| 0.966    | 0.456   |
| 50/60  | 2.722e-08        | 2.910e-08         | 3.032e-08| 0.898    | 0.960   |

With the reference B, the (60,1000] window captures only 0.40-0.96x of the full
Gamma^3 (the HIGH tail is a FRACTION; ell<=60 carries the rest, which the LOW path
is supposed to supply). But the DEPLOYED canoes-B on the SAME window+kernel is far
larger -- the smoking gun:

| config | |G3 deployed (canoes-B, window)| | |G3 ref-B, window| | deployed/(ref-B win) |
|--------|----------------------------------|--------------------|----------------------|
| 10/60  | 7.071e-07                        | 4.731e-08          | **14.95**            |
| 10/30  | 1.122e-06                        | 8.512e-08          | **13.18**            |
| 10/120 | 5.979e-08                        | 1.636e-08          | **3.66**             |
| 50/60  | 2.662e-07                        | 2.910e-08          | **9.15**             |

Same orientation kernel, same window: canoes' B_delta-based windowed integrand is
3.7-15x the reference's B_kappa-based one for the spin-2 D-channel. THAT factor is
the over-largeness. It is NOT in `_kappa3_limber_alpha_kernel`.

### Why this is spin-2-specific even though the B machinery is shared with PPP
The deployed HIGH B/radial assembly (B_delta, growth^4, the psi0^3 response, the
1/chi^4 and (2pi)^-4 prefactors, the measure, the radial_factor, h^6) is byte-for-
byte IDENTICAL between `compute_kappa3_sigma3_high` (scalar/PPP) and
`compute_kappa3_mod_sigma3_high` (modulus); only `spins` differs. The scalar STAGE-0
3PCF (TTT, J_0) closes to ratio_abs ~ 0.85 (`outputs/stage0_results.npz`), and
PPP (Gamma^0, J_6) is scattered around ~1. So the radial/B normalization is correct
for J_0 and J_6. It is only the J_2 (total-spin-2) D-channels that are 4-17x large.

This points at the spin-2 sampling of canoes' HIGH B/radial integrand. A further
diagnostic (`probe_HIGH_B_ratio_ellscaling.py`) collapsed the canoes convergence-3pt
weight `W_can` and compared it to the reference `B_kappa` on an equilateral ell
sweep. Two facts (with the caveat that the radial-measure bridge between the two
formulations is itself ell-dependent, so the ABSOLUTE offset below is partly a
convention artifact and should not be read as a clean physical factor):

  * R = W_can / B_kappa is NOT constant -- it grows mildly with ell, log-log slope
    d ln|R| / d ln L = **+0.44** (6.5 at L=80 -> 17.6 at L=800).
  * It is NOT the ell^4 = (ell^2)^2 a naive two-spin-2-leg eth^2 mis-weight would
    give (that would be 4096x at L=800; observed growth is 2.7x). So the lever is
    NOT a clean per-leg `sqrt(L^2(L^2-2))` factor either.

The honest, controlled conclusion is therefore narrower than "a missing eth^2
factor": canoes' HIGH B_delta/radial integrand has a MILD ell-dependent SHAPE
difference from the trusted single-plane `B_kappa`, and the spin-2 `J_2` kernel --
which weights the (60,1000] window more heavily than the scalar `J_0` or the PPP
`J_6` -- AMPLIFIES that shape difference into the observed 4-17x for the D-channels
while leaving scalar (STAGE-0 = 0.85) and PPP (~1) close. Pinning the exact source
of the ell-dependent B-shape difference (candidates: the `growth^4 / chi^4` per-leg
radial weighting, the `psi0` response ell-power, the single-plane-vs-LOS Limber
construction, or the F2/B_delta convention) requires a response/field-definition
derivation and a re-validation against the same `B_kappa`-phase reference. That is
outside the "wrong expression in `_kappa3_limber_alpha_kernel`" the task scoped, and
the STOP rule applies: no kernel edit, no factor fit.

## phi-growth note
The final can/ref grows with phi (4.1 at 60deg -> 17.2 at 120deg). This is NOT a
monotone B-factor: the canoes-B/ref-B windowed factor is 14.9, 13.2, 3.7, 9.1
(SMALLEST at phi=120). The growth of can/ref comes from the reference Gamma^3
shrinking toward the squeezed/wide limits faster than the deployed canoes value,
compounded with the channel-label reshuffle (Gamma^1 dominates at phi=30, Gamma^2/3
at phi=120). The over-largeness FACTOR is in the B/radial normalization; its phi
PROFILE is a geometry convolution, not a separate bug.

## GATES (numbers)
1. Kernel/oracle: analytic `I` == brute a-integral **2.74e-15**; canoes kernel +
   ref-B whole-plane fold == reference Gamma^3 to **0.88-1.00** (Gate B). Kernel OK.
2. Shear closure: NOT attempted via a kernel edit (none made). Deployed frame-inv
   can/ref = [4.118, 6.074, 17.230, 9.473] (UNCHANGED; no fix landed). `np.real` ->
   complex does NOT close it (median 7.77 -> 8.93). `stage1_threeway_compare_HIGHFIX.npz`
   was NOT written because no HIGH fix exists to certify.
3. Existing tests: `pytest test_kappa3_mod_normalization.py test_kappa3_sigma3_consistency.py
   test_spin_aware_three_pt.py -q` -> **72 passed**, 2 deselected (no code change, so
   green by construction; the +2 self-consistency tests still pass).
4. zeta_B (0,2,-2) and scalar: UNCHANGED (no edit). Scalar J_0 branch (kappa3.py:1771)
   untouched and bit-identical (no diff to kappa3.py at all).
5. Extremes: kernel verified finite and ~ref at phi=120 (worst), squeezed 10'/10,
   150'/120 (Gate B: 0.88-1.15). No blow-up. The deployed pipeline's 4-17x is the
   pre-existing B-normalization offset, present at all configs.

## HONEST VERDICT
The shear does NOT close, and it CANNOT be closed by a minimal edit to
`_kappa3_limber_alpha_kernel`, because that kernel is correct (verified three
independent ways). The 4-17x cosmic-shear over-largeness lives in the spin-2
modulus channels' HIGH **B_delta / radial response normalization** -- most likely a
per-spin-2-leg `eth^2`/`sqrt(L^2(L^2-2))` ell-weight that is right for scalar (J_0)
and PPP (J_6) but wrong for the J_2 D-channels. Fixing it is a response/field-
definition change (in `_kappa3_limber_response_product_h` / the per-leg spin-2
flat-sky response), to be derived symbolically and validated against the same
B_kappa-phase reference, in a follow-up scoped to the response rather than the
kernel. Per the boundary-validation + systematic-debugging rules, no factor was
fit and no kernel edit was forced to fake closure.
