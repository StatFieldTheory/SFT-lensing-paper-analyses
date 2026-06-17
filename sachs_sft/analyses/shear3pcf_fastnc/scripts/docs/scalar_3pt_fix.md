# Scalar 3-point normalization: root-cause power-count + scope verdict

Builder probes (DELIVERABLES):
- `scripts/probe_3pt_powercount.py` -> `outputs/powercount_3pt.npz` (both surpluses,
  machine-precision residuals)
- `scripts/probe_hpower_fast.py`    -> `outputs/hpower_fast.npz` (h-power, HIGH-only, fast)
- `scripts/probe_1pz4_fast.py`      -> `outputs/powercount_1pz4_fast.npz` ((1+z)^4 analytic)
- oracle: `scripts/probe_scalar_vs_pyccl_tree.py` -> `outputs/scalar_vs_pyccl_tree.npz`

Interpreters: sft-wick env for canoes builds, PyCCL env for the oracle (ABSOLUTE paths
in each script header).

## PHASE-1 RESULT (clean power-count, machine precision)

The deployed canoes scalar conv-3PCF (zeta_TTT fold) is ~10x (raw/dimensionless) over
the PyCCL standard at small gamma (oracle `realcan/std` ~ 9.95 at gamma<=12').  It
decomposes EXACTLY into two factors:

    realcan/std  =  (1+z_eff)^4   x   h^2
                    [surplus A]       [surplus B]
                    radial, per shell  units h-footing

Signal regime gamma<=12': (1+z_eff)^4 ~ 23.6, h^2 ~ 0.45, product ~ 10.6 vs measured
9.95 (Jensen gap <(1+z)^4 W>/<W> vs (1+z_eff)^4 ~ 6%).  They PARTIALLY CANCEL (one >1,
one <1), so the net is ~10x not ~50x.  Verified: after artificially removing the h^2
(h^6->h^4 trial), `realcan/std` rose 9.95 -> 22.08 ~ (1+z_eff)^4 EXACTLY (the trial was
then REVERTED -- see "Scope verdict" below.  2026-06-09 UPDATE: the h^6->h^4
piece was later RE-APPLIED to canoes and deployed once confirmed to be a genuine
dimensional correction; see the "SCOPE VERDICT STATUS (corrected 2026-06-09)"
note in the Scope verdict section.  Removing h^2 alone exposing the full
(1+z)^4 is expected: the radial (1+z)^4 piece is handled analysis-side.).

### SURPLUS (A) -- (1+z)^4 per shell  [machine precision]

`probe_3pt_powercount.py` SURPLUS (A) table: the per-shell radial-response ratio
deployed-canoes / PyCCL-standard:

    S_can = a^8 * K_geo^3 * [A(a)(1+z)^4]^3        (deployed: INT dlambda K_lambda^3
                                                    zeta, K_lambda=a^2 K_geo -> chi: a^8;
                                                    zeta carries per-leg Sachs (1+z)^4)
    S_std = g^3,  g = (3/2) Om H0^2 (1+z) K_geo    (PyCCL standard single-plane window)
    S_can/S_std = a^8 [A(a)(1+z)^4]^3 / [3/2 Om H0^2 (1+z) K_geo]^3 = -(1+z)^4
    (max|ratio/(-(1+z)^4) - 1| ~ 1e-15)

### SURPLUS (B) -- h^2  [machine precision]

`probe_3pt_powercount.py` SURPLUS (B) table: the canoes equal-shell zeta physical/h
ratio = h^6 EXACTLY (empirical h-power = 6.0000).  A dimensionless convergence 3PCF is
h-invariant; the dimensionless-correct physical footing is h^4 (the PyCCL oracle's own
h-test: a dimensionless conv-3PCF needs h^4).  So the canoes "physical" zeta is h^2
(= 1/h^2 ~ 2.22) larger than a dimensionless conv-3PCF.

## CORRECTED CODE MAP (why neither surplus is a single-line bug)

ALL FOUR equal-shell kappa3 builders apply the SAME h^6 "physical" footing AND the SAME
per-leg (1+z)^4 Sachs response -- CONSISTENTLY:

| builder (equal-shell, deployed)            | h-power line | (1+z)^4 per-leg |
|---|---|---|
| `compute_kappa3_zeta_table`     (scalar LOW)  | kappa3.py:4872 (h^6) | :511 (1+z)^12 = ((1+z)^4)^3 |
| `compute_kappa3_sigma3_high`    (scalar HIGH) | kappa3.py:3642 (h^6) | :1708 (1+z)^4 per leg |
| `compute_kappa3_mod_zeta_equal_time` (mod LOW)| kappa3.py:3906-07 (h^6) | shares TPP/PPP response |
| `compute_kappa3_mod_sigma3_high` (mod HIGH)   | kappa3.py:4050-51 (h^6) | shares TPP/PPP response |

The h^4 literals at kappa3.py:3151 and :3344 are NOT in the equal-shell path -- they
belong to `compute_kappa3_sigma3` (the CROSS-SHELL direct-sigma3 object, which has
integrated two radial lengths internally, so h^4 is correct for THAT object).  An
INITIAL hypothesis that the HIGH h^6 was inconsistent with a LOW h^4 was FALSIFIED by a
direct empirical h-power test: `compute_kappa3_zeta_table(tree_phi, sequential)`
physical/h = h^6 EXACTLY (empirical 6.0).  So LOW and HIGH are CONSISTENT at h^6; the
h^6 is the uniform canoes "physical" units convention (docstrings: "physical zeta = h^6
* h-native value"), NOT an internal inconsistency.

Likewise the per-leg (1+z)^4 is the paper's Sachs response: cosmology.tex eq.
"Phi00 scalar" gives Phi00 = E0^2 (1+z)^2 / a^2 * [transverse-Laplacian], = (1+z)^4 *
operator with E0=1 (local photon energy kept PRIMITIVE per
feedback_photon_energy_bookkeeping).  The a^2-full K = a^2 chi(chi_s-chi)/chi_s fold
(cosmology.tex eq. "K comoving kernel") cubed, with the LAMBDA-space equal-shell
collapse delta(lambda'_1-lambda'_2)delta(lambda'_1-lambda'_3) (appendix eq. "appendix
equal shell approx"), is implemented EXACTLY by the deployed stage0 fold.

## WHY THE 2-PT IS IMMUNE (3-point-specific, but as a CONVENTION not a bug)

The PyCCL-validated 2-pt (`compute_kappa2_grid`, a SEPARATE code path) matches PyCCL
to ~2% because the paper PROVES (insights.tex line 44-46, eq. "xi kappa comoving" + eq.
"BL kkkernel") its Order-0 reduces to the standard window W = (3/2)Om H0^2 (1+z)
chi(chi_s-chi)/chi_s -- i.e. the 2-pt projection reduces the per-leg (1+z)^4 photon-
energy DOWN to the standard (1+z)^1 in the final kernel, and uses h^4 (kappa2.py:364,
the dimensionless-correct footing).  The 3-pt equal-shell vertex (kappa3.py) does NOT
perform the same reduction -- it keeps the per-leg (1+z)^4 and the h^6 footing -- which
is why the surplus is 3-point-specific.

## SCOPE VERDICT -- STOP, REPORT, DO NOT GUESS-FIX

Per the boundary/guardrail rule ("If you cannot cleanly power-count the surplus to
specific lines -> STOP and report; do NOT guess-fix; do NOT fit a factor to PyCCL;
paper wins"):

- The surplus POWER is cleanly pinned to (1+z)^4 x h^2 (machine precision).
- But NEITHER factor is a single WRONG line:
  * (B) h^6 is the UNIFORM canoes "physical" units convention (all 4 equal-shell
    builders + documented as such).  Changing it to h^4 (the dimensionless-correct
    footing) is a CONVENTION change touching all 4 builders, and it ALONE makes the
    scalar oracle WORSE (9.95 -> 22.08, exposing the full (1+z)^4) -- so it cannot be
    landed in isolation.
  * (A) (1+z)^4 is the paper's per-leg local-E Sachs response (cosmology.tex eq.
    "Phi00 scalar") + the paper's lambda-space equal-shell collapse (appendix eq.
    "appendix equal shell approx").  The canoes code is FAITHFUL to the paper.  Per
    CLAUDE.md "the paper wins", auto-dividing (1+z)^4 would make the CODE contradict
    the PAPER.

- The real question is a PAPER/CONVENTION decision: does the STF convergence 3PCF adopt
  the standard dimensionless Born normalization (drop the local-E (1+z)^4 the way the
  2-pt projection does, and use h^4), or does it keep the local-E photon-energy
  bookkeeping?  This is the same physics the 2-pt resolves by reducing to the Born W;
  the 3-pt equal-shell collapse would need the analogous reduction DERIVED (the
  lambda-space delta-collapse vs the standard chi-space Limber collapse) -- a paper-
  level derivation, NOT a one-line canoes bugfix.

SCOPE VERDICT STATUS (corrected 2026-06-09): the verdict above (STOP / do-not-
guess-fix) applied to the RADIAL `(1+z)^4` piece (surplus A), which remains
analysis-side-only and was NOT folded into canoes.  The `units="physical"` h-power
piece (surplus B), however, WAS subsequently resolved as a genuine dimensional
correction and LANDED in canoes: native conv-3PCF units are `(h/Mpc)^4 =
A(a)^3 (P·P)/χ^4`, so the physical footing is h^4 not h^6.  The **h^6 -> h^4**
edit was applied to canoes `kappa3.py` at the 5 equal-shell sites (scalar
TTT/TTP/TPP/PPP `h4` block, the Bmod/Dmod LOW pair, the Bmod/Dmod HIGH pair, and
the per-shell `post` factor) on branch `fix/spin2-zetaD-2_2_-2-norm`, and the
DEPLOYED equal_time_limber table + FK sweep were REBUILT with h^4.  FK κκ shifted
by x2.2204 to **+1.219e-4** (γ=0.5'); Order-0 CCL κκ unchanged (0.998).  This
SUPERSEDES the "GATES (recorded for the trial h^6->h^4 before it was reverted)"
section below: the h^4 trial was NOT permanently reverted — it was re-applied and
deployed after the units argument was confirmed dimensionally correct.  The
branch also retains the prior LOW-kernel spin-2 fix
(`_spin_aware_three_pt.py`, uncommitted).

## SPIN-2 (subset) -- expected behaviour under any (A)+(B) fix

The spin-2 zeta_D (2,2,-2) frame-invariant canoes/ref = [4.118, 6.074, 17.230, 9.473]
(at the 4 configs; from `outputs/stage1_threeway_compare.npz`).  The scalar (1+z)^4 x
h^2 is channel-INDEPENDENT (radial + units), so removing it shifts the spin-2 by the
SAME factor -- but the spin-2 also has a phi-DEPENDENT residual (4.1 at phi=60deg ->
17.2 at phi=120deg) that is NOT a constant radial/units factor (documented in
zetaD_HIGH_fix.md: a mild B-shape ell-slope amplified by the J_2 kernel).  That
phi-dependent residual is the GENUINE spin-2 question on clean footing and would
SURVIVE the scalar (A)+(B) normalization fix.

## GATES (recorded for the trial h^6->h^4 before it was reverted)

- 3-pt scalar oracle (trial B-only, REVERTED): realcan/std 9.95 -> 22.08 ~ (1+z_eff)^4
  (B alone does NOT close; it exposes the full (1+z)^4).
- 2-pt no-regression: structurally untouched (kappa2.py separate, h^4; never edited).
- Extremes (B-trial HIGH): gamma=150' squeezed + phi=120deg + 0.5'-squeezed all FINITE,
  no blow-up (max|zeta| 1.76e-11, zero NaN/Inf).
- pytest: not completed (slow); the trial was reverted before landing, so no pins changed.
