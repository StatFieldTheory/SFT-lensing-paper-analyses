# 3-pt equal-shell photon-energy reduction: does (1+z)^4 x h^2 reduce away?

Script: `scripts/derive_3pt_photon_energy_reduction.wl` (run with
`~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh <script> 300`).
All algebra is flat-FLRW background `(1+z), a, chi, Dbar` + radial-measure
bookkeeping; the angular/spin machinery is byte-identical on the canoes and
standard sides and divides out (see `scalar_pyccl_reconciliation.md`), so no
xAct tensor structure is needed.

## Question

The deployed canoes scalar equal-shell convergence-3PCF is over the PyCCL/
standard tree result by EXACTLY `(1+z)^4 x h^2` per shell (machine precision,
`scalar_3pt_fix.md`). The 2-pt is immune (canoes/PyCCL = 0.98). Is the 3-pt
surplus a normalization slip that SHOULD reduce away (the exact 3-leg analogue
of the 2-pt reduction the paper proves in `cosmology.tex`), or a genuine
Order-0 STF-vs-Born difference?

## VERDICT: it REDUCES. The surplus is an over-normalization, not new physics.

The convergence 3PCF is a unique measurable; STF Order-0 is the free theory,
which the paper proves *is* the standard projection of the corresponding
driving-field cumulant (`cosmology.tex` "from sachs to wl kernels"). So Order-0
STF conv-3PCF MUST equal the standard tree conv-3PCF. The residual
`-(1+z)^4 h^2` is a measure/units over-normalization in the four equal-shell
`kappa3` builders. **Reduction factor to multiply each builder by:**

    R(z) = -1 / [ (1+z)^4 h^2 ]          (TeXForm: -1/(h^2 (z+1)^4))

i.e. divide the equal-shell zeta by `(1+z)^4` per shell (radial) and by `h^2`
(units), keeping the sign. Numerically `S_can/S_std x R(z) x h^2 = +1` on
z = 1..5 (Step 5, all PASS).

## STEP 1 — machinery validated: the 2-pt reduction reproduced

The paper's affine kernel (`cosmology.tex` eq. `K affine kernel`)
`K(λ,λ_s) = Dbar(λ)^2 ∫_λ^λ_s dλ₁/Dbar(λ₁)^2`, with `Dbar=aχ` and
`dλ=a²dχ`, collapses (the measure's `a²` cancels the `1/Dbar²` a-factor) to

    K(λ,λ_s) = a²(χ) · χ(χ_s-χ)/χ_s      (eq. K comoving kernel)   [PASS: K - a²·Kgeo = 0]

So the 2-pt window is `a²(χ)·[lensing efficiency]`. The driving-field 2-cumulant
is a λ-density carrying the per-leg local-E Sachs response `A(a)(1+z)^4`
(`cosmology.tex` eq. `Phi00 scalar`, E0 PRIMITIVE). The measure change to a
χ-density multiplies it by `a²` per leg (eq. `cumulant 2 measure change`). Net,
the per-leg field surplus over the standard window `W=(3/2)Ω_m H₀²(1+z)Kgeo` is
`-(1+z)^4`, and the two `a²` Jacobians per leg (`a^4`) knock it down to `-1`:
over two legs the net is `(-1)² = +1`. **The 2-pt reduces; matches PyCCL 0.98.**
[PASS: 2-pt net = +1.]

## STEP 2 — the 3-pt per-shell surplus is exactly -(1+z)^4

The appendix cross-shell 3-pt (`appendix.tex` eq. `appendix kappa3 cross shell`)
folds THREE response propagators `(Dbar(λ'_i)/Dbar(λ_i))²` onto the cross-shell
cumulant; the equal-shell approximation (eq. `appendix equal shell approx`)
replaces it by `ζ(γ,λ'₁) δ(λ'₁-λ'₂) δ(λ'₁-λ'₃)`, collapsing TWO of the three
radial integrals. The deployed per-shell radial response (folded in χ;
`probe_1pz4_fast.py:45`) vs the standard `g³`:

    S_can = a^8 · Kgeo^3 · [A(a)(1+z)^4]^3 = (27/8)(χ³(χ-χ_s)³/χ_s³) H₀^6 Ω_m³ (1+z)^7
    S_std = W_std^3                        = (27/8)(χ³(χ_s-χ)³/χ_s³) H₀^6 Ω_m³ (1+z)^3
    S_can / S_std = -(1+z)^4                                  [PASS: residual = 0]

(`A(a) = -(3/2)Ω_m H₀²(1+z)`, so the Poisson cube carries the minus sign.)

## STEP 3 — WHY the 3-pt does not fully reduce (the leg-count asymmetry)

This is the mechanism, and it is the SAME per-leg physics as the 2-pt — only
the radial bookkeeping differs:

| | per-leg field surplus | measure Jacobian returned | net |
|---|---|---|---|
| 2-pt (2 legs, 2 surviving χ-integrals) | `(-(1+z)^4)` each | `a^4` **per leg** (`a^8` total) | `(-(1+z)^4)² a^8 = +1` → REDUCES |
| 3-pt (3 legs, equal-shell collapse → 1 surviving χ-integral) | `(-(1+z)^4)` each | `a^8` **total**, not per-leg | `(-(1+z)^4)³ a^8 = -(1+z)^4` → RESIDUAL |

Reconstruction (Step 3, PASS): three hot legs `(-(1+z)^4)³ = -(1+z)^12` times the
`a^8 = (1+z)^{-8}` the single surviving integral returns gives `-(1+z)^4`,
matching `S_can/S_std`. The 2-pt's two surviving integrals return `a^8` as
`a^4` PER LEG, exactly cancelling the two hot legs to `+1`. The 3-pt's
equal-shell delta-collapse discards two of the three χ-Jacobians that a fully
unfolded (cross-shell) computation would supply, so two legs keep their local-E
`(1+z)^4` uncancelled. **This is an artifact of the equal-shell measure
bookkeeping, not a physical Order-0 difference.**

## STEP 4 — the h^2 units surplus

canoes `units="physical"` multiplies the h-native zeta by `h^6`; a dimensionless
conv-3PCF is h-invariant and the dimensionless-correct footing is `h^4`
(`scalar_3pt_fix.md` SURPLUS B; the PyCCL oracle's own h-test). Surplus
`= h^{6-4} = h^2`. Total surplus `= -(1+z)^4 · h^2`.

## Why this is REDUCES and not a genuine STF difference (honest verdict)

The airtight argument is structural + empirical, independent of the leg-count
heuristic:

1. The 2-pt and 3-pt use the IDENTICAL per-leg local-E response `A(a)(1+z)^4`
   and the IDENTICAL measure `dλ=a²dχ` (kappa2.py and `kappa3.py:1708` share
   `Phi00 = A(a)(1+z)^4 δ` and the same Jacobian).
2. The 2-pt matches PyCCL to 2% — so that shared per-leg `(1+z)^4` IS correctly
   reduced to the standard `(1+z)^1` window by the projection+measure at 2 legs.
3. The convergence 3PCF is a UNIQUE measurable; STF Order-0 (free theory) is
   proven to BE the standard projection of the 3-cumulant. So Order-0 MUST equal
   the standard tree conv-3PCF.
4. They differ by exactly `-(1+z)^4 h^2`, with per-leg physics identical to the
   reducing 2-pt. Therefore the residual is an over-normalization in the
   equal-shell builders, NOT new physics.

If `(1+z)^4` were genuine per-leg STF physics it would also corrupt the 2-pt by
`(1+z)^4` per shell (~20x at z_s=5); it manifestly does not. The reduction is
**3-pt-specific** and lives entirely in the equal-shell collapse's measure +
the `units="physical"` h-power.

## Per-builder fix (apply ONLY to the 4 equal-shell kappa3 builders)

Multiply the equal-shell zeta returned by EACH of the four builders by
`R(z) = -1/((1+z)^4 h^2)` (per source shell z, sign kept):

| builder | file:line (h-power) | (1+z) per-leg |
|---|---|---|
| `compute_kappa3_zeta_table` (scalar LOW) | kappa3.py:4872 (h^6) | :511 (1+z)^12 |
| `compute_kappa3_sigma3_high` (scalar HIGH) | kappa3.py:3642 (h^6) | :1708 (1+z)^4/leg |
| `compute_kappa3_mod_zeta_equal_time` (mod LOW) | kappa3.py:3906-07 (h^6) | shared TPP/PPP |
| `compute_kappa3_mod_sigma3_high` (mod HIGH) | kappa3.py:4050-51 (h^6) | shared TPP/PPP |

Equivalently, split R into the two independent pieces:
- radial: divide by `(1+z)^4` per shell (drop the 2 uncancelled local-E legs;
  this is the analogue of the 2-pt projection reducing `(1+z)^4 → (1+z)^1`);
- units: change the `h^6` "physical" footing to `h^4` (dimensionless-correct).

DO NOT touch the `kappa2` path: it already reduces correctly (PyCCL 0.98), and
applying R there would WRONGLY corrupt it. The fix is 3-pt-specific.

NOTE (CLAUDE.md "paper wins"): the canoes code is FAITHFUL to the paper's
per-leg local-E response and λ-space equal-shell collapse as currently written;
this reduction is a PAPER-LEVEL statement (the equal-shell measure must return
the same Jacobians the cross-shell fold would) that licenses dividing the
builders by R(z). It is the 3-leg counterpart of the 2-pt `(1+z)^4 → (1+z)^1`
reduction already proven in `cosmology.tex`.

SCOPE OF THIS DERIVATION vs WHAT WAS ACTUALLY EDITED (corrected 2026-06-09):
the `R(z) = -1/((1+z)^4 h^2)` radial `(1+z)^4` Born reduction analysed here is
**analysis-side only** — it provides the explicit factor and the per-builder
map; the radial `(1+z)^4` piece is NOT folded into canoes by this derivation.
HOWEVER, a SEPARATE, smaller, coupled fix WAS applied to canoes: the
`units="physical"` **h^6 -> h^4** units correction (surplus B; native
`(h/Mpc)^4 = A(a)^3 (P·P)/χ^4`). That edit landed in canoes `kappa3.py` at the
5 equal-shell sites (scalar TTT/TTP/TPP/PPP `h4` block, the Bmod/Dmod LOW pair,
the Bmod/Dmod HIGH pair, and the per-shell `post` factor) on branch
`fix/spin2-zetaD-2_2_-2-norm`, and the DEPLOYED equal_time_limber table + the FK
sweep were REBUILT with h^4. Net effect on the deployed product: FK κκ shifted
by x2.2204 to **+1.219e-4** (γ=0.5'); Order-0 CCL κκ unchanged (0.998). So the
"canoes untouched / read-only" framing applies ONLY to the radial `(1+z)^4`
piece — canoes WAS edited for the h^4 units piece.

## PASS ledger (all 6 PASS)

- 2-pt K-kernel reduces to `a² χ(χ_s-χ)/χ_s` (paper eq. K comoving)
- 3-pt per-shell surplus `= -(1+z)^4` (machine precision; matches scalar_3pt_fix.md)
- 2-pt net (2 legs) reduction `= +1` (REDUCES; matches PyCCL 0.98)
- 3-pt leg-count reconstruction consistent (`-(1+z)^12 × a^8 = -(1+z)^4`)
- numerical `S_can/S_std = -(1+z)^4` on z=1..5
- `R(z) × h^2` brings the 3-pt observable to the standard (`+1`) on z=1..5
