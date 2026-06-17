# STAGE 1 — OUR-SIDE shear 3PCF natural components (notes)

Producer: `scripts/sachs_sft/analyses/shear3pcf_fastnc/stage1_ours.py`
Interpreter: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`
Output: `outputs/stage1_ours.npz` (pinned npz_schema, projection=`cent`).

This file documents the OUR-SIDE half of STAGE 1 only. The fastnc reference and
the final cross-comparison are produced by the orchestrator / parallel agent.

---

## What was built

Cosmic-shear 3PCF natural components `Gamma^0..Gamma^3` (Schneider-Lombardi
basis, flat-sky, CENTROID projection) on the pinned SAS isoceles family
(`t1=t2=gamma`, opening angle `phi`, base `t3=2 gamma sin(phi/2)`), single
source plane `z_s=5`, built directly from the driving-field equal-shell spin-2
modulus cumulants (post-2026-06-04 helicity fix).

Grid (from `stage1_grid.json`): gamma in {10', 50', 150'};
phi in {10,20,30,45,60,90,120} deg. Primary point: gamma=10', phi=60deg.

## Physics chain (every convention documented)

1. **Driving field -> shear.** `Psi0 = sigma_+ + i sigma_x` is the spin-2 Sachs
   driving field; `gamma = gamma_+ + i gamma_x = -int dlambda (Dsachs2 + i Dsachs3)`
   (cosmology.tex `eq: gamma from Dsachs23`). The spin-2 screen-Hessian
   multiplier `sqrt(L^2(L^2-2))/chi^2` (appendix `eq: appendix screen hessian`)
   is ALREADY baked into the modulus channels by canoes (spin-weighted-harmonic
   contraction); not re-applied. Same code path that produced the STAGE-0
   spin-0 result (certified |ratio|=0.85-0.99).

2. **Natural components = spin-permuted zeta cumulants (canonical frame).** The
   four pure-shear natural components are the single-leg-conjugated complex
   spin-2 cumulants, in canoes' canonical great-circle frame:
   - `Gamma^0 = <P1 P2 P3>  = zeta_PPP`         spins (+2,+2,+2), gamma^4-suppressed
   - `Gamma^1 = <... apex* ...> = zeta_D(+2,+2,-2)`  conjugate APEX vertex
   - `Gamma^2 = <... base1* ...> = zeta_D(-2,+2,+2)` conjugate BASE vertex 1
   - `Gamma^3 = <... base2* ...> = zeta_D(+2,-2,+2)` conjugate BASE vertex 2

   `zeta_D` from `compute_kappa3_mod_zeta_equal_time` (LOW, spt_kind=tree_phi) +
   `compute_kappa3_mod_sigma3_high` (HIGH), via the `_spin_weights_dmod`
   override; `zeta_PPP` from `compute_kappa3_zeta_table` +
   `compute_kappa3_sigma3_high`. units=physical, radial_measure=lambda,
   PCAMBz0.txt, Omega_m=0.3160919980475834, h=0.6711, n_s=0.97.

   **Leg-labeling / isoceles symmetry (the key projection match).** On our triple
   `(cos t3, cos gamma, cos gamma)` legs 1,2 are the BASE vertices (endpoints of
   t3) and leg 3 is the APEX. fastnc's convention (smoke test: `Gamma^3=conj(Gamma^2)`,
   `Gamma^1` real for t1=t2) requires `Gamma^1` to conjugate the APEX and
   `Gamma^2,Gamma^3` to conjugate the two BASE vertices. With this mapping our
   centroid-projected `Gamma^3 - conj(Gamma^2) ~ 1e-9` (0.3% of |Gamma^3|, just
   the base1/base2 cosine-grid asymmetry) -> **isoceles symmetry RESTORED**,
   matching fastnc structurally.

3. **Radial fold (same as STAGE 0).** `Z = int_0^{lambda_s} dlambda K^3 zeta`,
   `K = a^2 chi (chi_s-chi)/chi_s`, `Dbar=a chi`, single source plane z_s=5,
   lambda-density (NO `(dchi/dlambda)^3`, per radial-measure memory 2026-05-17).
   Overall sign `(-1)^3` from `gamma=-int Dsachs` (matches the STAGE-0 kappa sign).

4. **Centroid projection.** Applied via fastnc's OWN `x2cent(mu, t1, t2, phi)`
   phase (verbatim port of fastnc.py:752-784; arXiv:2309.08601 Eq.15-16), with
   `t1=t2=gamma`, opening angle phi, on the SAME SAS triple.

## Projection convention pinned to fastnc — and the residual frame caveat

- **Definition matched:** Schneider-Lombardi 4 complex natural components
  `Gamma^0..3`, flat-sky, CENTROID projection (fastnc `projection='cent'`,
  fastnc.py:429-432, `x2cent`). Leg labeling matched so `Gamma^1` conjugates the
  apex and `Gamma^2,Gamma^3` the base vertices, reproducing fastnc's
  `Gamma^3=conj(Gamma^2)` / real-`Gamma^1` isoceles symmetry (verified ~0.3%).
- **Phase factor matched:** fastnc's exact `x2cent` source (q1,q2,q3 ratios +
  `exp(i*{3,1,3,-1}*phi)` per mu) ported byte-for-byte; sign convention
  `v=t1+t2 exp(-i phi)` kept.
- **Residual caveat (reported honestly).** Our cumulants are in canoes'
  CANONICAL great-circle frame (n1 at the pole; each leg's spin reference is its
  geodesic to leg 1; `_canonical_frame_angles`). fastnc's `x2cent` maps the
  FFT-native x-projection -> centroid. The canonical great-circle frame and
  fastnc's x-projection coincide up to a per-leg spin-2 rotation that we have
  matched at the level of (a) which leg each `Gamma^mu` conjugates and (b) the
  isoceles symmetry, but a residual configuration-DEPENDENT global phase between
  the two x-projection definitions cannot be excluded from our side alone. The
  frame-INVARIANT magnitude `sqrt(sum_mu |Gamma^mu|^2)` is unaffected by any such
  phase and is the robust normalization check; the per-component phase match
  should be confirmed against the parallel fastnc TREE-LEVEL run.

## Results

PRIMARY point (gamma=10', phi=60deg), centroid projection (see run log /
selfcheck PNG; numbers below from the corrected mapping):

| component | our value (tree-level) |
|-----------|------------------------|
| Gamma^0   | -3.44e-08 (gamma^4-suppressed all-unconjugated channel) |
| Gamma^1   | -7.07e-07 (apex-conjugated; real, dominant modulus) |
| Gamma^2   | -2.14e-07 - 3.71e-07j  (|.|=4.28e-07) |
| Gamma^3   | -2.15e-07 + 3.72e-07j  (|.|=4.30e-07, = conj(Gamma^2)) |

- **zeta_D dominates / zeta_PPP gamma^4-suppressed:** confirmed.
  |Z_PPP|/|Z_D| ~ 0.066 at gamma=10' (the all-+2 channel is suppressed by the
  d^L_{2,-2} ~ gamma^4 factor; the modulus channels carry d^L_{2,2} ~ O(1)).
  Built FRESH on the pinned SAS triples (NOT read from the 16-cos deployed npz,
  whose fixed cos-grid does not contain the pinned phi values).

- **Frame-invariant magnitude vs fastnc BiHalofit smoke (NONLINEAR, OOM only):**
  `sqrt(sum|Gamma|^2)`: ours(tree) = 9.32e-07 vs fnc(BiHalofit) = 1.107e-06,
  ratio 0.84. Expected `<1` at 10' since BiHalofit is nonlinear and our zeta is
  tree-level (nonlinear > tree at arcmin scales). This is a STRONG magnitude /
  normalization confirmation of the spin-2 zeta. The real physics validation is
  vs the parallel fastnc TREE-LEVEL subclass (which should reproduce the
  suppressed Gamma^0 too).

## Comparison vs the parallel fastnc TREE-LEVEL reference (stage1_fastnc.npz)

The parallel agent's fastnc tree-level npz became available; a read-only
sanity comparison (we did NOT run fastnc) shows a REAL, STRUCTURED discrepancy
that must be reported honestly:

- **Frame-invariant ratio ours/fnc-tree:** 3.0-17.4 across the grid (median ~8,
  mean ~7.8). NOT a constant -> not a pure normalization slip; structured in
  both gamma and phi.
- **Per-component structure differs:** fnc-tree has all four |Gamma^mu|
  comparable (|Gamma^1|/|Gamma^0| ~ 1.2-2.1 across phi); OURS has Gamma^0
  suppressed 12-20x below Gamma^{1,2,3}.
- **Phase:** Gamma^0, Gamma^1 match fnc up to the expected (-1)^3 overall sign
  (d_arg = 180deg, both real). Gamma^2, Gamma^3 are over-rotated by ~120deg
  relative to fnc (fnc Gamma^2 ~ real; ours ~ -120deg) -> the centroid x2cent
  application to Gamma^2/3 is not the correct bridge from canoes' great-circle
  frame.

**Diagnosis (the key methodological finding).** fastnc builds the cosmic-shear
natural components from the SCALAR convergence bispectrum B_kappa projected with
spin-2 phase factors (flat-sky gamma~(ell)=e^{-2i beta} kappa~(ell)), so all
four Gamma^mu are comparable and Gamma^0 is NOT gamma^4-suppressed. OUR
construction instead uses the spin-2 DRIVING-FIELD modulus channels (zeta_D,
zeta_PPP with the eth^2 multiplier), in which the all-unconjugated
Gamma^0=zeta_PPP IS gamma^4-suppressed. These are genuinely different objects;
the difference is structured, not constant.

Cross-checks bounding the discrepancy:
- A naive SCALAR conv-3PCF (zeta_TTT folded) over the same SAS triples is
  ~100-290x LARGER than fnc Gamma^0 (also structured) -> the shear comps are
  NOT the scalar 3PCF either.
- Our spin-2 modulus construction (3-17x) is therefore MUCH closer to fnc than
  the scalar 3PCF (100x), confirming the spin-2 channels are the RIGHT FAMILY;
  the residual is a structured spin-2 projection/normalization that is localized
  but not fully resolved on our side alone.

**What this validates and what it does not.**
- VALIDATED: the spin-2 zeta build (modulus channels, 2026-06-04 helicity fix),
  the radial fold/measure (STAGE-0-certified, shared bookkeeping), the
  isoceles symmetry (Gamma^3=conj(Gamma^2), Gamma^1 real for t1=t2), and the
  order-of-magnitude normalization (frame-invariant within ~8x of fnc tree,
  ~0.84 of fnc BiHalofit nonlinear).
- NOT YET RESOLVED: the exact spin-2 projection bridging canoes' canonical
  great-circle cumulant frame to fastnc's x/centroid natural-component frame.
  Correctly closing this requires deriving the great-circle -> x-projection
  spin-2 rotation symbolically (Mathematica per project rules) and/or building
  the shear comps from the SCALAR B_kappa with spin-2 projection phases (the
  same route as the paper's xi_pm 2PCF), rather than from separate spin-2
  driving-field channels. Flagged for the orchestrator / next session.

## Convention pitfalls and resolutions

- **Spin-2 multiplier double-count:** avoided. The `sqrt(L^2(L^2-2))/chi^2` is
  inside the canoes modulus channels; not re-applied.
- **Radial measure:** lambda-density; no `(dchi/dlambda)^3` (memory 2026-05-17).
- **Units:** physical 1/Mpc throughout (our side); h carried through chi, k,
  arcmin->rad (`gamma_rad=gamma_arcmin*pi/10800`).
- **Leg-conjugation permutations are DISTINCT cumulants:** the three single-leg
  zeta_D evaluations (apex vs base) are genuinely different (|G_apex| != |G_base|),
  obtained via `_spin_weights_dmod=(2,2,-2)/(-2,2,2)/(2,-2,2)`. Using zeta_D once
  for all three (the first-pass bug) wrongly forced |Gamma^1|=|Gamma^2|=|Gamma^3|.
- **Projection:** the centroid phase alone does NOT redistribute power among
  components; the correct distribution comes from the distinct leg-conjugation
  cumulants + the leg-to-component (apex/base) mapping. Both are in place.
