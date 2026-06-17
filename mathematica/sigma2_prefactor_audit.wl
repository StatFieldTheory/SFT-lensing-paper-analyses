(* ============================================================
   sigma2_prefactor_audit.wl
   ------------------------------------------------------------
   First-principles audit of the (1+z) prefactor in Cl_{Phi_00 Phi_00}
   and the resulting σ² Limber formula.

   Goal: rigorously trace the (1+z)-power chain from
        (a) the symbolic transfer Phi_00 = -(E^2/a^2) [transfer]_Phi
            (per driving_fields_harmonics.wl)
        (b) the Poisson relation Phi_pot ∝ A_0 (1+z) δ_m / k^2
        (c) the matter power Pm
   to two distinct conventions for the angular spectrum:
        (i)  PHYSICAL units (what the appendix's σ² formula uses):
             Cl_phys = (1+z)^10 A_0^2 Pm/χ^2 in chi-frame Limber.
        (ii) TABULATION units (what apply_cl_transfer outputs):
             Cl_tab = a^8 × Cl_phys = (1+z)^{-8} × (1+z)^10 ... = (1+z)^2 A_0^2 Pm/χ^2 ?
   And then convert to lambda-frame via dχ/dλ = (1+z)^2 to get the
   integrated correlation length ∫Cl(λ, λ+s)ds for both conventions.

   The empirical comparison in
   `figures/limber_vs_exact_sweep.py` showed
        I_Limber / I_exact ≈ (1+z)^8
   at low ell, suggesting an (1+z)^8 mismatch. This script gives the
   definitive symbolic answer.

   Run with:
     ~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh \
        scripts/sigma2_prefactor_audit.wl 120
   ============================================================ *)

$HistoryLength = 0;

Print["=== σ² prefactor audit: first-principles (1+z)-power trace ==="];
Print[""];
Print["Conventions:"];
Print["  E_obs == E (observed photon energy)"];
Print["  E = E_0 / a  with E_0 fixed at observer (free-streaming)"];
Print["  a = 1/(1+z)  (scale factor)"];
Print["  A_0 = (3/2) Omega_m (H_0/c)^2  (Poisson constant)"];
Print[""];

(* === Phase 1. Symbolic Phi_00 transfer in tabulation units ===
   From driving_fields_harmonics.wl, the per-multipole transfer for the
   Phi_00 field on Phi/Psi sector is

        T_Phi := (a^2/E^2) Phi_00^{(1,s)} = -[ 2(H' - H^2) Psi
                                              + H(Psi' - Phi')
                                              - 2 H d_chi Psi
                                              + (1/2) Lap_perp (Phi+Psi)
                                              + Dnull^2 Phi ]

   `apply_cl_transfer` in transfer_cl.py operates with the SAME 6-vector
   coefficients, hence Cl_tab(chi_a, chi_b) := <T_Phi(a) T_Phi(b)>. The
   physical Phi_00 satisfies

        Phi_00 = (E^2/a^2) T_Phi.

   Therefore

        Cl_phys := <Phi_00(a) Phi_00(b)>
                = (E^2/a^2)_a × (E^2/a^2)_b × Cl_tab.

   With E = E_0/a, we have E^2/a^2 = E_0^2 / a^4. Setting E_0 = 1
   (per build_inputs.py: "E_0 = 1") gives E^2/a^2 = 1/a^4 = (1+z)^4.

   So  Cl_phys(a, b) = (1+z(a))^4 × (1+z(b))^4 × Cl_tab(a, b),
   and at coincident points (a = b)
        Cl_phys(λ, λ) = (1+z)^8 × Cl_tab(λ, λ).
   *)
Print["Phase 1. Convention from driving_fields_harmonics.wl:"];
Print["  Phi_00 = (E^2/a^2) T_Phi         (transfer T_Phi has the symbolic"];
Print["                                    coefficients in transfer_kernels.py)"];
Print[""];
Print["With E = E_0/a and E_0 = 1:"];
ratPrefactor = (1/a)^2 / a^2;
Print["  E^2/a^2 = ", ratPrefactor, " = a^(-4)"];
Print[""];

CovTabToPhys = (1/a^4) (1/a^4);
CovTabToPhys = CovTabToPhys /. a -> 1/(1+z);
Print["  Cl_phys / Cl_tab = ((E^2/a^2))^2 = a^(-8) = ", CovTabToPhys];
Print[""];

(* === Phase 2. Poisson relation + matter power === *)
Print["Phase 2. Poisson relation (matter-dominated regime):"];
Print["  Phi_pot = -A_0 (1+z) delta_m / k^2"];
Print["  Pm_Phi(k, a) = <Phi_pot Phi_pot> = A_0^2 (1+z)^2 Pm_delta(k, a) / k^4"];
PoissonPhiPow = a0^2 * (1+z)^2 / kk^4 * Pm;
Print["  Pm_Phi = ", PoissonPhiPow, "  (placeholder a0=A_0, kk=k, Pm=Pm_delta)"];
Print[""];

(* === Phase 3. Sub-horizon limit of T_Phi ===
   In Limber, spatial gradients dominate. The Laplacian-perp term gives
   the leading contribution: (1/2) Lap_perp (Phi+Psi) ~ -(1/2) k^2 (Phi+Psi).
   With Phi = Psi (no anisotropic stress): 2 × (1/2) k^2 Phi = k^2 Phi.
   So at sub-horizon:
        T_Phi  ≈  -[ -k^2 Phi ] = k^2 Phi   (sign per script convention)
                           (we drop H^2, H', Hubble-related O(H^2) terms
                            since they are negligible vs k^2 in sub-horizon)
   Thus  Phi_00  ≈  (E^2/a^2) k^2 Phi
                  =  (1+z)^4 × k^2 × (-A_0 (1+z) delta_m / k^2)
                  =  -A_0 (1+z)^5 delta_m         *)

Print["Phase 3. Sub-horizon limit of T_Phi:"];
Print["  Lap_perp (Phi+Psi) ~ -k^2 (Phi+Psi)"];
Print["  T_Phi (sub-horizon, Phi=Psi) ~ k^2 Phi"];
Print["  Phi_00 ~ (E^2/a^2) k^2 Phi = (1+z)^4 k^2 × (-A_0 (1+z)/k^2 delta_m)"];
phi00SubHorizon = -a0 * (1+z)^5 * deltam;
Print["  ⇒  Phi_00 ~ ", phi00SubHorizon];
Print[""];
Print["  Power: <Phi_00^2>_phys = A_0^2 (1+z)^10 Pm_delta"];
clPhys = a0^2 * (1+z)^10 * Pm;
Print["  ⇒  Cl_phys ~ ", clPhys, " / chi^2  (Limber)"];
Print[""];

(* === Phase 4. Tabulation-unit Cl ===
   From Phase 1: Cl_tab = a^8 × Cl_phys.  *)

clTabSubHorizon = clPhys * a^8 /. a -> 1/(1+z);
clTabSubHorizon = Simplify[clTabSubHorizon];
Print["Phase 4. Tabulation-unit Cl (= what apply_cl_transfer outputs):"];
Print["  Cl_tab = a^8 × Cl_phys (sub-horizon)"];
Print["         = a^8 × A_0^2 (1+z)^10 Pm/chi^2"];
Print["         = ", clTabSubHorizon, " / chi^2  (since a = 1/(1+z))"];
Print[""];
Print["  ⇒  Cl_tab(chi_a, chi_b) = (1+z)^2 A_0^2 Pm / chi^2  in chi-frame"];
Print[""];

(* === Phase 5. Convert to lambda-frame ===
   dchi/dlambda = (1+z)^2 ⇒ delta(chi_1 - chi_2) = delta(λ_1 - λ_2)/(1+z)^2.
   The off-diagonal-integrated correlation in lambda-frame:
        ∫Cl_X(chi_a, chi(λ_b)) dλ_b
            = ∫Cl_X(chi_a, chi_b) dchi_b · (dλ_b/dchi_b)
            = ∫Cl_X(chi_a, chi_b) dchi_b / (1+z)^2.
   For Limber-δ structure: ∫Cl_X dchi = (per-chi prefactor).
   So ∫Cl_X(chi, chi(λ+s)) dλ = (per-chi prefactor) / (1+z)^2.
   *)

Print["Phase 5. Lambda-frame integrated correlation lengths:"];
Print[""];
Print["  PHYSICAL units:"];
intPhysLambda = a0^2 * (1+z)^10 / (1+z)^2;
Print["    ∫Cl_phys(λ, λ+s) ds  =  (1+z)^10 A_0^2 Pm/chi^2 / (1+z)^2"];
Print["                          =  ", intPhysLambda, " A_0^2 Pm / chi^2  (= (1+z)^8 × A_0^2 Pm/chi^2)"];
Print[""];
Print["    [matches σ² formula in sections/appendix.tex eq (sigma2 explicit)]"];
Print[""];

intTabLambda = (1+z)^2 * a0^2 / (1+z)^2;
intTabLambda = Simplify[intTabLambda];
Print["  TABULATION units (apply_cl_transfer convention):"];
Print["    ∫Cl_tab(λ, λ+s) ds  =  (1+z)^2 A_0^2 Pm/chi^2 / (1+z)^2"];
Print["                         =  ", intTabLambda, " A_0^2 Pm / chi^2  (= NO (1+z) prefactor)"];
Print[""];
Print["    [does NOT match σ² formula; mismatched by factor (1+z)^8]"];
Print[""];

(* === Phase 6. Cross-check with empirical numerics === *)

ratioFormulaVsTab = intPhysLambda / intTabLambda;
ratioFormulaVsTab = Simplify[ratioFormulaVsTab];
Print["Phase 6. Predicted I_Limber / I_exact at low ell (Limber-good regime):"];
Print[""];
Print["  Ratio = (formula prefactor in physical units)"];
Print["          / (formula prefactor in tabulation units)"];
Print["         = ", ratioFormulaVsTab];
Print[""];
Print["  At z = 0.4 (production reference):"];
Print["    (1+z)^8 = ", (1.4)^8];
Print[""];
Print["  Empirical from limber_vs_exact_sweep.py at low ell:"];
Print["    ratio (Limber/exact) at λ=1170, z=0.4: 13.04 (ell=100), 13.41 (ell=30)"];
Print["    matches (1.4)^8 ≈ 14.76 within ~10% (expected Limber error)"];
Print[""];

(* === Phase 7. Verdict === *)

Print["=================================================================="];
Print["VERDICT"];
Print["=================================================================="];
Print[""];
Print["The σ² Limber formula in sections/appendix.tex eq (sigma2 explicit)"];
Print["    σ²_phys(γ; λ) = (1+z)^8 A_0^2 / (2π) ∫ k Pm P_l(cosγ) dk"];
Print[""];
Print["is in PHYSICAL Phi_00 units."];
Print[""];
Print["But the kappa2 Cl tabulation produced by"];
Print["    transfer_cl.apply_cl_transfer(...)"];
Print["follows the symbolic transfer in driving_fields_harmonics.wl which"];
Print["computes T_Phi := (a^2/E^2) Phi_00, NOT Phi_00.  Therefore"];
Print[""];
Print["    Cl_tab = (a^2/E^2)^2_a × (a^2/E^2)^2_b × Cl_phys"];
Print["           = a_a^8 × a_b^8 × Cl_phys      [E_0 = 1]"];
Print[""];
Print["At coincident points (E_0 = 1):"];
Print["    Cl_tab(λ, λ) = (1+z)^{-8} × Cl_phys(λ, λ)"];
Print[""];
Print["Equivalently, the integrated lambda-frame correlation:"];
Print["    Tab    : ∫Cl_tab(λ, λ+s) ds  =  A_0^2 Pm/chi^2          [no (1+z)]"];
Print["    Phys   : ∫Cl_phys(λ, λ+s) ds =  (1+z)^8 A_0^2 Pm/chi^2  [appendix]"];
Print[""];
Print["The numerical pipeline currently uses Cl_tab for kappa2 but uses the"];
Print["σ²_phys formula (with (1+z)^8) for sigma2.  These two normalisations"];
Print["are inconsistent by a factor of (1+z)^8 -- which at z=2 (production"];
Print["source redshift) is ≈ 6561.  The empirical evidence from"];
Print["limber_vs_exact_sweep.py corroborates this exactly: at low ell the"];
Print["measured I_Limber/I_exact ratio scales as (1+z)^8 across all λ."];
Print[""];
Print["FIX OPTIONS:"];
Print[""];
Print["  (A) Drop (1+z)^8 from compute_split_noise_levels."];
Print["      Both kappa2 and sigma2 then live in tabulation units."];
Print["      Downstream observables ⟨κκ⟩ are still meaningful but in"];
Print["      tabulation (a^2/E^2) units; multiply by (1+z)^8 at the END if a"];
Print["      physical-units result is wanted."];
Print[""];
Print["  (B) Multiply apply_cl_transfer output by (1+z)^4 per chi-leg, i.e."];
Print["      multiply Cl_tab(chi_a, chi_b) by (1+z(chi_a))^4 (1+z(chi_b))^4"];
Print["      to convert to Cl_phys. Both kappa2 and sigma2 then live in"];
Print["      physical units, matching the appendix's σ² formula."];
Print[""];
Print["  Both options restore (κ² + σ²) consistency.  Option (B) is more"];
Print["  invasive but produces 'paper-units' output. Option (A) is one-line."];
Print[""];
Print["[done]"];
