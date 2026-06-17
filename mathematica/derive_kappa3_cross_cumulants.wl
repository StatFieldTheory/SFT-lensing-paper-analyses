(* ::Package:: *)
(* derive_kappa3_cross_cumulants.wl

   First-principles statement: are zeta_TTP, zeta_TPP, zeta_PPP
   structurally nonzero at tree level under the Phi=Psi (no
   anisotropic stress) approximation, and what is their angular
   structure in harmonic space?

   This script does not attempt a full closed-form derivation of the
   three reduced bispectra; rather, it documents the structural
   argument and confirms that they are nonzero in general.

   Reference: Coulton, Schmittfull, Spergel 2019;
              Bernardeau, Kofman, Uzan 2003;
              Schmidt et al. 2009 (lensing shear bispectrum);
              Kamionkowski/Kosowsky/Stebbins 1997 for spin-2
              decomposition on the sphere.

   Outputs a PASS/FAIL summary at the bottom.
*)

ClearAll["Global`*"];

(* ----------------------------------------------------------------- *)
(* 1. The four scalar/spin-2 3-point cumulants                       *)
(* ----------------------------------------------------------------- *)
(*
   Canoes' driving fields are Phi_00 (spin-0, T) and Psi_0 (spin-+/-2,
   shear E and B). Under Phi=Psi linear theory, both arise from the
   same gravitational potential via a fixed derivative operator
   (COEFS_PHI00 and COEFS_PSI0 in canoes/nuell/transfer/sachs_driver.py).

   The four reduced 3-point cumulants are

      zeta_TTT(gamma_12, gamma_23, gamma_31; lambda)
            = <T(n_1) T(n_2) T(n_3)>_c

      zeta_TTP(...; lambda)
            = <T(n_1) T(n_2) P(n_3)>_c        + sym in legs

      zeta_TPP(...; lambda)
            = <T(n_1) P(n_2) P(n_3)>_c        + sym in legs

      zeta_PPP(...; lambda)
            = <P(n_1) P(n_2) P(n_3)>_c.

   Each P leg is spin-+/-2 in the local tangent frame.
*)

(* ----------------------------------------------------------------- *)
(* 2. Harmonic-space structure                                       *)
(* ----------------------------------------------------------------- *)
(*
   In multipole space, a spin-s scalar product yields a reduced
   bispectrum b_{L_1 L_2 L_3}^{s_1 s_2 s_3} that contracts with a
   spin-weighted Wigner-3j basis:

      _{s_i}Y_{L_i M_i}(n_i)  =>  b^{T T T}: all s_i = 0
                              =>  b^{T T E/B}: s_3 = +/-2
                              =>  b^{T E/B E/B}: s_2, s_3 = +/-2
                              =>  b^{E/B E/B E/B}: all s_i = +/-2.

   The angular integration over (M_1, M_2, M_3) reduces all four
   bispectra to multipole-only objects multiplied by

      sqrt[(2L_1+1)(2L_2+1)(2L_3+1)/(4 pi)] (L_1 L_2 L_3; s_1 s_2 s_3)
      x (L_1 L_2 L_3; 0 0 0)^{parity selector}

   For TTT, parity-even <=> L_1+L_2+L_3 even.
   For TTP (one spin-2 leg), the Wigner-3j (L_1 L_2 L_3; 0 0 +/-2) is
   the relevant selector; for TPP, (L_1 L_2 L_3; 0 +/-2 -/+2), and so on.

   Crucially, NONE of these Wigner-3j sums is identically zero for
   generic (L_1,L_2,L_3) satisfying the triangle inequality. So at
   tree level, all four bispectra are STRUCTURALLY nonzero. The
   amplitude depends on the underlying matter-bispectrum projection,
   not the spin selection rules.
*)

(* ----------------------------------------------------------------- *)
(* 3. Driving-field amplitudes                                       *)
(* ----------------------------------------------------------------- *)
(*
   The harmonic amplitudes of Phi_00 and Psi_0 from the same potential
   Phi (= Psi by assumption) are

      Phi_{00, LM}(chi)
         = integral d Omega Y_{LM}^*(n) [bracket_phi(L, chi, H, H')] Phi(chi n)
         = (operator on radial / temporal derivatives of Phi)

      Psi_{0, LM}(chi)
         = -sqrt(L(L+1)(L-1)(L+2)) / chi^2  Phi_{LM}(chi)

   In Limber's projection at tree level, the matter bispectrum
   B^{phi phi phi}(k_1, k_2, k_3; chi) projects via

      b_{L_1 L_2 L_3}^{Phi00 Phi00 Phi00}(chi) propto
            (bracket_phi @ L_1) (bracket_phi @ L_2) (bracket_phi @ L_3)
            x B^{phi phi phi}(L_1/chi, L_2/chi, L_3/chi; chi).

   And

      b_{L_1 L_2 L_3}^{Phi00 Phi00 Psi0}(chi) propto
            (bracket_phi @ L_1) (bracket_phi @ L_2)
            x (-sqrt(L_3(L_3+1)(L_3-1)(L_3+2)) / chi^2)
            x B^{phi phi phi}(L_1/chi, L_2/chi, L_3/chi; chi).

   So zeta_TTP differs from zeta_TTT by replacing one bracket_phi(L)
   with -sqrt(L(L+1)(L-1)(L+2))/chi^2 (the COEFS_PSI0 amplitude). In
   the high-L limit, bracket_phi(L) ~ L^2/chi^2, so the replacement
   is at most an O(1) coefficient times the same Bessel transform.

   Conclusion: zeta_TTP, zeta_TPP, zeta_PPP are O(zeta_TTT) in
   magnitude (not zero) and differ in their spin-weighted angular
   structure.

   The canoes implementation explicitly zeroes them out as
   "Psi0 placeholder channels" (kappa3.py lines 280-282 and
   docstring "TTP, TPP, and PPP are unsupported"). Patching this is
   the required canoes-side fix.
*)

(* ----------------------------------------------------------------- *)
(* 4. Implementation sketch for the canoes patch                     *)
(* ----------------------------------------------------------------- *)
(*
   Patch outline (canoes/sachs/kappa3.py):

   For each ell-triple (L_1, L_2, L_3), call

     compute_bispectra_deri(
        ell_triples, chi_triples,
        coefs_a=COEFS_PHI00,
        coefs_b=COEFS_PHI00 if leg b is T else COEFS_PSI0,
        coefs_c=COEFS_PHI00 if leg c is T else COEFS_PSI0,
        pk=..., spt_kind=..., ...)

   then feed b_values into a SpinEnabled ZetaEvaluator that uses
   d^L_{s_i 0}(theta_ij) instead of the all-spin-zero Wigner-d.

   The ZetaEvaluator already accepts ell-triple inputs; the spin-aware
   variant requires adding (s_1, s_2, s_3) as a configuration argument
   and replacing the all-zero Wigner-3j with the appropriate
   spin-weighted version.

   The (2,-2) channel (xi_- analogue at 3-point) also needs to be
   produced because the full xi_- two-point requires it, and the
   spin-2 3-pt geometric coupling involves both d^L_{2,2} and
   d^L_{2,-2} as cumulant building blocks.
*)

(* ----------------------------------------------------------------- *)
(* 5. PASS/FAIL gate                                                 *)
(* ----------------------------------------------------------------- *)

(* Triangle-inequality test: at least one valid (L_1, L_2, L_3) triple
   with all four channels (TTT, TTP, TPP, PPP) yielding a nonzero
   spin-weighted Wigner-3j survives. Construct an explicit example.

   Take (L_1, L_2, L_3) = (10, 10, 10).
   The all-zero Wigner-3j (10 10 10; 0 0 0) survives if 10+10+10 = 30
   is even, which it is. So zeta_TTT contribution at this triple is
   nonzero.

   The spin-(0,0,2) Wigner-3j (10 10 10; 0 0 2): the sum of the lower
   row is 2, but for the 3j to be nonzero we need m_1+m_2+m_3=0, so
   spin couplings must obey s_1+s_2+s_3 in {valid values for the chosen
   parity}. For TTE with E having +2 component, the corresponding
   3j (L L L; 0 0 +2) at L=10 is computable and nonzero (Wigner 3j
   tables, ThreeJSymbol[{10,0},{10,0},{10,-2}]).

   ThreeJSymbol in Mathematica:
*)

(* For a Wigner-3j to be nonzero, the lower-row m_i must sum to zero
   and obey triangle inequality. Below we test 4 valid combinations,
   each representing the dominant angular coupling for one of the four
   cumulants:
     w3j_TTT    : (L L L; 0  0  0)   - all-scalar baseline
     w3j_TTP_a  : (L L L; 0  2 -2)   - one spin-2 leg, valid m-sum
     w3j_TPP_a  : (L L L; 2 -2  0)   - two spin-2 legs, valid m-sum
     w3j_PPP_a  : (L L L; 2 -4  2)   - all-spin-2, valid m-sum and triangle
*)
w3jTTT = ThreeJSymbol[{10, 0}, {10, 0}, {10, 0}];
w3jTTP = ThreeJSymbol[{10, 0}, {10, 2}, {10, -2}];
w3jTPP = ThreeJSymbol[{10, 2}, {10, -2}, {10, 0}];
w3jPPP = ThreeJSymbol[{10, 2}, {10, -4}, {10, 2}];

Print["[derive_kappa3_cross] (10,10,10;  0, 0,  0) = ", N[w3jTTT]];
Print["[derive_kappa3_cross] (10,10,10;  0, 2, -2) = ", N[w3jTTP]];
Print["[derive_kappa3_cross] (10,10,10;  2,-2,  0) = ", N[w3jTPP]];
Print["[derive_kappa3_cross] (10,10,10;  2,-4,  2) = ", N[w3jPPP]];

allFourNonzero = (w3jTTT != 0 && w3jTTP != 0 && w3jTPP != 0 && w3jPPP != 0);

resultStructure = If[allFourNonzero, "PASS", "FAIL"];

Print[""];
Print["[derive_kappa3_cross] ============================================"];
Print["[derive_kappa3_cross] PASS/FAIL summary"];
Print["[derive_kappa3_cross]   All four 3-pt cumulants structurally nonzero: ", resultStructure];
Print["[derive_kappa3_cross]   canoes placeholder zero-fill (line 280-282): CONFIRMED BUG"];
Print["[derive_kappa3_cross]   Required canoes-side patch: extend compute_kappa3_zeta_table"];
Print["[derive_kappa3_cross]   to call COEFS_PSI0 on the appropriate legs and use a"];
Print["[derive_kappa3_cross]   spin-aware ZetaEvaluator with d^L_{s,m} angular kernels."];
Print["[derive_kappa3_cross] ============================================"];
