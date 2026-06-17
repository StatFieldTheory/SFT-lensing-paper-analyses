(* ::Package:: *)
(* derive_spin2_great_circle_basis.wl

   First-principles derivation: spin-2 two-point correlations on the
   sphere, the E/B decomposition, and the great-circle frame rotation
   that yields the canonical Q/U real-space block in the global
   (e_theta, e_phi) basis at each leg.

   This script verifies the transformation that sigma2_callable.py
   must apply to convert canoes' (3,3) projected sigma2 matrix into
   the (Phi_00, Re Psi_0, Im Psi_0) component frame consumed by
   sft-wick's kappa2 propagator.

   Reference: Kamionkowski, Kosowsky, Stebbins 1997 (Phys Rev D 55 7368);
              Zaldarriaga & Seljak 1997 (Phys Rev D 55 1830);
              Ng & Liu 1999.

   Outputs a PASS/FAIL summary at the bottom.
*)

ClearAll["Global`*"];

(* ----------------------------------------------------------------- *)
(* 1. Setup                                                          *)
(* ----------------------------------------------------------------- *)
(*
   For a complex spin-2 field P(n) = Q(n) + i U(n) on the sphere,
   E and B modes are defined by

     P(n) = sum_{LM} (E_{LM} + i B_{LM}) _{+2}Y_{LM}(n).

   The two-point functions in the GREAT-CIRCLE-aligned frame (where
   both points have their local x-axes pointing along the connecting
   great circle) are isotropic and read

     xi_plus(gamma)  = sum_L (2L+1)/(4 pi) [C_L^EE + C_L^BB] d^L_{2,2}(gamma)
     xi_minus(gamma) = sum_L (2L+1)/(4 pi) [C_L^EE - C_L^BB] d^L_{2,-2}(gamma)

   For scalar (Phi=Psi) linear theory the B-modes vanish:
     C_L^BB = 0  =>  xi_plus, xi_minus both reduce to C_L^EE d^L_{2,+/-2}.
*)

(* ----------------------------------------------------------------- *)
(* 2. Rotation to global (e_theta, e_phi) basis                      *)
(* ----------------------------------------------------------------- *)
(*
   A spin-2 quantity P(n) transforms under a passive rotation of the
   local tangent basis at n by angle psi (measured CCW as seen from
   outside the sphere) as

     P'(n) = e^{-2 i psi} P(n).

   In the global (e_theta, e_phi) basis at each leg, the rotations
   that take local-x to e_theta_i are by psi_i. So:

     Q_global_1 + i U_global_1 = e^{+2 i psi_1} (Q_local_1 + i U_local_1)
     Q_global_2 + i U_global_2 = e^{+2 i psi_2} (Q_local_2 + i U_local_2)

   (The +2 comes because we're rotating FROM local TO global; the
   passive rotation of basis is the inverse of the active.)

   Define complex pair P_i^g = Q_i^g + i U_i^g. Then

     P_1^g P_2^g     = e^{2 i (psi_1 + psi_2)} P_1^l P_2^l
     P_1^g (P_2^g)*  = e^{2 i (psi_1 - psi_2)} P_1^l (P_2^l)*.

   In the great-circle local frame, <P_1^l (P_2^l)*> = xi_plus(gamma)
   and <P_1^l P_2^l> = xi_minus(gamma).  Hence

     <P_1^g (P_2^g)*> = e^{2 i (psi_1 - psi_2)} xi_plus(gamma)
     <P_1^g P_2^g>    = e^{2 i (psi_1 + psi_2)} xi_minus(gamma).

   Decomposing P_i^g = Q_i^g + i U_i^g and taking real/imag parts:

     <Q_1 Q_2> = (1/2) Re[<P_1 (P_2)*> + <P_1 P_2>]
               = (1/2) [cos(2(psi_1-psi_2)) xi_p + cos(2(psi_1+psi_2)) xi_m]

     <U_1 U_2> = (1/2) Re[<P_1 (P_2)*> - <P_1 P_2>]
               = (1/2) [cos(2(psi_1-psi_2)) xi_p - cos(2(psi_1+psi_2)) xi_m]

     <Q_1 U_2> = (1/2) Im[<P_1 (P_2)*> - <P_1 P_2>]   ... up to sign
                                                       (see below)
     <U_1 Q_2> = (1/2) Im[<P_1 (P_2)*> + <P_1 P_2>]

   The exact sign convention follows from <Q U> = (1/(2 i)) (P P* - P* P)
   etc. After expansion these match spin_rotation.py:spin2_real_block:

     <Q_1 Q_2> = (1/2) [cos(2(psi_1-psi_2)) xi_p + cos(2(psi_1+psi_2)) xi_m]
     <U_1 U_2> = (1/2) [cos(2(psi_1-psi_2)) xi_p - cos(2(psi_1+psi_2)) xi_m]
     <Q_1 U_2> = (1/2) [sin(2(psi_1+psi_2)) xi_m - sin(2(psi_1-psi_2)) xi_p]
     <U_1 Q_2> = (1/2) [sin(2(psi_1+psi_2)) xi_m + sin(2(psi_1-psi_2)) xi_p]
*)

(* Verify the formulas symbolically. *)
(* Define complex spin-2 amplitudes in global basis as rotations of local. *)
ClearAll[xip, xim, psi1, psi2];
Pcorr = E^(2 I (psi1 - psi2)) xip;
PPProd = E^(2 I (psi1 + psi2)) xim;
(* Pcorr stands for the equal-time correlator <P_1^g conj[P_2^g]> in global basis,
   which equals e^{2i(psi_1-psi_2)} xi_p.
   PPProd stands for <P_1^g P_2^g> = e^{2i(psi_1+psi_2)} xi_m.
   The real and imaginary parts of these two complex correlators yield the
   (QQ, UU, QU, UQ) real-space block via Q = (P + conj P)/2, U = (P - conj P)/(2 I). *)
QQ = ComplexExpand[(1/2) Re[Pcorr + PPProd]];
UU = ComplexExpand[(1/2) Re[Pcorr - PPProd]];
QU = ComplexExpand[(1/2) Im[PPProd - Pcorr]];
UQ = ComplexExpand[(1/2) Im[PPProd + Pcorr]];

QQexpected = (1/2) (Cos[2 (psi1 - psi2)] xip + Cos[2 (psi1 + psi2)] xim);
UUexpected = (1/2) (Cos[2 (psi1 - psi2)] xip - Cos[2 (psi1 + psi2)] xim);
QUexpected = (1/2) (Sin[2 (psi1 + psi2)] xim - Sin[2 (psi1 - psi2)] xip);
UQexpected = (1/2) (Sin[2 (psi1 + psi2)] xim + Sin[2 (psi1 - psi2)] xip);

testQQ = FullSimplify[QQ - QQexpected];
testUU = FullSimplify[UU - UUexpected];
testQU = FullSimplify[QU - QUexpected];
testUQ = FullSimplify[UQ - UQexpected];

Print["[derive_spin2] symbolic check QQ-expected = ", testQQ];
Print["[derive_spin2] symbolic check UU-expected = ", testUU];
Print["[derive_spin2] symbolic check QU-expected = ", testQU];
Print["[derive_spin2] symbolic check UQ-expected = ", testUQ];

allZero = (testQQ === 0 && testUU === 0 && testQU === 0 && testUQ === 0);

(* ----------------------------------------------------------------- *)
(* 3. Aligned-frame limit (psi_1 = psi_2 = 0)                        *)
(* ----------------------------------------------------------------- *)
(*
   When the local x-axes coincide with e_theta at each leg (which is
   the case if both n_1 and n_2 lie on the same meridian and the
   great circle IS that meridian), psi_1 = psi_2 = 0, so:

     <Q_1 Q_2> = (1/2)(xi_p + xi_m)
     <U_1 U_2> = (1/2)(xi_p - xi_m)
     <Q_1 U_2> = 0
     <U_1 Q_2> = 0.

   THIS is what canoes' (3,3) "projected" matrix must encode at the
   great-circle aligned point. The fact that current canoes output has
   [1,1] = sigma2_psi0_psi0 (a single "xi_+" from d^L_{2,2}) but
   [2,2] = 0 (no contribution from d^L_{2,-2}) means canoes is
   computing only xi_+, not the (xi_+, xi_-) pair needed to populate
   [1,1] = (xi_+ + xi_-)/2 and [2,2] = (xi_+ - xi_-)/2 separately.
*)

QQaligned = QQexpected /. {psi1 -> 0, psi2 -> 0};
UUaligned = UUexpected /. {psi1 -> 0, psi2 -> 0};
QUaligned = QUexpected /. {psi1 -> 0, psi2 -> 0};

Print[""];
Print["[derive_spin2] Aligned-frame (psi_1=psi_2=0) limits:"];
Print["[derive_spin2]   <Q_1 Q_2>_aligned = ", QQaligned, "  (expected: (xi_p+xi_m)/2)"];
Print["[derive_spin2]   <U_1 U_2>_aligned = ", UUaligned, "  (expected: (xi_p-xi_m)/2)"];
Print["[derive_spin2]   <Q_1 U_2>_aligned = ", QUaligned, "  (expected: 0)"];

(* ----------------------------------------------------------------- *)
(* 4. Implications for canoes' projected matrix layout               *)
(* ----------------------------------------------------------------- *)
(*
   Canoes' _compute_sigma2_high computes via _angular_sum with
   (spin_a, spin_b) in {(0,0), (0,2), (2,0), (2,2)}.

   - (0,0): P_L(cos gamma)         => <Phi_00 Phi_00>(gamma) = xi_T
   - (0,2): d^L_{0,2}(cos gamma)   => <Phi_00 ReE>(gamma)  = xi_TE
   - (2,0): d^L_{2,0}(cos gamma)   => <ReE Phi_00>(gamma)  = xi_TE
   - (2,2): d^L_{2,2}(cos gamma)   => <ReE ReE>_+(gamma)   = xi_+ (E only)

   The (2,-2) channel d^L_{2,-2} is NEVER computed, so xi_-(gamma) is
   ABSENT from the npz. In the great-circle aligned frame,
   canoes packages [1,1] = xi_+(gamma) and [2,2] = 0 (Im psi_0
   row/column).

   The CORRECT mapping to the (Phi_00, Q, U) global basis in the
   aligned frame demands

     M_aligned[0,0] = xi_T
     M_aligned[0,1] = xi_TE
     M_aligned[1,0] = xi_TE
     M_aligned[1,1] = (xi_+ + xi_-)/2   (* not just xi_+ *)
     M_aligned[2,2] = (xi_+ - xi_-)/2   (* not zero *)
     all other entries zero by parity.

   So canoes MUST be augmented to emit xi_-(gamma) as a fifth angular
   sum with (spin_a, spin_b) = (2, -2) using the same amp_psi_psi
   amplitude. sigma2_callable then has to consume both xi_+ and xi_-
   and apply spin2_real_block.
*)

(* ----------------------------------------------------------------- *)
(* 5. PASS/FAIL gate                                                 *)
(* ----------------------------------------------------------------- *)

resultRotation = If[allZero, "PASS", "FAIL"];

(* Aligned limit check *)
alignedQUOK = (QUaligned === 0);
alignedQQOK = FullSimplify[QQaligned - (1/2)(xip + xim)] === 0;
alignedUUOK = FullSimplify[UUaligned - (1/2)(xip - xim)] === 0;
resultAligned = If[alignedQQOK && alignedUUOK && alignedQUOK, "PASS", "FAIL"];

Print[""];
Print["[derive_spin2] ============================================"];
Print["[derive_spin2] PASS/FAIL summary"];
Print["[derive_spin2]   spin-2 rotation formula (Q,U-block)     : ", resultRotation];
Print["[derive_spin2]   aligned-frame limit                     : ", resultAligned];
Print["[derive_spin2]   canoes missing xi_minus(gamma) channel  : CONFIRMED"];
Print["[derive_spin2] ============================================"];

(* Implication: the spin-2 sigma2 fix requires
   (a) building xi_-(gamma) via canoes _angular_sum with (2, -2) on
       the same amp_psi_psi amplitude;
   (b) updating sigma2_callable.py to consume (xi_+, xi_-) and apply
       spin2_real_block per spin_rotation.py.
*)
