(* ::Package:: *)
(* derive_kappa3_spin2_helicity.wl

   First-principles derivation of the CORRECT spin-2 helicity / conjugation
   structure for the four driving-field three-point cumulant channels
   {TTT, TTP, TPP, PPP}, and the exact linear map from the computable
   helicity cumulants to the local real-component tensor

       K_{abc} = <phi_a phi_b phi_c>,   a,b,c in {0,1,2}
                 (0 = Phi00, 1 = sigma_+ = Re Psi0, 2 = sigma_x = Im Psi0)

   that sft-wick contracts (sft-wick works in this LOCAL real basis).

   WHY: canoes computes each spin-2 leg at helicity +2 only
   (spin_weights (0,2,2) for TPP, (2,2,2) for PPP) and takes Re(...).
   This yields the "++"/difference combination, which:
     (a) is the WRONG combination for the Sachs F-vertex contraction
         (the F-vertex couples sigma sigma-bar = |Psi0|^2 = sigma_+^2 + sigma_x^2,
          i.e. the SUM, which is the helicity (+2,-2) MODULUS channel), and
     (b) is gamma^4-suppressed at small angle (d^L_{2,-2} ~ gamma^4),
         whereas the needed modulus channel ~ d^L_{2,2} ~ O(1).

   This is the 3-point extension of derive_spin2_great_circle_basis.wl
   (2-point: <P P*> = xi_+ ~ d^L_{2,2}, <P P> = xi_- ~ d^L_{2,-2}).

   Convention: P = sigma_+ + i sigma_x (helicity +2), Pb = conjugate
   (helicity -2). sigma_+ = (P + Pb)/2, sigma_x = (P - Pb)/(2 i).
   Phi = Phi00 (spin 0, real). Statistical parity invariance about the
   separation axis swaps P <-> Pb, so every helicity cumulant equals its
   P<->Pb image (hence is real); cumulants odd in sigma_x vanish.

   Reference: Kamionkowski, Kosowsky, Stebbins 1997 (PRD 55 7368);
              Schmidt, Rozo, Dodelson, Hui, Sheldon 2009 (lensing shear bispectrum);
              and the 2-point sibling derive_spin2_great_circle_basis.wl.

   Outputs a PASS/FAIL summary at the bottom.
*)

ClearAll["Global`*"];

pass = {}; fail = {};
check[name_, cond_] := If[TrueQ[cond], AppendTo[pass, name], AppendTo[fail, name]];

(* ----------------------------------------------------------------- *)
(* 1. Helicity <-> local real component map                          *)
(* ----------------------------------------------------------------- *)
sigPlus  = (P + Pb)/2;           (* sigma_+ = Re Psi0 *)
sigCross = (P - Pb)/(2 I);       (* sigma_x = Im Psi0 *)

(* Expectation as a linear functional: expand a field product in
   {Phi, P, Pb} and map each monomial Phi^a P^p Pb^q to the helicity
   cumulant cum[a,p,q] = <Phi^a P^p Pb^q>_c. *)
cumOf[expr_] := Module[{cr},
  cr = CoefficientRules[Expand[expr], {Phi, P, Pb}];
  Total[(#[[2]] (cum @@ #[[1]])) & /@ cr]
];

(* Parity invariance about the separation axis swaps P<->Pb, so
   cum[a,p,q] = cum[a,q,p] and both are real. Impose by mapping each
   helicity cumulant to a canonical symmetric symbol csym[a,{min,max}]. *)
parity = cum[a_, p_, q_] :> csym[a, Sort[{p, q}]];

red[expr_] := Simplify[cumOf[expr] /. parity];

(* ----------------------------------------------------------------- *)
(* 2. The local real-component cumulant tensor entries               *)
(* ----------------------------------------------------------------- *)
K000 = red[Phi^3];
K001 = red[Phi^2 sigPlus];
K002 = red[Phi^2 sigCross];
K011 = red[Phi sigPlus^2];
K012 = red[Phi sigPlus sigCross];
K022 = red[Phi sigCross^2];
K111 = red[sigPlus^3];
K112 = red[sigPlus^2 sigCross];
K122 = red[sigPlus sigCross^2];
K222 = red[sigCross^3];

Print["== Local real-basis cumulant entries (parity-reduced) =="];
Print["  K000 (TTT)  = ", K000];
Print["  K001 (TTσ+) = ", K001];
Print["  K002 (TTσx) = ", K002, "    (must vanish: one sigma_x, parity-odd)"];
Print["  K011 (Tσ+σ+)= ", K011];
Print["  K012        = ", K012, "    (must vanish)"];
Print["  K022 (Tσxσx)= ", K022];
Print["  K111 (σ+^3) = ", K111];
Print["  K112        = ", K112, "    (must vanish)"];
Print["  K122 (σ+σxσx)=", K122];
Print["  K222 (σx^3) = ", K222, "    (must vanish)"];

(* ----------------------------------------------------------------- *)
(* 3. Name the helicity channels                                     *)
(* ----------------------------------------------------------------- *)
(* canoes-COMPUTED channels (all-+2 helicity, then Re):
     zetaTTP = Re<Phi^2 P>      = csym[2,{0,1}]
     zetaTPP = Re<Phi   P P>    = csym[1,{0,2}]
     zetaPPP = Re<       P P P> = csym[0,{0,3}]
   MISSING conjugate-paired (modulus) channels:
     Bmod = <Phi P Pb>   = csym[1,{1,1}]   (spin weights (0,2,-2))
     Dmod = <    P P Pb> = csym[0,{1,2}]   (spin weights (2,2,-2))            *)
zetaTTP = csym[2, {0, 1}];
zetaTPP = csym[1, {0, 2}];
zetaPPP = csym[0, {0, 3}];
Bmod    = csym[1, {1, 1}];
Dmod    = csym[0, {1, 2}];

Print["\n== Channel identification =="];
Print["  zetaTTP = ", zetaTTP, "   zetaTPP = ", zetaTPP, "   zetaPPP = ", zetaPPP];
Print["  Bmod = <Phi P Pb> = ", Bmod, "   (MISSING; spin (0,2,-2))"];
Print["  Dmod = <P P Pb>   = ", Dmod, "   (MISSING; spin (2,2,-2))"];

(* ----------------------------------------------------------------- *)
(* 4. Express the needed tensor entries in computable channels       *)
(* ----------------------------------------------------------------- *)
relK001 = Simplify[K001 - zetaTTP];        (* should be 0: TTP needs no fix *)
relK011 = Simplify[K011 - (zetaTPP + Bmod)/2];
relK022 = Simplify[K022 - (Bmod - zetaTPP)/2];
relK111 = Simplify[K111 - (zetaPPP + 3 Dmod)/4];
relK122 = Simplify[K122 - (Dmod - zetaPPP)/4];

check["K001 = zetaTTP (TTP unaffected)", relK001 === 0];
check["K011 = (zetaTPP + Bmod)/2",        relK011 === 0];
check["K022 = (Bmod - zetaTPP)/2",        relK022 === 0];
check["K111 = (zetaPPP + 3 Dmod)/4",      relK111 === 0];
check["K122 = (Dmod - zetaPPP)/4",        relK122 === 0];

(* ----------------------------------------------------------------- *)
(* 5. F-vertex consistency: the FK <kappa kappa> contraction needs   *)
(*    K011 + K022 = <Phi |Psi0|^2> = Bmod (the modulus channel).     *)
(* ----------------------------------------------------------------- *)
sumCheck = Simplify[(K011 + K022) - Bmod];
diffCheck = Simplify[(K011 - K022) - zetaTPP];
check["K011 + K022 = Bmod  (F-vertex needs the SUM = modulus)", sumCheck === 0];
check["K011 - K022 = zetaTPP  (what canoes currently computes)", diffCheck === 0];

(* The current callable error: it sets out[0,1,1]=zetaTPP, out[0,2,2]=0.
   The FK kappa-kappa contraction F[0,1,1] K[0,1,1] + F[0,2,2] K[0,2,2]
   with F[0,1,1]=F[0,2,2]=-1 then yields -(zetaTPP + 0) instead of
   -(K011 + K022) = -Bmod.  Error = correct - wrong = -(Bmod) + zetaTPP
   = -(Bmod - zetaTPP) = -2 K022. *)
errTPP = Simplify[(-(K011 + K022)) - (-(zetaTPP + 0))];
check["FK(kk) TPP-sector error = -2 K022 (order-1)", Simplify[errTPP - (-2 K022)] === 0];
Print["\n  FK <kk> TPP-sector error of current code = ", errTPP, "  ( = -2*K022 )"];

(* Parity selection: odd-sigma_x entries vanish identically. *)
check["K002 = 0 (parity)", K002 === 0];
check["K012 = 0 (parity)", K012 === 0];
check["K112 = 0 (parity)", K112 === 0];
check["K222 = 0 (parity)", K222 === 0];

(* ----------------------------------------------------------------- *)
(* 6. Small-angle suppression of the all-+2 channel (Wigner-d)       *)
(*    d^L_{2,-2}(theta) ~ theta^4  vs  d^L_{2,2}(theta) ~ 1 .         *)
(*    The all-+2 cumulant <P P> (= zetaTPP) projects on d_{2,-2};     *)
(*    the modulus <P Pb> (= Bmod) projects on d_{2,2}.                *)
(* ----------------------------------------------------------------- *)
dmm[L_, th_] := WignerD[{L, 2, -2}, 0, th, 0];   (* spin-flip: ~ th^4 *)
dpp[L_, th_] := WignerD[{L, 2,  2}, 0, th, 0];   (* aligned:   ~ 1    *)
Lt = 10;
serFlip  = Series[dmm[Lt, th], {th, 0, 5}] // Normal;
serAlign = Series[dpp[Lt, th], {th, 0, 2}] // Normal;
ordFlip  = Exponent[serFlip /. th -> th, th, Min];   (* leading power *)
leadAlign = Limit[dpp[Lt, th], th -> 0];
Print["\n== Small-angle Wigner-d (L=", Lt, ") =="];
Print["  d^L_{2,-2}(th) leading power in th = ", ordFlip,
      "   (all-+2 / zetaTPP channel -> gamma^4 suppressed)"];
Print["  d^L_{2, 2}(0) = ", leadAlign, "   (modulus / Bmod channel -> O(1))"];
check["all-+2 channel d_{2,-2} ~ theta^4 (suppressed)", ordFlip == 4];
check["modulus channel d_{2,2}(0) = 1 (unsuppressed)", leadAlign == 1];

(* ----------------------------------------------------------------- *)
(* 7. Fix recipe (printed for the canoes/callable patch)             *)
(* ----------------------------------------------------------------- *)
Print["\n== FIX RECIPE =="];
Print["  Add canoes channels with conjugate helicity:"];
Print["    Bmod  via spin_weights (0, 2, -2)  = <Phi00 Psi0 conj(Psi0)>"];
Print["    Dmod  via spin_weights (2, 2, -2)  = <Psi0 Psi0 conj(Psi0)>"];
Print["  Then the callable fills the LOCAL real-basis tensor as:"];
Print["    out[0,0,0]               = zetaTTT"];
Print["    out[0,0,1]+perm          = zetaTTP"];
Print["    out[0,1,1]+perm          = (zetaTPP + Bmod)/2     = K011"];
Print["    out[0,2,2]+perm          = (Bmod - zetaTPP)/2     = K022   (NEW)"];
Print["    out[1,1,1]               = (zetaPPP + 3 Dmod)/4   = K111"];
Print["    out[1,2,2]+perm          = (Dmod - zetaPPP)/4     = K122   (NEW)"];
Print["    out[*,*,2] with odd #2   = 0   (parity)"];

(* ----------------------------------------------------------------- *)
(* 8. PASS/FAIL summary                                              *)
(* ----------------------------------------------------------------- *)
Print["\n=================== SUMMARY ==================="];
Print["PASS (", Length[pass], "): ", pass];
If[Length[fail] > 0,
  Print["FAIL (", Length[fail], "): ", fail],
  Print["FAIL (0): {}   --  ALL CHECKS PASS"]
];
