(* ============================================================================
   derive_spin2_gaunt_norm.wl

   Pin the CORRECT rotational-invariant reduced->full anchor for a spin-weighted
   three-point cumulant, and quantify the spurious factor the deployed
   SpinAwareZetaEvaluator carries by using the (0,0,0) anchor for spin channels.

   STRUCTURE OF THE CLAIM
   ----------------------
   The spin-weighted-spherical-harmonic Gaunt integral is

     INT dOmega  {s1}Y_{L1 M1} {s2}Y_{L2 M2} {s3}Y_{L3 M3}
        = sqrt[(2L1+1)(2L2+1)(2L3+1)/(4 pi)]
          ( L1 L2 L3 ; -s1 -s2 -s3 ) ( L1 L2 L3 ; M1 M2 M3 ).

   So a rotationally invariant three-point cumulant of spin-(s1,s2,s3) fields,
   built from a reduced (isotropic) bispectrum b_{L1L2L3}, carries the geometric
   weight

     h_{L1L2L3} * ( L1 L2 L3 ; -s1 -s2 -s3 )    [the SPIN anchor]

   where h = sqrt[(2L1+1)(2L2+1)(2L3+1)/4pi].  For the SCALAR case
   (s=0,0,0) this is the (0,0,0) symbol -- which is what the scalar
   ZetaEvaluator correctly uses.  The deployed SpinAwareZetaEvaluator uses
   the SAME (0,0,0) symbol for the spin channels (kernel build line 604:
   coef = h*w0*sq_l1*(-1)^s1, w0 = wigner_3j(...;0,0,0)).  The spurious factor
   is therefore

     R(L1,L2,L3; s) = (L1 L2 L3 ; 0,0,0) / (L1 L2 L3 ; -s1,-s2,-s3).

   We tabulate R for the modulus channels D=(2,2,-2), B=(0,2,-2) and confirm
   it is ell-dependent (=> angle-structured after the b-sum), and check the
   small-angle anchor sqrt(L^2(L^2-2)) -> L^2.

   PASS/FAIL at the bottom.
*)

ClearAll["Global`*"];
pass = {}; fail = {};
check[name_, cond_] := If[TrueQ[cond], AppendTo[pass, name], AppendTo[fail, name]];

w3 = ThreeJSymbol;

(* note: ThreeJSymbol[{l1,m1},{l2,m2},{l3,m3}] *)
w0[l1_,l2_,l3_] := w3[{l1,0},{l2,0},{l3,0}];
ws[l1_,l2_,l3_,s1_,s2_,s3_] := w3[{l1,-s1},{l2,-s2},{l3,-s3}];

(* ---- modulus channel D = (2,2,-2): -s = (-2,-2,2) ---- *)
Print["%%SECTION%% R(L,L,L) = (000)/(-s) for D=(2,2,-2) [so -s=(-2,-2,2)]"];
Llist = {2,4,6,10,20,40,80,160};
Do[
  Module[{a = N[w0[L,L,L]], b = N[ws[L,L,L,2,2,-2]]},
   Print["  L=", L, "  (000)=", a, "  (-2,-2,2)=", b,
     "  R=", If[b == 0, "Inf", a/b]]],
  {L, Llist}];

Print["%%SECTION%% R(L,L,L) for B=(0,2,-2) [so -s=(0,-2,2)]"];
Do[
  Module[{a = N[w0[L,L,L]], b = N[ws[L,L,L,0,2,-2]]},
   Print["  L=", L, "  (000)=", a, "  (0,-2,2)=", b,
     "  R=", If[b == 0, "Inf", a/b]]],
  {L, Llist}];

(* ---- asymptotic R at large L on the diagonal ---- *)
Print["%%SECTION%% large-L asymptotics of R (numeric)"];
Module[{Lbig = 400, a, b},
  a = N[w0[Lbig,Lbig,Lbig]]; b = N[ws[Lbig,Lbig,Lbig,2,2,-2]];
  Print["  D channel L=400  R=", a/b];
];

(* ---- The small-angle / high-ell d-function anchor from the appendix ----
   The MODULUS spin-2 leg carries d^L_{2,2}(theta) -> 1; the projection
   constant relating |Psi0|^2 to |kappa|^2 in the small-angle limit is
   sqrt(L^2 (L^2 - 2)) -> L^2.  Confirm symbolically. *)
Print["%%SECTION%% appendix small-angle limit sqrt(L^2(L^2-2)) -> L^2"];
limFac = Limit[Sqrt[L^2 (L^2 - 2)]/L^2, L -> Infinity];
Print["  lim_{L->inf} sqrt(L^2(L^2-2))/L^2 = ", limFac];
check["sqrt(L^2(L^2-2))/L^2 -> 1", limFac === 1];

dpp[L_, th_] := WignerD[{L, 2, 2}, 0, th, 0];
dmm[L_, th_] := WignerD[{L, 2, -2}, 0, th, 0];
Print["  d^10_{2,2}(0) = ", Limit[dpp[10, th], th -> 0], " (modulus, O(1))"];
Print["  d^10_{2,-2}(th) leading power = ",
  Exponent[Normal[Series[dmm[10, th], {th, 0, 5}]], th, Min],
  " (un-conjugated, gamma^4-suppressed)"];
check["d_{2,2}(0)=1", Limit[dpp[10, th], th -> 0] === 1];
check["d_{2,-2}~th^4", Exponent[Normal[Series[dmm[10, th], {th, 0, 5}]], th, Min] == 4];

(* ---- VERDICT: is R ~ few (3.5-5x) and ell-dependent? ---- *)
Print["%%SECTION%% verdict numbers"];
rD = Table[Module[{a = N[w0[L,L,L]], b = N[ws[L,L,L,2,2,-2]]},
   If[b == 0, Missing[], a/b]], {L, {6,10,20,40,80,160}}];
rD = DeleteCases[rD, _Missing];
Print["  R(D) over L in {6,10,20,40,80,160} = ", rD];
Print["  R(D) range = [", Min[rD], ", ", Max[rD], "]  spread=", Max[rD]/Min[rD]];

Print["\n=================== SUMMARY ==================="];
Print["PASS (", Length[pass], "): ", pass];
If[Length[fail] > 0, Print["FAIL (", Length[fail], "): ", fail],
  Print["FAIL (0): {} -- ALL CHECKS PASS"]];
Print["%%DONE%%"];
