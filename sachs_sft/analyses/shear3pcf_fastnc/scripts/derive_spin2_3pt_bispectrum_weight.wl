(* derive_spin2_3pt_bispectrum_weight.wl
   ----------------------------------------------------------------------------
   Phase-1 symbolic derivation for the canoes spin-2 zeta_D HIGH/Limber fix.

   GOAL
   ----
   Determine the CORRECT per-leg ell-weighting that converts the scalar
   convergence (kappa) response into the spin-2 shear (gamma) response in a
   3-point reduced bispectrum, in the flat-sky / high-ell limit, and check
   whether the canoes HIGH Born/sub-Hubble response

       Phi00(l) = +A(a)(1+z)^4 delta(k),    Psi0(l) = -Phi00(l)   (ell-FLAT)

   reproduces it or DROPS a per-leg factor.

   PHYSICS
   -------
   On the sphere, the lensing potential psi_lm sources:
     kappa_lm  = -(1/2) Laplacian psi      eigenvalue  L2 = ell(ell+1)
                  => kappa_lm = +(1/2) L2 psi_lm
     gamma_lm  = (1/2) eth^2 psi  (spin +2)  -- spin-raising eth applied twice
                  the spin-weighted Laplacian eigenvalue of eth^2 on a spin-0
                  field is  sqrt( (ell-1)(ell+2) * ell(ell+1) ) = sqrt(L2(L2-2))
                  => |gamma_lm| = (1/2) sqrt(L2(L2-2)) |psi_lm|
   so the EXACT per-leg spin-2/scalar modulus ratio is
       W2(ell) = |gamma_lm|/|kappa_lm| = sqrt(L2(L2-2)) / L2 ,   L2 = ell(ell+1).

   This is exactly canoes' COEFS_PSI0 ID coefficient
       -sqrt(L2(L2-2))/chi^2     vs   COEFS_PHI00 ID  L2/chi^2 (high-ell lead).

   The eth-raising eigenvalue identity (Newman-Penrose / spin-weighted sphere):
     eth   {}_sY_lm  = +sqrt( (ell-s)(ell+s+1) ) {}_{s+1}Y_lm
     ethbar{}_sY_lm  = -sqrt( (ell+s)(ell-s+1) ) {}_{s-1}Y_lm
   so eth^2 on s=0 :
     eth^2 {}_0Y_lm = sqrt((ell-0)(ell+1)) * sqrt((ell-1)(ell+2)) {}_2Y_lm
                    = sqrt( ell(ell+1)(ell-1)(ell+2) ) {}_2Y_lm
                    = sqrt( L2 (L2 - 2) ) {}_2Y_lm           [since (ell-1)(ell+2)=L2-2]
   We verify (ell-1)ell(ell+1)(ell+2) = L2(L2-2) symbolically and the high-ell
   limit W2(ell) -> 1, and we EXPAND W2 - 1 to show its leading correction is
   O(1/ell^2) (NOT a power-law growth), which is the crux of the diagnosis.
   ---------------------------------------------------------------------------- *)

Print["%%RUN_START%%"];

(* ---- 1. eth-raising eigenvalues (closed form, no package needed) ---- *)
ethEig[ell_, s_]    := Sqrt[(ell - s)(ell + s + 1)];      (* eth : s -> s+1 *)
ethbarEig[ell_, s_] := -Sqrt[(ell + s)(ell - s + 1)];     (* ethbar : s -> s-1 *)

(* eth^2 acting on spin-0 field: product of the two raising eigenvalues *)
eth2Eig = ethEig[ell, 0] * ethEig[ell, 1];
eth2EigSimp = Simplify[eth2Eig, Assumptions -> ell >= 2];

Print["eth^2 (s=0->2) eigenvalue  = ", ToString[InputForm[eth2EigSimp]]];

(* ---- 2. Identity:  (ell-1) ell (ell+1) (ell+2) == L2 (L2 - 2) ---- *)
L2 = ell (ell + 1);
lhs = (ell - 1) ell (ell + 1) (ell + 2);
rhs = L2 (L2 - 2);
identityZero = Simplify[lhs - rhs];
Print["[T1] (ell-1)ell(ell+1)(ell+2) - L2(L2-2) = ", ToString[InputForm[identityZero]],
      "   -> ", If[identityZero === 0, "PASS", "FAIL"]];

(* eth^2 eigenvalue equals sqrt(L2(L2-2)) *)
eth2VsSqrt = Simplify[eth2EigSimp - Sqrt[rhs], Assumptions -> ell >= 2];
Print["[T2] eth^2eig - sqrt(L2(L2-2)) = ", ToString[InputForm[eth2VsSqrt]],
      "   -> ", If[eth2VsSqrt === 0, "PASS", "FAIL"]];

(* ---- 3. Per-leg spin-2/scalar modulus ratio and high-ell limit ---- *)
W2[ell_] := Sqrt[L2 (L2 - 2)] / L2 /. {L2 -> ell (ell + 1)};
W2sym = Sqrt[(ell (ell + 1)) ((ell (ell + 1)) - 2)] / (ell (ell + 1));

limHigh = Limit[W2sym, ell -> Infinity];
Print["[T3] lim_{ell->inf} W2(ell) = ", ToString[InputForm[limHigh]],
      "   -> ", If[limHigh === 1, "PASS (per-leg spin-2 weight -> 1)", "FAIL"]];

(* leading correction: series of W2 - 1 at large ell *)
ser = Series[W2sym - 1, {ell, Infinity, 4}] // Normal;
Print["[T4] W2(ell) - 1  (large-ell series) = ", ToString[InputForm[Simplify[ser]]]];
Print["     => leading correction is O(1/ell^2); W2 is FLAT-to-1, NOT power-law."];

(* numeric table of W2 over the (60,1000] window used by canoes HIGH *)
Print["%%W2_TABLE_START%%"];
Print["ell      W2(ell)=sqrt(L2(L2-2))/L2     1-W2"];
Do[
  wv = N[W2sym /. ell -> elv, 12];
  Print[ToString[elv], "    ", ToString[wv], "    ", ToString[N[1 - wv, 6]]],
  {elv, {60, 80, 100, 150, 200, 300, 500, 800, 1000}}
];
Print["%%W2_TABLE_END%%"];

(* ---- 4. Three-leg product weight (the 3-point reduced bispectrum) ---- *)
(* zeta_D / zeta_TTT (psi0^3 / phi00^3) per-leg-product spin-2/scalar weight,
   equilateral l1=l2=l3=L: should be W2(L)^3 -> 1 at high ell. *)
W2prod3[L_] := W2sym^3 /. ell -> L;
Print["[T5] equilateral 3-leg spin-2/scalar weight  W2(L)^3:"];
Print["     lim_{L->inf} = ", ToString[InputForm[Limit[W2sym^3, ell -> Infinity]]],
      "   -> ", If[Limit[W2sym^3, ell -> Infinity] === 1, "PASS", "FAIL"]];
Print["%%W2CUBE_TABLE_START%%"];
Print["L       W2(L)^3        1-W2(L)^3"];
Do[
  wv3 = N[W2sym^3 /. ell -> Lv, 12];
  Print[ToString[Lv], "    ", ToString[wv3], "    ", ToString[N[1 - wv3, 6]]],
  {Lv, {60, 80, 100, 150, 200, 300, 500, 800, 1000}}
];
Print["%%W2CUBE_TABLE_END%%"];

(* ---- 5. VERDICT on the canoes HIGH ell-flat response ---- *)
(* canoes HIGH: Psi0(l) = -Phi00(l), i.e. per-leg spin-2/scalar weight = 1 EXACTLY.
   Correct per-leg weight = W2(ell) = sqrt(L2(L2-2))/L2.
   The RATIO canoes/correct per leg is:  1 / W2(ell) = L2/sqrt(L2(L2-2)).      *)
canoesOverCorrect[ell_] := 1 / W2sym;
Print["%%VERDICT_START%%"];
Print["canoes-HIGH per-leg weight = 1 (ell-flat, Psi0=-Phi00)."];
Print["correct  per-leg weight    = W2(ell) = sqrt(L2(L2-2))/L2  (-> 1 high-ell)."];
Print["per-leg canoes/correct     = 1/W2(ell) = L2/sqrt(L2(L2-2)) (-> 1 high-ell)."];
Print["3-leg canoes/correct       = 1/W2(ell)^3 (-> 1 high-ell)."];
Print["Series of 1/W2 - 1 (large ell) = ",
      ToString[InputForm[Simplify[Normal[Series[1/W2sym - 1, {ell, Infinity, 4}]]]]]];
Print["=> At HIGH ell the canoes ell-flat response is CORRECT to O(1/ell^2)."];
Print["=> Dropping W2 CANNOT produce a +0.44 log-log slope or a 4-17x at ell~100-800."];
Print["=> Therefore the per-leg spin-2 eth^2 weight is NOT the +0.44 lever."];
Print["%%VERDICT_END%%"];

Print["%%RUN_END%%"];
