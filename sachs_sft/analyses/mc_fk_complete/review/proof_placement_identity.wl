(* ::Package:: *)

(* =====================================================================
   proof_placement_identity.wl

   SYMBOLIC PROOF (generic symbols, no numerics anywhere) of the
   "Q-placement" identity for the FK cross term:

        (T1 + T2) + (T1 + T2)^T  ==  TRUTH   ==   M + M^T

   where
        T1_AB = Cov( kappa2_A , kappa_d1_B )
        T2_AB = Cov( kappa1_A , kappa_d2_B )
   and TRUTH is built INDEPENDENTLY of T1/T2 from the induced O(Q)
   three-point function of the deformed driving field,
        Xi_{bcB}[m,n,p] = < f_b[m] f_c[n] f_B[p] > |_{O(Q)}
   contracted through the FK kernels.

   Everything (Wick contractions included) is done by one generic
   Isserlis engine on purely symbolic Gaussian covariances.  No numbers,
   no fits, no special values.

   Discrete model (see SPEC.md, section 1-2), with ncomp components and
   nnode nodes k = 0 .. nnode-1:

     z[k]                  Gaussian, zero mean, < z[k] z[l]^T > = R[k,l]
                           with R[k,l] = Transpose[R[l,k]] (imposed), all
                           independent entries generic symbols rr[...]

     df[k]_a = 1/2 Q[k]_{a,mn} ( z_m[k] z_n[k] - C[k]_{mn} ),  C[k] = R[k,k]
                           Q[k]_{a,mn} generic, symmetric in (m,n)

     f[k]    = z[k] + df[k]                       (deformed drive)

     s1[k]   = resp[k] s1[k-1] + z[k]  dlam ,     s1[-1]  = 0
     sd1[k]  = resp[k] sd1[k-1] + df[k] dlam ,    sd1[-1] = 0

     c_j     = dlam * Wd[j]      (generic weights ww[j])
     resp[k] = generic weights rho[k]
     F_{abc} = generic ff[a,b,c]   (NO symmetry assumed)

     kappa1_A   = Sum_j c_j z_A[j]
     kappa2_A   = Sum_j c_j F_{Abc} s1_b[j-1] s1_c[j-1]
     kappa_d1_B = Sum_p c_p df_B[p]
     kappa_d2_B = Sum_j c_j (F_{Bbc} + F_{Bcb}) s1_b[j-1] sd1_c[j-1]

   All weights are generic symbols; the claimed identity is a polynomial
   identity in them, so "generic symbol" is strictly stronger than
   "generic positive number".

   ---------------------------------------------------------------------
   WHY IT HOLDS (index-level argument; the script below verifies every
   step of it, and the leg-by-leg checks "placement: ..." are the machine
   statement of the two bullet points).

   Write G[k,m] for the linear response of s1 at node k to z at node m
   (G[k,m] = prod_{i=m+1}^k resp[i], zero for m > k), so

       s1_b[j-1] = dlam Sum_m G[j-1,m] z_b[m].

   Because the drive is a *quadratic* deformation of a Gaussian, one Wick
   contraction of <df_a[k] z_c[n] z_B[p]> gives two equal pairings which,
   by the (m,n)-symmetry of Q, fold into the single term
   Q[k]_{a,uv} R[k,n]_{uc} R[k,p]_{vB} with coefficient 1; the disconnected
   pairing is cancelled exactly by the C[k] = R[k,k] subtraction.  Hence

       Xi_{bcB}[m,n,p] = (leg-1 term) + (leg-2 term) + (leg-3 term)

   with the Q sitting on m, on n, or on p respectively.  Folding Xi through
   the FK kernels (F vertex + two G legs on the kappa2 side, one c_q leg on
   the kappa1 side) therefore splits M_AB = <kappa2_A kappa1_B>_Xi into
   exactly three pieces M1 + M2 + M3, one per Q placement.  Then:

   * M3 == T1.  The third placement puts Q on the leg that is contracted
     with kappa1_B.  Replacing z_B[q] by df_B[q] in kappa1_B is literally
     kappa_d1_B, so this placement is generated verbatim by
     T1 = Cov(kappa2_A, kappa_d1_B).  Nothing is symmetrised here.

   * M1 + M2 == T2^T (i.e. M1_AB + M2_AB = T2_BA).  The first two placements
     put Q on one of the two F-vertex input legs.  The pair
     (b,m) <-> (c,n) is a relabelling symmetry of s1_b[j-1] s1_c[j-1], so
     the two placements collapse onto ONE placement carrying the combination
     F_{Abc} + F_{Acb}.  That is exactly kappa_d2 with one s1 leg promoted to
     sd1, and the free index A now sits on the F vertex while B sits on the
     bare linear leg -- i.e. it is T2 with A and B exchanged, hence T2^T.

   Therefore  M = T1 + T2^T  identically, and
       TRUTH = M + M^T = (T1 + T2^T) + (T1^T + T2) = (T1+T2) + (T1+T2)^T.

   The argument uses only: (i) Q symmetric in its last two indices,
   (ii) C[k] = R[k,k], (iii) the two F-vertex legs carrying the SAME time
   argument j-1 and the SAME kernel, (iv) kappa_d1 / kappa_d2 reusing the
   kernels of kappa1 / kappa2.  It never uses the values of n_comp or N, the
   number of components appearing only as a dummy summation range and the
   number of nodes only through the range of the m,n,p,j,q sums.  So the
   identity is (n_comp, N)-independent; the explicit cases below are a
   machine check of the bookkeeping, not the source of its generality.
   ===================================================================== *)

ClearAll["Global`*"];
$HistoryLength = 0;

(* ------------------------------------------------------------------ *)
(* 1.  Generic Gaussian covariance with the only structural constraint *)
(*     R[k,l] = Transpose[R[l,k]]  (and hence R[k,k] symmetric).       *)
(* ------------------------------------------------------------------ *)

rE[k_Integer, l_Integer, m_Integer, n_Integer] := Which[
   k < l, rr[k, l, m, n],
   k > l, rr[l, k, n, m],
   True , rr[k, k, Min[m, n], Max[m, n]]
];

(* ------------------------------------------------------------------ *)
(* 2.  Generic Isserlis / Wick engine.                                 *)
(*     zz[k,a] are the formal Gaussian variables.  EW[.] is the exact  *)
(*     Gaussian expectation of an arbitrary polynomial in them.        *)
(* ------------------------------------------------------------------ *)

cov[zz[k_, m_], zz[l_, n_]] := rE[k, l, m, n];

isser[{}] = 1;
isser[l_List] := isser[l] =
  If[OddQ[Length[l]], 0,
   Module[{a = First[l], rest = Rest[l]},
    Sum[cov[a, rest[[i]]] isser[Delete[rest, i]], {i, Length[rest]}]]];

wickMono[t_] := Module[{facs, zf = {}, co = 1},
  facs = If[Head[t] === Times, List @@ t, {t}];
  Do[
   Which[
    MatchQ[f, zz[_, _]],
      AppendTo[zf, f],
    MatchQ[f, Power[zz[_, _], _Integer?Positive]],
      Do[AppendTo[zf, f[[1]]], {f[[2]]}],
    True,
      co = co f
    ],
   {f, facs}];
  co isser[zf]];

EW[expr_] := Module[{ex = Expand[expr]},
  If[ex === 0, 0,
   Total[wickMono /@ If[Head[ex] === Plus, List @@ ex, {ex}]]]];

(* --- engine self-tests (generic, symbolic) ------------------------- *)
engineSelfTest[] := Module[{ok = True, e1, e2, e3, e4},
  e1 = Expand[EW[zz[0, 1] zz[1, 2]] - rE[0, 1, 1, 2]];
  e2 = Expand[EW[zz[1, 2] zz[0, 1]] - rE[0, 1, 1, 2]];       (* symmetry *)
  e3 = Expand[EW[zz[0, 1] zz[1, 1] zz[2, 2]]];               (* odd -> 0 *)
  e4 = Expand[EW[zz[0,1] zz[1,2] zz[2,1] zz[0,2]]
        - ( rE[0,1,1,2] rE[2,0,1,2] + rE[0,2,1,1] rE[1,0,2,2]
          + rE[0,0,1,2] rE[1,2,2,1] )];
  ok = (e1 === 0) && (e2 === 0) && (e3 === 0) && (e4 === 0);
  Print["  Wick engine self-test (2pt, 2pt-sym, 3pt-odd, 4pt-Isserlis): ",
        If[ok, "PASS", "FAIL " <> ToString[{e1, e2, e3, e4}]]];
  ok];

(* ------------------------------------------------------------------ *)
(* 3.  One full case (ncomp, nnode).                                   *)
(* ------------------------------------------------------------------ *)

runCase[nc_Integer, nd_Integer] := Module[
  {qE, dfv, s1, sdv, cw, gcT, k1, k2, kd1, kd2, k2alt, chkKer, chkComp,
   T1, T2, Xi, XiL, XiM, XiR, fold, M, M1, M2, M3, lhs, truth,
   pl1, pl2, pl3, pl4, kd2bad, T2bad, lhsbad, ctlA, kd1bad, T1bad, ctlB,
   d0, d1, dA, dB, dC, nzTruth, nTerms, t0},

  Print["==================================================================="];
  Print["CASE  (n_comp, N) = (", nc, ", ", nd, ")"];
  Print["==================================================================="];
  t0 = AbsoluteTime[];

  (* Q[k]_{a,mn}, generic, symmetric in (m,n) *)
  qE[k_, a_, m_, n_] := qq[k, a, Min[m, n], Max[m, n]];

  (* df[k]_a  ---- array index: dfv[[k+1, a]] is node k *)
  dfv = Table[
    (1/2) Sum[qE[k, a, u, v] (zz[k, u] zz[k, v] - rE[k, k, u, v]),
       {u, 1, nc}, {v, 1, nc}],
    {k, 0, nd - 1}, {a, 1, nc}];

  (* s1 and sd1 by the stated recursions.
     array index: s1[[k+2, a]] is node k ; s1[[1, a]] is node -1 (= 0). *)
  s1  = ConstantArray[0, {nd + 1, nc}];
  sdv = ConstantArray[0, {nd + 1, nc}];
  Do[
   Do[
    s1[[k + 2, a]]  = Expand[rho[k] s1[[k + 1, a]]  + zz[k, a] dl];
    sdv[[k + 2, a]] = Expand[rho[k] sdv[[k + 1, a]] + dfv[[k + 1, a]] dl];
    , {a, 1, nc}],
   {k, 0, nd - 1}];

  cw[j_] := dl ww[j];

  (* FK linear kernel read off the s1 recursion (NOT assumed in closed
     form): gcT[[kk, m+1]] = d s1[node kk-2]_a / d z_a[m]                *)
  gcT = Table[Coefficient[s1[[kk, 1]], zz[m, 1]], {kk, 1, nd + 1}, {m, 0, nd - 1}];

  (* sanity: s1 is component-diagonal and component-independent *)
  chkComp = Expand[
    Flatten[Table[
      {Coefficient[s1[[kk, a]], zz[m, a]] - gcT[[kk, m + 1]],
       If[a =!= b, Coefficient[s1[[kk, a]], zz[m, b]], 0]},
      {kk, 1, nd + 1}, {m, 0, nd - 1}, {a, 1, nc}, {b, 1, nc}]]];
  Print["  kernel extraction consistent (s1 component-diagonal): ",
        If[Union[chkComp] === {0}, "PASS", "FAIL"]];

  (* ---------------- the four FK arms ----------------------------- *)
  k1 = Table[Sum[cw[j] zz[j, A], {j, 0, nd - 1}], {A, 1, nc}];

  k2 = Table[
    Sum[cw[j] Sum[ff[A, b, c] s1[[j + 1, b]] s1[[j + 1, c]],
        {b, 1, nc}, {c, 1, nc}], {j, 0, nd - 1}],
    {A, 1, nc}];

  kd1 = Table[Sum[cw[p] dfv[[p + 1, B]], {p, 0, nd - 1}], {B, 1, nc}];

  kd2 = Table[
    Sum[cw[j] Sum[(ff[B, b, c] + ff[B, c, b]) s1[[j + 1, b]] sdv[[j + 1, c]],
        {b, 1, nc}, {c, 1, nc}], {j, 0, nd - 1}],
    {B, 1, nc}];

  (* kappa2 rewritten purely through the extracted kernels -- used below
     to build TRUTH.  Check it reproduces the recursion-built kappa2. *)
  k2alt = Table[
    Sum[cw[j] Sum[ff[A, b, c] Sum[
        gcT[[j + 1, m + 1]] gcT[[j + 1, n + 1]] zz[m, b] zz[n, c],
        {m, 0, nd - 1}, {n, 0, nd - 1}], {b, 1, nc}, {c, 1, nc}],
      {j, 0, nd - 1}],
    {A, 1, nc}];
  chkKer = Expand[k2 - k2alt];
  Print["  kappa2 via kernels == kappa2 via recursion: ",
        If[Union[Flatten[{chkKer}]] === {0}, "PASS", "FAIL"]];

  (* ---------------- T1, T2 by Wick ------------------------------- *)
  T1 = Table[Expand[EW[k2[[A]] kd1[[B]]] - EW[k2[[A]]] EW[kd1[[B]]]],
     {A, 1, nc}, {B, 1, nc}];
  T2 = Table[Expand[EW[k1[[A]] kd2[[B]]] - EW[k1[[A]]] EW[kd2[[B]]]],
     {A, 1, nc}, {B, 1, nc}];
  Print["  T1, T2 built.  t = ", Round[AbsoluteTime[] - t0, 0.01], " s"];

  (* ---------------- TRUTH, built without reference to T1/T2 ------- *)
  (* Xi_{bcB}[m,n,p] = <f_b[m] f_c[n] f_B[p]>|_{O(Q)} : exactly the three
     terms in which one leg carries the Q, each evaluated by the same
     generic Wick engine (so the C = R[k,k] subtraction is handled, not
     assumed).  Kept SEPARATE so that the placement bookkeeping can be
     verified leg by leg.                                             *)
  XiL = Table[EW[dfv[[m + 1, b]] zz[n, c] zz[p, B]],
    {b, 1, nc}, {c, 1, nc}, {B, 1, nc},
    {m, 0, nd - 1}, {n, 0, nd - 1}, {p, 0, nd - 1}];
  XiM = Table[EW[zz[m, b] dfv[[n + 1, c]] zz[p, B]],
    {b, 1, nc}, {c, 1, nc}, {B, 1, nc},
    {m, 0, nd - 1}, {n, 0, nd - 1}, {p, 0, nd - 1}];
  XiR = Table[EW[zz[m, b] zz[n, c] dfv[[p + 1, B]]],
    {b, 1, nc}, {c, 1, nc}, {B, 1, nc},
    {m, 0, nd - 1}, {n, 0, nd - 1}, {p, 0, nd - 1}];
  Xi = XiL + XiM + XiR;
  Print["  Xi built.        t = ", Round[AbsoluteTime[] - t0, 0.01], " s"];

  (* fold[X]_AB = <kappa2_A kappa1_B>_X : contract through the FK kernels *)
  fold[X_] := Table[
    Expand[Sum[cw[j] Sum[ff[A, b, c] Sum[
        gcT[[j + 1, m + 1]] gcT[[j + 1, n + 1]]
          Sum[cw[q] X[[b, c, B, m + 1, n + 1, q + 1]], {q, 0, nd - 1}],
        {m, 0, nd - 1}, {n, 0, nd - 1}], {b, 1, nc}, {c, 1, nc}],
      {j, 0, nd - 1}]],
    {A, 1, nc}, {B, 1, nc}];

  M1 = fold[XiL];   (* Q on the FIRST  F-vertex leg *)
  M2 = fold[XiM];   (* Q on the SECOND F-vertex leg *)
  M3 = fold[XiR];   (* Q on the kappa1 (external) leg *)
  M  = Expand[M1 + M2 + M3];

  truth = Expand[M + Transpose[M]];
  lhs   = Expand[(T1 + T2) + Transpose[T1 + T2]];

  (* ---- leg-by-leg placement bookkeeping (the "why") ---------------- *)
  pl1 = Expand[M3 - T1];
  pl2 = Expand[(M1 + M2) - Transpose[T2]];
  pl3 = Expand[M1 - Transpose[T2]];   (* must be NONzero: needs BOTH legs *)
  pl4 = Expand[M2 - Transpose[T2]];   (* must be NONzero *)
  Print["  placement: M3 (Q on kappa1 leg)      == T1    : ",
        If[Union[Flatten[{pl1}]] === {0}, "PASS", "FAIL"]];
  Print["  placement: M1+M2 (Q on F-vertex legs)== T2^T  : ",
        If[Union[Flatten[{pl2}]] === {0}, "PASS", "FAIL"]];
  Print["  placement: M1 alone != T2^T, M2 alone != T2^T : ",
        Union[Flatten[{pl3}]] =!= {0} && Union[Flatten[{pl4}]] =!= {0}];

  (* ---------------- non-triviality --------------------------------- *)
  nTerms = Total[Flatten[Map[Length[If[Head[#] === Plus, #, {#}]] &, truth, {2}]]];
  nzTruth = Union[Flatten[truth]] =!= {0};
  Print["  TRUTH nonzero: ", nzTruth, "   (total monomials in TRUTH: ", nTerms, ")"];

  (* ---------------- THE IDENTITY ----------------------------------- *)
  d0 = Simplify[Together[Expand[lhs - truth]]];
  Print[">>> MAIN IDENTITY  (T1+T2)+(T1+T2)^T - TRUTH  ==  ",
        If[Union[Flatten[{d0}]] === {0}, "0   ***PASS***",
           "NONZERO   ***FAIL***"]];

  (* sharper, unsymmetrised statement: M = T1 + T2^T exactly *)
  d1 = Simplify[Together[Expand[M - (T1 + Transpose[T2])]]];
  Print[">>> SHARP FORM     M - (T1 + T2^T)             ==  ",
        If[Union[Flatten[{d1}]] === {0}, "0   ***PASS***",
           "NONZERO   ***FAIL***"]];

  (* ---------------- teeth: the test must be able to fail ------------ *)
  dA = Expand[(T1 + Transpose[T1]) - truth];
  dB = Expand[(T2 + Transpose[T2]) - truth];
  dC = Expand[(T1 + T2) - truth];                 (* unsymmetrised LHS *)
  Print["  teeth: T1-only  != TRUTH : ", Union[Flatten[{dA}]] =!= {0}];
  Print["  teeth: T2-only  != TRUTH : ", Union[Flatten[{dB}]] =!= {0}];
  Print["  teeth: (T1+T2) unsymmetrised != TRUTH : ",
        Union[Flatten[{dC}]] =!= {0}];

  (* --- falsification control A: drop the F_{Bcb} half of kappa_d2 --- *)
  kd2bad = Table[
    Sum[cw[j] Sum[ff[B, b, c] s1[[j + 1, b]] sdv[[j + 1, c]],
        {b, 1, nc}, {c, 1, nc}], {j, 0, nd - 1}],
    {B, 1, nc}];
  T2bad = Table[Expand[EW[k1[[A]] kd2bad[[B]]]], {A, 1, nc}, {B, 1, nc}];
  lhsbad = Expand[(T1 + T2bad) + Transpose[T1 + T2bad]];
  ctlA = Union[Flatten[{Expand[lhsbad - truth]}]] =!= {0};
  Print["  control A: unsymmetrised F in kappa_d2 BREAKS identity : ", ctlA];

  (* --- falsification control B: mismatch kappa_d1 kernel vs kappa1 -- *)
  kd1bad = Table[Sum[cw[p] ww[p] dfv[[p + 1, B]], {p, 0, nd - 1}], {B, 1, nc}];
  T1bad = Table[Expand[EW[k2[[A]] kd1bad[[B]]]], {A, 1, nc}, {B, 1, nc}];
  ctlB = Union[Flatten[{Expand[
      ((T1bad + T2) + Transpose[T1bad + T2]) - truth]}]] =!= {0};
  Print["  control B: mismatched kappa_d1 kernel BREAKS identity   : ", ctlB];

  Print["  elapsed: ", Round[AbsoluteTime[] - t0, 0.01], " s"];
  Print[""];

  {Union[Flatten[{d0}]] === {0}, Union[Flatten[{d1}]] === {0},
   nzTruth, Union[Flatten[{dA}]] =!= {0}, Union[Flatten[{dB}]] =!= {0},
   Union[Flatten[{dC}]] =!= {0},
   Union[chkComp] === {0}, Union[Flatten[{chkKer}]] === {0},
   Union[Flatten[{pl1}]] === {0}, Union[Flatten[{pl2}]] === {0},
   Union[Flatten[{pl3}]] =!= {0}, Union[Flatten[{pl4}]] =!= {0},
   ctlA, ctlB}
];

(* ------------------------------------------------------------------ *)
(* 4.  Run                                                             *)
(* ------------------------------------------------------------------ *)

Print["############ SYMBOLIC PROOF OF THE FK Q-PLACEMENT IDENTITY ############"];
Print["Wolfram Language version: ", $Version];
Print[""];
Print["Engine self-tests:"];
selfOK = engineSelfTest[];
Print[""];

cases = {{2, 2}, {2, 3}, {2, 4}, {3, 3}, {3, 4}};
results = Table[runCase[cs[[1]], cs[[2]]], {cs, cases}];

Print["==================================================================="];
Print["SUMMARY"];
Print["==================================================================="];
Do[
 Print["  (n_comp,N) = ", cases[[i]],
   "   main=", If[results[[i, 1]], "PASS", "FAIL"],
   "  sharp=", If[results[[i, 2]], "PASS", "FAIL"],
   "  placement(M3=T1, M1+M2=T2^T)=",
     {If[results[[i, 9]], "PASS", "FAIL"], If[results[[i, 10]], "PASS", "FAIL"]},
   "  nontrivial=", results[[i, 3]],
   "  teeth=", {results[[i, 4]], results[[i, 5]], results[[i, 6]]},
   "  single-leg-insufficient=", {results[[i, 11]], results[[i, 12]]},
   "  controls=", {results[[i, 13]], results[[i, 14]]},
   "  checks=", {results[[i, 7]], results[[i, 8]]}],
 {i, Length[cases]}];

allOK = selfOK && And @@ Flatten[results];
Print[""];
Print[If[allOK,
  "OVERALL: ***PASS*** - the identity holds identically in generic symbols.",
  "OVERALL: ***FAIL***"]];
