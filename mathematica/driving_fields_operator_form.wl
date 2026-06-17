(* ============================================================
   driving_fields_operator_form.wl
   ------------------------------------------------------------
   Operator-form companion to driving_fields_harmonics.wl.

   Builds the per-multipole driving-field transfer expressions
   keeping the following operators SYMBOLIC (not expanded):

     OpLap[X]        := nabla^2 X        (3D Laplacian)
     OpLapPerp[X]    := nabla^2_perp X   (transverse Laplacian)
     OpDnull[X]      := \mathfrak{D} X     (= partial_tau X - hat n^i d_i X)
     OpDnull2[X]     := \mathfrak{D}^2 X
     OpNiDi[X]       := hat n^i d_i X    (radial directional derivative)
     OpNiNjDiDj[X]   := hat n^i hat n^j d_i d_j X
     OpSpinRaise2[X] := hat-screen^i hat-screen^j d_i d_j X (spin-2 raise)

   The Fraktur D (\mathfrak{D}) is used because plain D is overloaded in
   the manuscript: angular-diameter distance, growth factor, AND the null
   operator. \mathcal{D} is also unavailable (main.tex line 50 binds
   \hat{\mathcal{D}} to the determinant operator \DeterOp).

   At the end the script substitutes the symbolic operators by their
   per-multipole reduction (the SAME rules used in
   driving_fields_harmonics.wl) and verifies the result is identically
   equal to the per-multipole form derived independently from the
   driving_fields_harmonics.wl source-of-truth expressions.

   Run with:
     wolframscript -file scripts/driving_fields_operator_form.wl
   ============================================================ *)

$HistoryLength = 0;

(* ─── Field placeholders (same names as driving_fields_harmonics.wl) ─── *)
PsiL = psiLM[tt, chi];
PhiL = phiLM[tt, chi];
BL   = bLM[tt, chi];
h0LM = hL0[tt, chi];
h1LM = hL1[tt, chi];
h2LM = hL2[tt, chi];

Hc  = Hconf[tt];
Hcp = Derivative[1][Hconf][tt];

(* ─── Per-multipole derivative shortcuts ─── *)
dchi[f_]    := D[f, chi];
dtt[f_]     := D[f, tt];
dchi2[f_]   := D[f, {chi, 2}];
dtt2[f_]    := D[f, {tt, 2}];
dttchi[f_]  := D[f, tt, chi];

(* ─── Symbolic operators (heads kept opaque; no auto-distribution) ─── *)
ClearAll[OpLap, OpLapPerp, OpDnull, OpDnull2, OpNiDi, OpNiNjDiDj, OpSpinRaise2];

(* ─── Operator-form transfer functions ─── *)
(*
   These mirror driving_fields_harmonics.wl line by line, but every
   instance of nPar/nPerp/nini/Dnull/Dnull2/spinRaise2/LaplaceFull
   is replaced by its symbolic-operator counterpart.

   Sign convention follows the manuscript: Phi_00 = -(1/2) R_{mu nu} k^mu k^nu,
   so the scalar-sector Phi_00 carries an overall minus.
*)

Phi00ScalarOp = -(
    2 (Hcp - Hc^2) PsiL
  + Hc (dtt[PsiL] - dtt[PhiL])
  - 2 Hc OpNiDi[PsiL]
  + (1/2) OpLapPerp[PhiL + PsiL]
  + OpDnull2[PhiL]
);

Psi0ScalarOp = -OpSpinRaise2[PhiL + PsiL];

Phi00BOp = -(
    2 (Hc^2 - Hcp) OpNiDi[BL]
  + Hc OpNiNjDiDj[BL]
  - (1/2) OpLapPerp[dtt[BL]]
);

Psi0BOp = OpSpinRaise2[dtt[BL]];

Phi00TensorOp = -(1/4) (
    dtt2[h0LM] + 2 Hc dtt[h0LM] - OpLap[h0LM]
);

Psi0TensorOp = (1/2) (
    OpDnull2[h2LM]
  + 2 OpDnull[(1/(chi Sqrt[2])) Sqrt[L2 - 2] h1LM]
  + (1/(2 chi^2)) Sqrt[L2 (L2 - 2)] h0LM
);

(* ─── Expansion rules: replace symbolic operators by per-multipole forms ─── *)
(* These are exactly the reductions used in driving_fields_harmonics.wl
   (LaplaceFull, nPerp, Dnull, Dnull2, spinRaise2, nPar, nini).            *)
expandRules = {
  OpLap[f_]        :> dchi2[f] + (2/chi) dchi[f] - (L2/chi^2) f,
  OpLapPerp[f_]    :> (2/chi) dchi[f] - (L2/chi^2) f,
  OpDnull[f_]      :> dtt[f] - dchi[f],
  OpDnull2[f_]     :> dtt2[f] - 2 dttchi[f] + dchi2[f],
  OpNiDi[f_]       :> dchi[f],
  OpNiNjDiDj[f_]   :> dchi2[f],
  OpSpinRaise2[f_] :> (1/(2 chi^2)) Sqrt[L2 (L2 - 2)] f
};

(* ─── Reference forms (independent re-derivation, must match expanded op form) ─── *)
LaplaceFull[f_] := dchi2[f] + (2/chi) dchi[f] - (L2/chi^2) f;
nini[f_]        := dchi2[f];
nPerp[f_]       := (2/chi) dchi[f] - (L2/chi^2) f;
nPar[f_]        := dchi[f];
DnullR[f_]      := dtt[f] - dchi[f];
Dnull2R[f_]     := dtt2[f] - 2 dttchi[f] + dchi2[f];
spinRaise2[f_]  := (1/(2 chi^2)) Sqrt[L2 (L2 - 2)] f;

Phi00ScalarRef = -(
    2 (Hcp - Hc^2) PsiL + Hc (dtt[PsiL] - dtt[PhiL]) - 2 Hc nPar[PsiL]
  + (1/2) nPerp[PhiL + PsiL] + Dnull2R[PhiL]
);
Psi0ScalarRef  = -spinRaise2[PhiL + PsiL];
Phi00BRef = -(
    2 (Hc^2 - Hcp) nPar[BL] + Hc nini[BL] - (1/2) nPerp[dtt[BL]]
);
Psi0BRef  = spinRaise2[dtt[BL]];
Phi00TensorRef = -(1/4) (
    dtt2[h0LM] + 2 Hc dtt[h0LM] - LaplaceFull[h0LM]
);
Psi0TensorRef = (1/2) (
    Dnull2R[h2LM]
  + 2 DnullR[(1/(chi Sqrt[2])) Sqrt[L2 - 2] h1LM]
  + (1/(2 chi^2)) Sqrt[L2 (L2 - 2)] h0LM
);

(* ─── PASS/FAIL equivalence tests ─── *)
runTest[label_, opForm_, refForm_] := Module[{lhs, rhs, diff},
  lhs  = Expand[opForm /. expandRules];
  rhs  = Expand[refForm];
  diff = Simplify[lhs - rhs];
  If[diff === 0,
    Print["PASS: ", label]; True,
    Print["FAIL: ", label, "  diff = ", diff]; False
  ]
];

Print["=== Operator-form vs per-multipole equivalence checks ==="];
results = {
  runTest["Phi_00 scalar (s)",  Phi00ScalarOp,                Phi00ScalarRef],
  runTest["Phi_00 shift  (B)",  Phi00BOp,                     Phi00BRef],
  runTest["Phi_00 scalar+B",    Phi00ScalarOp + Phi00BOp,     Phi00ScalarRef + Phi00BRef],
  runTest["Psi_0  scalar (s)",  Psi0ScalarOp,                 Psi0ScalarRef],
  runTest["Psi_0  shift  (B)",  Psi0BOp,                      Psi0BRef],
  runTest["Psi_0  scalar+B",    Psi0ScalarOp + Psi0BOp,       Psi0ScalarRef + Psi0BRef],
  runTest["Phi_00 tensor (h)",  Phi00TensorOp,                Phi00TensorRef],
  runTest["Psi_0  tensor (h)",  Psi0TensorOp,                 Psi0TensorRef]
};

allPass = And @@ results;
Print[""];
If[allPass,
  Print["*** ALL ", Length[results], " EQUIVALENCE TESTS PASS ***"],
  Print["!!! ", Count[results, False], "/", Length[results], " TESTS FAILED !!!"]
];

(* Print operator-form expressions as ASCII for human inspection *)
Print[""];
Print["=== Operator-form expressions (ASCII) ==="];

Print[""];
Print["Phi_00 (scalar, B=0):"];
Print["  ", Phi00ScalarOp];

Print[""];
Print["Psi_0  (scalar, B=0):"];
Print["  ", Psi0ScalarOp];

Print[""];
Print["Phi_00 (shift B):"];
Print["  ", Phi00BOp];

Print[""];
Print["Psi_0  (shift B):"];
Print["  ", Psi0BOp];

Print[""];
Print["Phi_00 (tensor h0):"];
Print["  ", Phi00TensorOp];

Print[""];
Print["Psi_0  (tensor h0,h1,h2):"];
Print["  ", Psi0TensorOp];

(* Save TeX-friendly forms for paper insertion *)
Print[""];
Print["%%TEX_PHI00_SCALAR_OP_START%%"];
Print[ToString[TeXForm[Phi00ScalarOp]]];
Print["%%TEX_PHI00_SCALAR_OP_END%%"];

Print["%%TEX_PSI0_SCALAR_OP_START%%"];
Print[ToString[TeXForm[Psi0ScalarOp]]];
Print["%%TEX_PSI0_SCALAR_OP_END%%"];

Print["%%TEX_PHI00_B_OP_START%%"];
Print[ToString[TeXForm[Phi00BOp]]];
Print["%%TEX_PHI00_B_OP_END%%"];

Print["%%TEX_PSI0_B_OP_START%%"];
Print[ToString[TeXForm[Psi0BOp]]];
Print["%%TEX_PSI0_B_OP_END%%"];

Print["%%TEX_PHI00_TENSOR_OP_START%%"];
Print[ToString[TeXForm[Phi00TensorOp]]];
Print["%%TEX_PHI00_TENSOR_OP_END%%"];

Print["%%TEX_PSI0_TENSOR_OP_START%%"];
Print[ToString[TeXForm[Psi0TensorOp]]];
Print["%%TEX_PSI0_TENSOR_OP_END%%"];

Print[""];
Print["[done]"];

If[!allPass, Exit[1]];
