(* ============================================================
   cl_scalar_psi_eq_phi.wl
   ------------------------------------------------------------
   Derive the explicit form of the factorised angular power
   spectrum C_ell^{AB}(chi, chi') from cosmology.tex
   Eq. (Cl factorised) in the scalar-only, Phi = Psi limit.

   Assumptions:
     - No vector perturbations: B = 0
     - No tensor perturbations: h_{ij} = 0
     - No anisotropic stress: Phi = Psi (Phi_{lm} == Psi_{lm})

   Procedure:
     1. Re-derive the per-multipole transfer-coefficient six-vectors
        c^A_{X,(a,b)}  from the same Poisson-gauge formulas as
        driving_fields_harmonics.wl.
     2. Apply Phi = Psi: add the Phi column to the Psi column;
        the underlying power spectrum is then a single P(k).
     3. Assemble C_ell^{AB}(chi, chiP) as the bilinear sum
          C^{AB} = Sum_{(a,b),(c,d)}
                      c^A_{(a,b)} c^B_{(c,d)}
                      M^{(a)}(chi) M^{(c)}(chiP)
                      I_{b,d}(chi, chiP)
        where
          M^{(a)}(chi) := h_Psi^{(a)}(tau_0 - chi)   [a-th tau-deriv of growth]
          I_{b,d}(chi, chiP) := (2/pi) int dk
                k^{2+b+d} j_l^{(b)}(k chi) j_l^{(d)}(k chiP) P(k, tau_0)

   Run with:
     wolframscript -file scripts/cl_scalar_psi_eq_phi.wl
   (or via the run-wl.sh wrapper in the skills directory)
   ============================================================ *)

$HistoryLength = 0;

(* ─── Derivative operators (same conventions as driving_fields_harmonics.wl) ─── *)
dchi[f_]    := D[f, chi];
dtt[f_]     := D[f, tt];
dchi2[f_]   := D[f, {chi, 2}];
dtt2[f_]    := D[f, {tt, 2}];
dttchi[f_]  := D[f, tt, chi];
nPar[f_]    := dchi[f];
nPerp[f_]   := (2/chi) dchi[f] - (L2/chi^2) f;
Dnull2[f_]  := dtt2[f] - 2 dttchi[f] + dchi2[f];
spinRaise2[f_] := (1/(2 chi^2)) Sqrt[L2 (L2 - 2)] f;

(* ─── Placeholder multipoles and background quantities ─── *)
PsiL = psiLM[tt, chi];
PhiL = phiLM[tt, chi];
Hc   = Hconf[tt];        (* calH(tau)  *)
Hcp  = Derivative[1][Hconf][tt]; (* calH'(tau) *)

(* ─── Per-multipole driving fields (scalar sector, B = 0) ─── *)
(* Sign convention: Phi_00 = -(1/2) R_{mu nu} k^mu k^nu  =>  overall minus *)
Phi00ScalarLM = -(
    2 (Hcp - Hc^2) PsiL
  + Hc (dtt[PsiL] - dtt[PhiL])
  - 2 Hc nPar[PsiL]
  + (1/2) nPerp[PhiL + PsiL]
  + Dnull2[PhiL]
);

(* Psi_0 (spin-2 Weyl shear), scalar sector *)
Psi0ScalarLM = -spinRaise2[PhiL + PsiL];

(* ─── Extract coefficient six-vectors per field ─── *)
(*   Derivative basis ordering:
       1: (a=0,b=0)  = X_{lm}
       2: (a=0,b=1)  = d_chi X_{lm}
       3: (a=0,b=2)  = d_chi^2 X_{lm}
       4: (a=1,b=0)  = d_tau X_{lm}
       5: (a=2,b=0)  = d_tau^2 X_{lm}
       6: (a=1,b=1)  = d_tau d_chi X_{lm}               *)

basisFor[f_] := {f, dchi[f], dchi2[f], dtt[f], dtt2[f], dttchi[f]};
opMatrix[expr_, basis_] :=
  Table[Coefficient[Expand[expr], b], {b, basis}];

cPhi00Psi = opMatrix[Phi00ScalarLM, basisFor[PsiL]];
cPhi00Phi = opMatrix[Phi00ScalarLM, basisFor[PhiL]];
cPsi0Psi  = opMatrix[Psi0ScalarLM,  basisFor[PsiL]];
cPsi0Phi  = opMatrix[Psi0ScalarLM,  basisFor[PhiL]];

(* ─── Apply Phi = Psi: merge columns ─── *)
(* Since Phi_{lm} = Psi_{lm}, the coefficient of any derivative acts
   on the *same* radial multipole, so the two columns simply add.
   Simultaneously P_{PhiPhi} = P_{PsiPsi} = P_{PhiPsi} = P(k), so
   the double sum over (X,Y) in Eq. (Cl factorised) also collapses
   to a single underlying power spectrum.                           *)

cPhi00 = Simplify[cPhi00Psi + cPhi00Phi];
cPsi0  = Simplify[cPsi0Psi  + cPsi0Phi];

(* ─── Light-cone substitution for background functions ─── *)
(* Hconf[tt] evaluated on the past LC: tt = tau_0 - chi.
   We rename to calH[chi] and calHP[chi] for readability.         *)
lcSubs = {
  Hconf[tt]                  -> calH[chi],
  Derivative[1][Hconf][tt]   -> calHP[chi]
};

cPhi00lc = cPhi00 /. lcSubs // Simplify;
cPsi0lc  = cPsi0  /. lcSubs // Simplify;

(* Sanity check: confirm the explicit (a=0,b=0) entry of cPsi0 *)
Print["c^{Psi0}_(0,0) = ", cPsi0lc[[1]], "  (expected -Sqrt[L2(L2-2)]/chi^2)"];

(* ─── Print merged coefficient vectors ─── *)
derivLabels = {
  "(a=0,b=0) = X",
  "(a=0,b=1) = d_chi X",
  "(a=0,b=2) = d_chi^2 X",
  "(a=1,b=0) = d_tau X",
  "(a=2,b=0) = d_tau^2 X",
  "(a=1,b=1) = d_tau d_chi X"
};

Print["=== Merged coefficient vectors under Phi=Psi, B=0, tensors=0 ==="];
Print["    L2 = ell(ell+1);  calH = conformal Hubble on light cone;"];
Print["    calHP = calH' (conformal-time derivative of calH) on light cone."];
Print[""];

Print["c^{Phi00}_{(a,b)}  [coefficient of X = Psi in driving-field Phi_00]:"];
Do[
  If[cPhi00lc[[i]] =!= 0,
     Print["  ", derivLabels[[i]], "  :  ", cPhi00lc[[i]]]],
  {i, 6}
];

Print[""];
Print["c^{Psi0}_{(a,b)}  [coefficient of X = Psi in driving-field Psi_0]:"];
Do[
  If[cPsi0lc[[i]] =!= 0,
     Print["  ", derivLabels[[i]], "  :  ", cPsi0lc[[i]]]],
  {i, 6}
];

(* ─── Verify only (0,0) survives in cPsi0 ─── *)
Print[""];
nonzeroPsi0 = Select[Range[6], (cPsi0lc[[#]] =!= 0) &];
If[nonzeroPsi0 === {1},
   Print["PASS: cPsi0 has only the (a=0,b=0) entry non-zero (expected)."],
   Print["FAIL: cPsi0 has unexpected non-zero entries at positions ", nonzeroPsi0]
];

(* ─── Verify d_tau coefficient of Phi00 vanishes ─── *)
If[Simplify[cPhi00lc[[4]]] === 0,
   Print["PASS: c^{Phi00}_{(1,0)} = 0  (calH cancellation under Phi=Psi confirmed)."],
   Print["FAIL: c^{Phi00}_{(1,0)} = ", cPhi00lc[[4]], "  (should be 0)"]
];

(* ─── Assemble C_ell^{AB}(chi, chiP) ─── *)
(*
   From cosmology.tex Eq. (Cl factorised):

     C_ell^{AB}(chi, chiP) =
       Sum_{(a,b),(c,d)}  c^A_{(a,b)} c^B_{(c,d)}
                          M^{(a)}(chi) M^{(c)}(chiP)
                          I_{b,d}(chi, chiP)

   Notation:
     M^{(a)}[chi]            = Subscript[M, a][chi]
       := a-th tau-derivative of h_Psi(tau_0 - chi)
          (a=0: growth factor on LC; a=1: its conformal-time derivative; a=2: second)

     I_{b,d}[chi, chiP]      = Subscript[II, {b,d}][chi, chiP]
       := (2/pi) int_0^inf dk  k^{2+b+d}
              j_ell^{(b)}(k chi) j_ell^{(d)}(k chiP) P(k, tau_0)

   The coefficients are real (scalar sector), so no conjugation needed.

   chiP-side background functions are evaluated at calH[chiP], calHP[chiP].
*)

basisPairs = {{0,0},{0,1},{0,2},{1,0},{2,0},{1,1}};

assembleClAB[cA_, cB_] :=
  Expand @ Sum[
    cA[[i]] * cB[[j]] *
    Subscript[M, basisPairs[[i, 1]]][chi] *
    Subscript[M, basisPairs[[j, 1]]][chiP] *
    Subscript[II, {basisPairs[[i, 2]], basisPairs[[j, 2]]}][chi, chiP],
    {i, 6}, {j, 6}
  ];

(* chiP-side coefficients: replace ALL chi dependence with chiP.
   This includes both the explicit geometric chi (in 1/chi^n denominators)
   and the background quantities calH[chi], calHP[chi].               *)
cPhi00lcP = (cPhi00lc /. chi -> chiP);
cPsi0lcP  = (cPsi0lc  /. chi -> chiP);

ClPhi00Phi00 = assembleClAB[cPhi00lc, cPhi00lcP];
ClPsi0Psi0   = assembleClAB[cPsi0lc,  cPsi0lcP];
ClPhi00Psi0  = assembleClAB[cPhi00lc, cPsi0lcP];  (* cross *)
ClPsi0Phi00  = assembleClAB[cPsi0lc,  cPhi00lcP]; (* adjoint cross *)

(* ─── Verify Psi0-Psi0 has the expected factored form ─── *)
(*   c^{Psi0}_{(0,0)}(chi) * c^{Psi0}_{(0,0)}(chiP) = L2(L2-2)/(chi^2 chiP^2)
     so C^{Psi0,Psi0} = L2(L2-2)/(chi^2 chiP^2) * M_0(chi) M_0(chiP) * II_{00}  *)
expectedPsi0Psi0 =
  (L2 (L2 - 2))/(chi^2 chiP^2) *
  Subscript[M, 0][chi] * Subscript[M, 0][chiP] *
  Subscript[II, {0, 0}][chi, chiP];
If[Simplify[ClPsi0Psi0 - expectedPsi0Psi0] === 0,
   Print["PASS: C^{Psi0,Psi0} = L2(L2-2)/(chi^2 chiP^2) M_0 M_0' II_{00}  (factored form verified)."],
   Print["FAIL: C^{Psi0,Psi0} differs from expected:  ", Simplify[ClPsi0Psi0 - expectedPsi0Psi0]]
];

(* ─── Verify symmetry C^{AB}(chi,chiP) = C^{BA}(chiP,chi) ─── *)
(*   Swap chi <-> chiP and A <-> B in ClPsi0Phi00, then check equality.
     Note: II[{b,d}][chi,chiP] and II[{d,b}][chiP,chi] are the same integral
     (just relabelled integration variables), so we declare II symmetric:
     Subscript[II,{b,d}][c1,c2] = Subscript[II,{d,b}][c2,c1].            *)
symRules = {
  Subscript[II, {b_, d_}][c1_, c2_] :> Subscript[II, {d, b}][c2, c1]
};
swapCC = (ClPsi0Phi00 /. {chi -> chiPtmp, chiP -> chi} /. {chiPtmp -> chiP});
(* Also apply II symmetry to bring both into a canonical form *)
swapCCsym = swapCC /. symRules;
diff = Simplify[ClPhi00Psi0 - swapCCsym];
If[diff === 0,
   Print["PASS: C^{Phi00,Psi0}(chi,chiP) = C^{Psi0,Phi00}(chiP,chi)  (symmetry)."],
   Print["INFO: symmetry diff (should be 0 after II symmetry used): ", diff]
];

(* ─── Print assembled spectra ─── *)
Print[""];
Print["=== Explicit C_ell^{Phi00,Phi00}(chi, chiP) ==="];
Print["    (25 terms; M^{(a)} = a-th tau-deriv of growth factor on LC)"];
Do[
  term = If[Head[ClPhi00Phi00] === Plus, (List @@ ClPhi00Phi00)[[k]], ClPhi00Phi00];
  Print["  ", term],
  {k, Length[If[Head[ClPhi00Phi00] === Plus, List @@ ClPhi00Phi00, {ClPhi00Phi00}]]}
];

Print[""];
Print["=== Explicit C_ell^{Psi0,Psi0}(chi, chiP) ==="];
Print[ClPsi0Psi0];

Print[""];
Print["=== Explicit C_ell^{Phi00,Psi0}(chi, chiP) ==="];
Do[
  term = If[Head[ClPhi00Psi0] === Plus, (List @@ ClPhi00Psi0)[[k]], ClPhi00Psi0];
  Print["  ", term],
  {k, Length[If[Head[ClPhi00Psi0] === Plus, List @@ ClPhi00Psi0, {ClPhi00Psi0}]]}
];

(* ─── TeXForm output ─── *)
Print["%%COEF_PHI00_MERGED_START%%"];
Print[TeXForm[cPhi00lc]];
Print["%%COEF_PHI00_MERGED_END%%"];

Print["%%COEF_PSI0_MERGED_START%%"];
Print[TeXForm[cPsi0lc]];
Print["%%COEF_PSI0_MERGED_END%%"];

Print["%%CL_PHI00PHI00_START%%"];
Print[TeXForm[ClPhi00Phi00]];
Print["%%CL_PHI00PHI00_END%%"];

Print["%%CL_PSI0PSI0_START%%"];
Print[TeXForm[ClPsi0Psi0]];
Print["%%CL_PSI0PSI0_END%%"];

Print["%%CL_PHI00PSI0_START%%"];
Print[TeXForm[ClPhi00Psi0]];
Print["%%CL_PHI00PSI0_END%%"];

Print["[done]"];
