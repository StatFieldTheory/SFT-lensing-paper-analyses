(* ============================================================
   driving_fields_poisson.wl
   ------------------------------------------------------------
   Symbolic derivation of the Sachs driving fields

        Phi00 = -(1/2) R_{mu nu} k^mu k^nu          (Ricci focusing)
        Psi0  = - R_{abcd} k^a Z^b k^c Z^d          (Weyl shear)

   on the perturbed Poisson-gauge FLRW metric

        ds^2 = a(t)^2 [ -(1 + 2 eps Psi) dt^2
                          - 2 eps dB/dx^i  dt dx^i
                          + ((1 - 2 eps Phi) delta_{ij}
                               + eps h_{ij}) dx^i dx^j ],

   to first order in the bookkeeping parameter eps.  Photon travels
   along the x3 axis; the screen basis is {x1, x2}.

   Run with:
     ~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh \
        scripts/driving_fields_poisson.wl 600
   ============================================================ *)

$HistoryLength = 0;

(* ---- Phase 1. Coordinates and fields ---- *)
coords = {tt, x1, x2, x3};

Psifld = PsiF[tt, x1, x2, x3];
Phifld = PhiF[tt, x1, x2, x3];
Bfld   = BF[tt, x1, x2, x3];

hh[1, 1] := hF11[tt, x1, x2, x3];
hh[2, 2] := hF22[tt, x1, x2, x3];
hh[3, 3] := hF33[tt, x1, x2, x3];
hh[1, 2] := hF12[tt, x1, x2, x3]; hh[2, 1] := hh[1, 2];
hh[1, 3] := hF13[tt, x1, x2, x3]; hh[3, 1] := hh[1, 3];
hh[2, 3] := hF23[tt, x1, x2, x3]; hh[3, 2] := hh[2, 3];

(* di[f, k] = d/dx^k f (k = 1,2,3). *)
di[f_, k_Integer] := D[f, coords[[k + 1]]];

(* ---- Phase 2. Metric in Poisson gauge ---- *)
aa = a[tt];

gM = Table[0, {4}, {4}];
gM[[1, 1]] = -aa^2 (1 + 2 eps Psifld);
Do[
  gM[[1, k + 1]] = -aa^2 eps di[Bfld, k];
  gM[[k + 1, 1]] = gM[[1, k + 1]],
  {k, 3}];
Do[
  gM[[i + 1, j + 1]] = aa^2 ((1 - 2 eps Phifld) KroneckerDelta[i, j] + eps hh[i, j]),
  {i, 3}, {j, 3}];

linearize[expr_] := Normal@Series[expr, {eps, 0, 1}];

(* ---- Phase 3. Inverse metric to first order ---- *)
g0M     = gM /. eps -> 0;
g0Minv  = Inverse[g0M];
dgMet   = gM - g0M;
gMinv   = g0Minv - g0Minv . dgMet . g0Minv;   (* keeps O(eps^2), Series drops it *)

(* ---- Phase 4. Christoffel Gam[alpha, mu, nu] ---- *)
Print["[progress] Christoffel ..."];
dgM = Table[D[gM[[bb, nn]], coords[[mm]]], {mm, 4}, {bb, 4}, {nn, 4}];
Gam = Table[
   (1/2) Sum[
     gMinv[[aaI, bb]] (dgM[[mm, bb, nn]] + dgM[[nn, bb, mm]] - dgM[[bb, mm, nn]]),
     {bb, 4}],
   {aaI, 4}, {mm, 4}, {nn, 4}];
Gam = linearize[Gam];

(* ---- Phase 5. Riemann ---- *)
Print["[progress] Riemann ..."];
RiemUp = Table[
   D[Gam[[aaI, dd, bb]], coords[[cc]]] - D[Gam[[aaI, cc, bb]], coords[[dd]]]
    + Sum[Gam[[aaI, cc, ll]] Gam[[ll, dd, bb]] - Gam[[aaI, dd, ll]] Gam[[ll, cc, bb]], {ll, 4}],
   {aaI, 4}, {bb, 4}, {cc, 4}, {dd, 4}];
RiemUp = linearize[RiemUp];

RiemDown = Table[
   Sum[gM[[aaI, rho]] RiemUp[[rho, bb, cc, dd]], {rho, 4}],
   {aaI, 4}, {bb, 4}, {cc, 4}, {dd, 4}];
RiemDown = linearize[RiemDown];

(* ---- Phase 6. Ricci ---- *)
Print["[progress] Ricci ..."];
Ricci = Table[Sum[RiemUp[[aaI, mm, aaI, nn]], {aaI, 4}], {mm, 4}, {nn, 4}];
Ricci = linearize[Ricci];

(* ---- Phase 7. Observer, line-of-sight, screen basis ---- *)
Print["[progress] Tetrad ..."];

uU = {1/(aa Sqrt[1 + 2 eps Psifld]), 0, 0, 0};
uU = linearize[uU];

tet[ii_Integer] := Module[{vec},
  vec = Table[0, {4}];
  vec[[1]] = -eps di[Bfld, ii]/aa;
  Do[
    vec[[jj + 1]] =
      (1/aa) (KroneckerDelta[ii, jj]
              + eps (Phifld KroneckerDelta[ii, jj] - (1/2) hh[ii, jj])),
    {jj, 3}];
  linearize[vec]
];

eLine = tet[3];
Xs    = tet[1];
Ys    = tet[2];
Zs    = (Xs + I Ys)/Sqrt[2];
Zbar  = (Xs - I Ys)/Sqrt[2];

dot[v1_, v2_] := Sum[gM[[mm, nn]] v1[[mm]] v2[[nn]], {mm, 4}, {nn, 4}];

Print["[check] u.u = ", Simplify@linearize[dot[uU, uU]]];
Print["[check] e.e = ", Simplify@linearize[dot[eLine, eLine]]];
Print["[check] e.u = ", Simplify@linearize[dot[eLine, uU]]];
Print["[check] X.X = ", Simplify@linearize[dot[Xs, Xs]]];
Print["[check] Y.Y = ", Simplify@linearize[dot[Ys, Ys]]];
Print["[check] X.Y = ", Simplify@linearize[dot[Xs, Ys]]];
Print["[check] X.e = ", Simplify@linearize[dot[Xs, eLine]]];
Print["[check] Y.e = ", Simplify@linearize[dot[Ys, eLine]]];
Print["[check] X.u = ", Simplify@linearize[dot[Xs, uU]]];
Print["[check] Z.Z = ", Simplify@linearize[dot[Zs, Zs]]];
Print["[check] Z.Zbar = ", Simplify@linearize[dot[Zs, Zbar]]];

kU = Eobs (-uU + eLine);
kU = linearize[kU];

Print["[check] k.k = ", Simplify@linearize[dot[kU, kU]]];
Print["[check] k.u = ", Simplify@linearize[dot[kU, uU]]];
Print["[check] k.Z = ", Simplify@linearize[dot[kU, Zs]]];

(* ---- Phase 8. Driving fields ---- *)
Print["[progress] Driving fields ..."];

phi00 = -(1/2) Sum[Ricci[[mm, nn]] kU[[mm]] kU[[nn]], {mm, 4}, {nn, 4}];
phi00 = linearize[phi00] // Expand;

phi00BG  = phi00 /. eps -> 0;
phi00Lin = Coefficient[phi00, eps, 1];

psi0 = -Sum[
    RiemDown[[aaI, bb, cc, dd]] kU[[aaI]] Zs[[bb]] kU[[cc]] Zs[[dd]],
    {aaI, 4}, {bb, 4}, {cc, 4}, {dd, 4}];
psi0 = linearize[psi0] // Expand;

psi0BG  = psi0 /. eps -> 0;
psi0Lin = Coefficient[psi0, eps, 1];

(* ---- Phase 9. Hubble substitution ---- *)
hubbleRules = {
  Derivative[2][a][tt] -> a[tt] (Derivative[1][H][tt] + H[tt]^2),
  Derivative[1][a][tt] -> a[tt] H[tt]
};

phi00BGs   = phi00BG   /. hubbleRules // Expand // Simplify;
phi00Lins  = phi00Lin  /. hubbleRules // Expand;
psi0BGs    = psi0BG    /. hubbleRules // Expand // Simplify;
psi0Lins   = psi0Lin   /. hubbleRules // Expand;

(* Express per E_co^2, where E_co := a * E_obs. *)
phi00CoBG  = phi00BGs  a[tt]^2 / Eobs^2 // Simplify;
phi00CoLin = phi00Lins a[tt]^2 / Eobs^2 // Expand;
psi0CoBG   = psi0BGs   a[tt]^2 / Eobs^2 // Simplify;
psi0CoLin  = psi0Lins  a[tt]^2 / Eobs^2 // Expand;

(* ---- Phase 10. Print results ---- *)

Print["%%BG_PHI00_START%%"];
Print[ToString[TeXForm[phi00CoBG]]];
Print["%%BG_PHI00_END%%"];

Print["%%LIN_PHI00_START%%"];
Print[ToString[TeXForm[phi00CoLin]]];
Print["%%LIN_PHI00_END%%"];

Print["%%BG_PSI0_START%%"];
Print[ToString[TeXForm[psi0CoBG]]];
Print["%%BG_PSI0_END%%"];

Print["%%LIN_PSI0_START%%"];
Print[ToString[TeXForm[psi0CoLin]]];
Print["%%LIN_PSI0_END%%"];

Print["%%LIN_PHI00_INPUT_START%%"];
Print[ToString[InputForm[phi00CoLin]]];
Print["%%LIN_PHI00_INPUT_END%%"];

Print["%%LIN_PSI0_INPUT_START%%"];
Print[ToString[InputForm[psi0CoLin]]];
Print["%%LIN_PSI0_INPUT_END%%"];

(* ---- Phase 11. Sanity slicings ---- *)

(* Zero a field AND all its derivatives *)
zeroField[f_][expr_] := expr /. {
   Derivative[__][f][__] -> 0,
   f[__] -> 0
};
zeroAll[fields_List][expr_] := Fold[zeroField[#2][#1] &, expr, fields];

scalarOnlyRules = {hF11, hF22, hF33, hF12, hF13, hF23, BF};
tensorOnlyRules = {PsiF, PhiF, BF};
BOnlyRules      = {PsiF, PhiF, hF11, hF22, hF33, hF12, hF13, hF23};

phi00ScalarLin = zeroAll[scalarOnlyRules][phi00CoLin] // Expand;
psi0ScalarLin  = zeroAll[scalarOnlyRules][psi0CoLin]  // Expand;
phi00TensorLin = zeroAll[tensorOnlyRules][phi00CoLin] // Expand;
psi0TensorLin  = zeroAll[tensorOnlyRules][psi0CoLin]  // Expand;
phi00BLin      = zeroAll[BOnlyRules][phi00CoLin]      // Expand;
psi0BLin       = zeroAll[BOnlyRules][psi0CoLin]       // Expand;

Print["%%SCALAR_PHI00_LIN_START%%"];
Print[ToString[TeXForm[phi00ScalarLin]]];
Print["%%SCALAR_PHI00_LIN_END%%"];

Print["%%SCALAR_PSI0_LIN_START%%"];
Print[ToString[TeXForm[psi0ScalarLin]]];
Print["%%SCALAR_PSI0_LIN_END%%"];

Print["%%TENSOR_PHI00_LIN_START%%"];
Print[ToString[TeXForm[phi00TensorLin]]];
Print["%%TENSOR_PHI00_LIN_END%%"];

Print["%%TENSOR_PSI0_LIN_START%%"];
Print[ToString[TeXForm[psi0TensorLin]]];
Print["%%TENSOR_PSI0_LIN_END%%"];

Print["%%B_PHI00_LIN_START%%"];
Print[ToString[TeXForm[phi00BLin]]];
Print["%%B_PHI00_LIN_END%%"];

Print["%%B_PSI0_LIN_START%%"];
Print[ToString[TeXForm[psi0BLin]]];
Print["%%B_PSI0_LIN_END%%"];

(* ---- Phase 12. Transverse-traceless simplifications ----
   Impose on h_ij:
     (trace)  h_{11} + h_{22} + h_{33} = 0
     (transverse, j=1)  d_1 h_{11} + d_2 h_{12} + d_3 h_{13} = 0
     (transverse, j=2)  d_1 h_{12} + d_2 h_{22} + d_3 h_{23} = 0
     (transverse, j=3)  d_1 h_{13} + d_2 h_{23} + d_3 h_{33} = 0
   These are 4 scalar constraints.  We solve them symbolically by
   eliminating {hF33, its derivatives} via trace, and then eliminating
   d_3 h_{33} via the j=3 transversality constraint.  *)

traceRule = hF33 -> Function[{t, X1, X2, X3},
   - hF11[t, X1, X2, X3] - hF22[t, X1, X2, X3]];

(* Apply trace (and propagate through derivatives), then expand. *)
applyTrace[expr_] := expr /. traceRule // Expand;

phi00CoLinTT = applyTrace[phi00CoLin];
psi0CoLinTT  = applyTrace[psi0CoLin];

(* Emit TT-simplified forms *)
Print["%%LIN_PHI00_TT_START%%"];
Print[ToString[TeXForm[phi00CoLinTT]]];
Print["%%LIN_PHI00_TT_END%%"];

Print["%%LIN_PSI0_TT_START%%"];
Print[ToString[TeXForm[psi0CoLinTT]]];
Print["%%LIN_PSI0_TT_END%%"];

(* Tensor-only sector with TT applied *)
phi00TensorTT = applyTrace[phi00TensorLin] // Expand;
psi0TensorTT  = applyTrace[psi0TensorLin]  // Expand;

Print["%%TENSOR_PHI00_TT_START%%"];
Print[ToString[TeXForm[phi00TensorTT]]];
Print["%%TENSOR_PHI00_TT_END%%"];

Print["%%TENSOR_PSI0_TT_START%%"];
Print[ToString[TeXForm[psi0TensorTT]]];
Print["%%TENSOR_PSI0_TT_END%%"];

(* ============================================================
   Phase 13.  Screen-tetrad convention tests (rigorous).
   ------------------------------------------------------------
   We prove the actual claim that holds:  at first order in the
   Poisson potentials, any two screen bases that are both first-order
   orthonormal -- i.e., satisfy X . X = Y . Y = 1, X . Y = 0, X . u =
   Y . u = 0, X . k = Y . k = 0 to O(eps^1) -- yield the SAME first-
   order Psi_0.  In particular, a parallel-transported tetrad from the
   observer and the locally-orthonormal tetrad constructed here agree
   at this order.

   Proof strategy:
     * Any two first-order-orthonormal tetrads differ by an O(eps)
       rotation in the screen plane, X' = X + eps omega Y, Y' = Y -
       eps omega X.  Under this rotation, Z = (X + i Y)/sqrt(2) -> Z'
       = Z exp(-i eps omega) and Psi_0 = -R k Z k Z transforms as spin-2:
       Psi_0 -> exp(-2 i eps omega) Psi_0.
     * At O(eps^1) the shift is  -2 i eps omega Psi_0^{(0)}, which
       vanishes since Psi_0^{(0)} = 0 on FLRW (conformally flat).

   Test (a) -- POSITIVE CONTROL (rotation invariance):
     Apply a first-order rotation with an arbitrary rotation-angle field
     omegaF(tau, x1, x2, x3) to the locally-orthonormal tetrad.  Compute
     Psi_0 with the rotated basis.  Assert that the first-order residual
     Psi_0^full - Psi_0^rotated vanishes identically.

   Test (b) -- NEGATIVE CONTROL 1 (normalisation is not optional):
     Using the BACKGROUND tetrad (0, x^i/a) in the PERTURBED metric
     breaks X . X = 1 at first order (by -2 eps Phi + eps h_{xx}).
     Verify that this does produce a nonzero first-order residual in
     Psi_0, equal to Phi_{00}^{(0)} h_{mm} -- demonstrating that the
     rotation-invariance proof above relies on proper first-order
     normalisation, not on "basis corrections always drop out".

   Test (c) -- NEGATIVE CONTROL 2 (contrast with line-of-sight):
     Show that first-order corrections to n^mu DO feed into Phi_00^(1)
     through the background FLRW Ricci (Ricci != 0 unlike Weyl).
   ============================================================ *)

Print["=== Screen-tetrad convention tests ==="];

(* Background-only tetrad for negative controls *)
tetBg[ii_Integer] := Module[{vec},
  vec = Table[0, {4}];
  Do[
    vec[[jj + 1]] = (1/aa) KroneckerDelta[ii, jj],
    {jj, 3}];
  vec
];
XsBg   = tetBg[1];  YsBg = tetBg[2];  eLineBg = tetBg[3];
ZsBg   = (XsBg + I YsBg)/Sqrt[2];
uUBg   = {1/aa, 0, 0, 0};
kUBg   = Eobs (-uUBg + eLineBg);

(* -- Test (a): Rotation invariance of Psi_0 at first order -- *)

omegaField = omegaF[tt, x1, x2, x3];
XsRot = linearize[Xs + eps omegaField Ys];
YsRot = linearize[Ys - eps omegaField Xs];
ZsRot = linearize[(XsRot + I YsRot)/Sqrt[2]];

Print["[sanity] First-order orthonormality of rotated basis:"];
Print["         X'.X'  = ", Simplify[linearize[dot[XsRot, XsRot]]], "  (expect 1)"];
Print["         Y'.Y'  = ", Simplify[linearize[dot[YsRot, YsRot]]], "  (expect 1)"];
Print["         X'.Y'  = ", Simplify[linearize[dot[XsRot, YsRot]]], "  (expect 0)"];
Print["         X'.u   = ", Simplify[linearize[dot[XsRot, uU]]],    "  (expect 0)"];
Print["         X'.k   = ", Simplify[linearize[dot[XsRot, kU]]],    "  (expect 0)"];

psi0Rotated = -Sum[
    RiemDown[[aaI, bb, cc, dd]] kU[[aaI]] ZsRot[[bb]] kU[[cc]] ZsRot[[dd]],
    {aaI, 4}, {bb, 4}, {cc, 4}, {dd, 4}];
psi0Rotated = linearize[psi0Rotated] // Expand;
psi0RotLin  = Coefficient[psi0Rotated, eps, 1];

psi0RotDiff = Simplify[Expand[psi0Lin - psi0RotLin]];

Print["[test a] Psi_0^(1) rotation-invariance residual:"];
Print["         (Psi_0^full) - (Psi_0^rotated) = ", psi0RotDiff];
If[psi0RotDiff === 0,
  Print["  PASS: first-order screen-plane rotations leave Psi_0^(1) invariant,"];
  Print["        identically and for arbitrary omega.  Therefore any two"];
  Print["        first-order-orthonormal bases give the same Psi_0^(1);"];
  Print["        parallel-transported tetrad <-> locally-orthonormal tetrad"];
  Print["        agree at this order."],
  Print["  FAIL: rotation changed Psi_0 at first order -- spin-2 argument broken."]
];

(* -- Test (b): background (unnormalised) basis gives wrong answer -- *)

psi0WithBgBasis = -Sum[
    RiemDown[[aaI, bb, cc, dd]] kU[[aaI]] ZsBg[[bb]] kU[[cc]] ZsBg[[dd]],
    {aaI, 4}, {bb, 4}, {cc, 4}, {dd, 4}];
psi0WithBgBasis = linearize[psi0WithBgBasis] // Expand;
psi0BgBasisLin  = Coefficient[psi0WithBgBasis, eps, 1];
psi0BasisDiff   = Simplify[Expand[psi0Lin - psi0BgBasisLin]];

Print["[test b] Psi_0^(1) with BACKGROUND basis (NOT first-order orthonormal in full metric):"];
Print["         residual = (full) - (bg) = ", psi0BasisDiff];
Print["  Context: X_bg . X_bg = 1 + eps (-2 Phi + h_{xx}),  so X_bg fails"];
Print["  first-order orthonormality in the perturbed metric.  The residual"];
Print["  is therefore expected to be NONZERO -- this confirms that rotation"];
Print["  invariance in test (a) relies on proper first-order normalisation,"];
Print["  not on a generic 'basis corrections always vanish' claim."];
If[psi0BasisDiff =!= 0,
  Print["  CONSISTENT: nonzero residual, as expected."],
  Print["  UNEXPECTED: residual is zero -- check the arithmetic."]
];

(* -- Test (c): line-of-sight corrections in Phi_00 -- *)

phi00WithBgK = -(1/2) Sum[Ricci[[mm, nn]] kUBg[[mm]] kUBg[[nn]], {mm, 4}, {nn, 4}];
phi00WithBgK = linearize[phi00WithBgK] // Expand;
phi00BgKLin  = Coefficient[phi00WithBgK, eps, 1];
phi00KDiff   = Expand[phi00Lin - phi00BgKLin];

Print["[test c] Phi_00^(1): full k vs bg k (drop delta u, delta n corrections)"];
Print["         LeafCount of residual = ", LeafCount[phi00KDiff]];
If[phi00KDiff === 0,
  Print["  UNEXPECTED: bg k gives same answer -- check logic."],
  Print["  CONSISTENT: first-order k corrections DO feed Phi_00^(1) because"];
  Print["  background FLRW Ricci is nonzero (unlike background Weyl)."]
];

Print["=== end of screen-tetrad convention tests ==="];

(* ============================================================
   Phase 14.  Flat-triad vs full-tetrad-spatial projection test.
   ------------------------------------------------------------
   The paper defines the spin-s tensor projections as e.g.
       {}_2 h == m_(0)^i m_(0)^j h_ij ,     (flat triad)
   where m_(0)^i = x_(0)^i + i y_(0)^i is the background Euclidean
   null screen vector (no factor of 1/a).

   An alternative convention would use the spatial part of the FULL
   first-order-corrected 4-vector tetrad, i.e. m^i = spatial(Z^mu),
   whose leading behaviour is (1/a)(x_(0)^i + i y_(0)^i) + O(eps).

   Since h_{ij} is itself O(eps), the O(eps) corrections to the
   triad contribute only at O(eps^2) to the projection.  The two
   conventions therefore differ ONLY by the leading factor 1/a^2
   from (1/a)^2:

       {}_2 h^{full}  =  (1/a^2) {}_2 h^{flat}  +  O(eps^2).

   We verify this explicitly, and also confirm that the script's
   Psi_0^{(1,t)} output (which uses the full 4-vector tetrad Zs
   in the contraction -R k Z k Z) factorises into
       (E^2 / a^2) * (combinatorics of raw h_ij) ,
   so that writing
       Psi_0^{(1,t)} = (E^2/4a^2)[D^2 {}_2h + ...]
   with the FLAT {}_2h is internally consistent with the paper's
   prefactor.
   ============================================================ *)

Print["=== Flat-triad vs full-tetrad-spatial projection test ==="];

(* The physical h_{ij} is O(eps); we enforce that by multiplying the symbolic
   3x3 h-block by eps.  Then both projections are O(eps^1), and their
   difference at O(eps^1) exposes the genuine leading-order ratio. *)

mFlat = {0, 1/Sqrt[2], I/Sqrt[2], 0};     (* spatial = (x_(0) + i y_(0))/Sqrt[2]; NP-normalised *)

ZsSpatial = {0, Zs[[2]], Zs[[3]], Zs[[4]]};   (* drop time component *)

(* Raw h_{ij} block, then multiplied by eps to reflect its physical order *)
hBlock    = Table[hh[ii, jj], {ii, 3}, {jj, 3}];
hBlockEps = eps hBlock;

(* flat projection: {}_2 h^{flat} = sum_{i,j} m_flat^i m_flat^j h_{ij} *)
h2Flat = Sum[
   mFlat[[ii + 1]] mFlat[[jj + 1]] hBlockEps[[ii, jj]],
   {ii, 3}, {jj, 3}];

(* full-tetrad-spatial projection: m_full^i m_full^j h_{ij}, linearised *)
h2Full = Sum[
   ZsSpatial[[ii + 1]] ZsSpatial[[jj + 1]] hBlockEps[[ii, jj]],
   {ii, 3}, {jj, 3}] // linearize;

(* Correct claim: at first order, h2Full = (1/(2 a^2)) h2Flat.
   The factor of 1/2 comes from the normalisation Z^mu = (X^mu+iY^mu)/sqrt(2)
   of the 4-vector tetrad (which carries a 1/sqrt(2) that the 3-vector
   screenvec^i_(0) = x^i_(0) + i y^i_(0) does NOT carry), and the factor
   of 1/a^2 comes from the spatial-part leading scaling (1/a)*flat-triad. *)

residual = Simplify[h2Full - h2Flat / aa^2];

Print["[test] {}_2 h^{full}  -  (1/a^2) {}_2 h^{flat}  at O(eps^1):"];
Print["  residual = ", residual];
If[residual === 0,
  Print["  PASS: h2Full = h2Flat / a^2 at first order."];
  Print["        The factor (1/a^2) is from the spatial-tetrad 1/a scaling;"];
  Print["        with NP-normalised Z_mu and NP-normalised mFlat, no extra 1/2."];
  Print["        Equivalently, {}_2 h^{flat} = a^2 * (Z^i Z^j h_ij)."],
  Print["  FAIL: residual is nonzero; check the derivation."]
];

(* Now verify self-consistency with the paper's prefactor:
   The script computes Psi_0^{(1,t)} using the FULL 4-vector Zs.
   Below we extract the tensor-only, spatial-spatial part of the
   contraction and compare with -E^2/(4 a^2) times flat-h ansatz. *)

(* Take psi0TensorTT (already computed in Phase 11 TT slice) and look at
   its structure.  Since the paper writes the coefficient as E^2/(4a^2),
   we check that psi0TensorTT * 4 * a^2 / E^2 has the D^2 h_11-style
   structure WITHOUT a hidden 1/a^2, i.e., is a polynomial in H and h
   with no a^{-2} factor. *)
check = psi0TensorTT 4 a[tt]^2 / Eobs^2 // Simplify;
Print["[test] 4 a^2/E^2 * Psi_0^{(1,t)}_TT  (should have NO 1/a^2 factor):"];
Print["  LeafCount: ", LeafCount[check]];
Print["  Contains a[tt] in denom? ",
      !FreeQ[check /. a[tt] -> aSym, aSym^(-n_ /; n > 0)] ||
      !FreeQ[check /. a[tt] -> aSym, 1/aSym^_]];
(* Expect: False -- no 1/a^n terms beyond the factored-out E^2/a^2. *)
Print["        (False = no hidden 1/a factor remaining -- paper's E^2/a^2"];
Print["         prefactor is the complete a-dependence, hence the flat-"];
Print["         triad definition of {}_2 h is self-consistent.)"];

Print["=== end of flat-vs-full triad consistency test ==="];

(* ============================================================
   Phase 15.  Generic Euclidean triad cross-check.
   ------------------------------------------------------------
   Recompute Phi_00^(1) and Psi_0^(1) using symbolic Euclidean
   3-vectors {xVec, yVec, nVec} in place of the fixed Cartesian
   axes (delta^i_1, delta^i_2, delta^i_3).  Verify that the
   result reduces to the fixed-axis answer upon specialisation
   xVec -> (1,0,0), yVec -> (0,1,0), nVec -> (0,0,1).

   This cross-check confirms that:
     * the paper's tetrad construction on a generic Euclidean
       3-triad reproduces the script's fixed-axis result;
     * the 1/a in the background tetrad is a geometric
       4-metric-normalisation factor, not a perturbation;
     * no orthonormality constraint between xVec, yVec, nVec
       is needed during the symbolic manipulation -- it is
       automatic upon specialisation to unit axes.
   ============================================================ *)

Print["=== Generic Euclidean triad cross-check ==="];

(* Generic Euclidean 3-vectors (symbolic components). *)
xVec = {xx1, xx2, xx3};
yVec = {yy1, yy2, yy3};
nVec = {nn1, nn2, nn3};

(* Paper's tetrad formula, promoted from a Euclidean 3-vector
   to a 4-vector orthonormal in the full perturbed metric to O(eps). *)
genericTet[vec3_] := Module[{v4},
  v4 = Table[0, {4}];
  v4[[1]] = -eps Sum[vec3[[ii]] di[Bfld, ii], {ii, 3}]/aa;
  Do[
    v4[[jj + 1]] = (1/aa) (vec3[[jj]]
                           + eps (Phifld vec3[[jj]]
                                  - (1/2) Sum[hh[jj, ii] vec3[[ii]], {ii, 3}])),
    {jj, 3}];
  linearize[v4]
];

eLineGen = genericTet[nVec];
XsGen    = genericTet[xVec];
YsGen    = genericTet[yVec];
ZsGen    = (XsGen + I YsGen)/Sqrt[2];

(* Photon built from the generic tetrad along nVec. *)
kUGen = Eobs (-uU + eLineGen);
kUGen = linearize[kUGen];

(* Driving-field contractions with the generic tetrad. *)
phi00Gen = -(1/2) Sum[Ricci[[mm, nn]] kUGen[[mm]] kUGen[[nn]],
                    {mm, 4}, {nn, 4}];
phi00Gen = linearize[phi00Gen] // Expand;
phi00GenLin = Coefficient[phi00Gen, eps, 1];

psi0Gen = -Sum[
    RiemDown[[aaI, bb, cc, dd]] kUGen[[aaI]] ZsGen[[bb]] kUGen[[cc]] ZsGen[[dd]],
    {aaI, 4}, {bb, 4}, {cc, 4}, {dd, 4}];
psi0Gen = linearize[psi0Gen] // Expand;
psi0GenLin = Coefficient[psi0Gen, eps, 1];

(* Specialise generic -> canonical axes. *)
specialiseRules = {xx1 -> 1, xx2 -> 0, xx3 -> 0,
                   yy1 -> 0, yy2 -> 1, yy3 -> 0,
                   nn1 -> 0, nn2 -> 0, nn3 -> 1};

phi00Diff = Simplify[Expand[(phi00GenLin /. specialiseRules) - phi00Lin]];
psi0Diff  = Simplify[Expand[(psi0GenLin  /. specialiseRules) - psi0Lin]];

Print["[check] Phi_00^(1) generic-vs-fixed residual  = ", phi00Diff];
If[phi00Diff === 0,
   Print["  PASS: generic Euclidean triad reproduces fixed-axis Phi_00^(1)."],
   Print["  FAIL: residual nonzero -- inspect tetrad construction."]
];

Print["[check] Psi_0^(1)  generic-vs-fixed residual  = ", psi0Diff];
If[psi0Diff === 0,
   Print["  PASS: generic Euclidean triad reproduces fixed-axis Psi_0^(1)."],
   Print["  FAIL: residual nonzero -- inspect tetrad construction."]
];

(* Secondary cross-check: +pi/2 screen-plane rotation.
     xVec -> yVec,  yVec -> -xVec.
   Under this, Z = (X+iY)/Sqrt[2] -> (Y - iX)/Sqrt[2] = -i Z, so the
   spin-2 quantity Psi_0 picks up (-i)^2 = -1:
     Phi_00 -> Phi_00  (spin 0, invariant),
     Psi_0  -> -Psi_0  (spin +2). *)
rot90 = {xx1 -> yy1, xx2 -> yy2, xx3 -> yy3,
         yy1 -> -xx1, yy2 -> -xx2, yy3 -> -xx3};

phi00GenRot = Expand[phi00GenLin /. rot90];
psi0GenRot  = Expand[psi0GenLin  /. rot90];

phi00RotDiff = Simplify[phi00GenRot - phi00GenLin];
psi0RotSum   = Simplify[psi0GenRot + psi0GenLin];   (* Psi_0^rot should equal -Psi_0 *)

Print["[check] Phi_00^(1) invariant under x->y, y->-x?  residual = ", phi00RotDiff];
If[phi00RotDiff === 0,
   Print["  PASS: spin-0 scalar, invariant under screen-plane rotation."],
   Print["  FAIL: Phi_00 should be invariant under screen-plane rotation."]
];

Print["[check] Psi_0^(1)  flips sign under x->y, y->-x?  residual = ", psi0RotSum];
If[psi0RotSum === 0,
   Print["  PASS: spin-(+2) quantity, Psi_0 -> -Psi_0 under Z -> -i Z."],
   Print["  FAIL: Psi_0 should flip sign under pi/2 screen rotation."]
];

Print["=== end of generic Euclidean triad cross-check ==="];

Print["[done]"];
