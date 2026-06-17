(* ============================================================
   driving_fields_harmonics.wl
   ------------------------------------------------------------
   Spherical-harmonic / spin-weighted transform of the Poisson-
   gauge driving fields

        Phi_00 = (1/2) R_{mu nu} k^mu k^nu
        Psi_0  = - R_{abcd} k^a Z^b k^c Z^d

   For a scalar field X(tau, chi * n_hat), expand on the sphere

        X(tau, chi n) = Sum_{ell, m} X_{ell m}(tau, chi) Y_{ell m}(n).

   The driving field Phi_00 (a scalar observable) is likewise expanded
   in Y_{ell m}, while Psi_0 (a complex spin-2 observable) is expanded
   in the spin-+2 harmonics {}_2 Y_{ell m}:

        Phi_00(tau, chi, n) = Sum  Phi00_{ell m}(tau, chi) Y_{ell m}(n),
        Psi_0 (tau, chi, n) = Sum  Psi0_{ell m} (tau, chi)
                                      {}_2 Y_{ell m}(n).

   Using radial/angular decomposition of partial derivatives

        d_i  = n_i d_r + (1/r) nabla_i^{(S)}
        n^i d_i f         = d_chi f
        n^i n^j d_i d_j f = d_chi^2 f
        Lap f_{ell m}     = [d_chi^2 + (2/chi) d_chi - L2/chi^2] f_{ell m}
        Lap_perp f_{ell m}= [(2/chi) d_chi - L2/chi^2] f_{ell m}
        m^i d_i f         = (1/chi) edth f           (spin-raising)
        m^i m^j d_i d_j f = (1/chi^2) edth^2 f

   with eigenvalues on Y_{ell m}:
        edth^2 Y_{ell m} = Sqrt[L2 (L2 - 2)]  {}_2 Y_{ell m},
        L2 := ell (ell + 1).

   The script produces the per-multipole transfer functions

        Phi00_{ell m}(tau, chi) = Sum_X  T^{Phi_00}_{X, ell}(tau, chi, d_tau, d_chi) X_{ell m}(tau, chi),
        Psi0_{ell m} (tau, chi) = Sum_X  T^{Psi_0}_{X, ell} (tau, chi, d_tau, d_chi) X_{ell m}(tau, chi),

   where X ranges over the scalar perturbations {Psi, Phi, B}.  The
   angular power spectra C_ell(chi, chi') of the driving fields then
   follow by acting with the operators T on the redshift-redshift
   angular power spectra C_ell^{XY}(chi, chi') of the perturbations.

   Run with:
     ~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh \
        scripts/driving_fields_harmonics.wl 300
   ============================================================ *)

$HistoryLength = 0;

(* --- Helpers (operator actions on a radial multipole f[tt, chi]) --- *)

(* L2 represents ell(ell+1); kept symbolic *)
dchi[f_]  := D[f, chi];
dtt[f_]   := D[f, tt];
dchi2[f_] := D[f, {chi, 2}];
dtt2[f_]  := D[f, {tt, 2}];
dttchi[f_] := D[f, tt, chi];

(* Physical 3D operators reduced to radial + ell-eigenvalue forms. *)
LaplaceFull[f_]  := dchi2[f] + (2/chi) dchi[f] - (L2/chi^2) f;
nini[f_]         := dchi2[f];
nPerp[f_]        := (2/chi) dchi[f] - (L2/chi^2) f;
nPar[f_]         := dchi[f];
Dnull[f_]        := dtt[f] - dchi[f];
Dnull2[f_]       := dtt2[f] - 2 dttchi[f] + dchi2[f];

(* Spin-raising (twice) coefficient, so that m^i m^j d_i d_j acting on
   a scalar multipole gives the coefficient in front of {}_2 Y_{ell m}.
   With NP-normalised m = (x + i y)/Sqrt[2], the factor of 1/2 appears. *)
spinRaise2[f_]   := (1/(2 chi^2)) Sqrt[L2 (L2 - 2)] f;

(* --- Scalar perturbation multipoles (placeholders) --- *)

PsiL = psiLM[tt, chi];
PhiL = phiLM[tt, chi];
BL   = bLM[tt, chi];

(* Conformal Hubble (background) *)
Hc  = Hconf[tt];
Hcp = Derivative[1][Hconf][tt];

(* --- Per-multipole driving fields ---
   Input: first-order Poisson-gauge expressions from driving_fields_poisson.wl.

   Scalar sector:
     (a^2/E^2) Phi_00^{(1,s)} = 2(H' - H^2) Psi + H(Psi' - Phi') - 2 H n^i d_i Psi
                              + (1/2) Lap_perp (Phi + Psi) + D^2 Phi
     (a^2/E^2) Psi_0 ^{(1,s)} = -(1/2) m^i m^j d_i d_j (Phi + Psi)

   B sector:
     (a^2/E^2) Phi_00^{(1,B)} = 2(H^2 - H') n^i d_i B + H n^i n^j d_i d_j B - (1/2) Lap_perp B'
     (a^2/E^2) Psi_0 ^{(1,B)} = (1/2) m^i m^j d_i d_j B'                                                     *)

(* Paper sign convention: Phi_00 = -(1/2) R_{mu nu} k^mu k^nu, so the
   scalar-sector Phi_00 expression carries an overall minus. *)
Phi00ScalarLM = -(
    2 (Hcp - Hc^2) PsiL
  + Hc (dtt[PsiL] - dtt[PhiL])
  - 2 Hc nPar[PsiL]
  + (1/2) nPerp[PhiL + PsiL]
  + Dnull2[PhiL]
);

Psi0ScalarLM = -spinRaise2[PhiL + PsiL];

Phi00BLM = -(
    2 (Hc^2 - Hcp) nPar[BL]
  + Hc nini[BL]
  - (1/2) nPerp[dtt[BL]]
);

Psi0BLM = spinRaise2[dtt[BL]];

(* --- Collect as operator-form transfer functions ---
   We isolate the coefficient of each scalar multipole (and each of its
   partial derivatives) so the user can read off the differential
   operator T^{driving}_{X, ell}(tau, chi, d_tau, d_chi). *)

(* (TeX paper-notation replacement table removed; raw TeXForm suffices.) *)

(* --- Emit --- *)
Print["%%TEX_PHI00_SCALAR_LM_START%%"];
Print[ToString[TeXForm[Phi00ScalarLM]]];
Print["%%TEX_PHI00_SCALAR_LM_END%%"];

Print["%%TEX_PSI0_SCALAR_LM_START%%"];
Print[ToString[TeXForm[Psi0ScalarLM]]];
Print["%%TEX_PSI0_SCALAR_LM_END%%"];

Print["%%TEX_PHI00_B_LM_START%%"];
Print[ToString[TeXForm[Phi00BLM]]];
Print["%%TEX_PHI00_B_LM_END%%"];

Print["%%TEX_PSI0_B_LM_START%%"];
Print[ToString[TeXForm[Psi0BLM]]];
Print["%%TEX_PSI0_B_LM_END%%"];

(* --- Transfer-operator coefficients ---
   For each field X, collect Phi00 and Psi0 per multipole as a sum of
   derivative operators applied to X_{ell m}.  We expose these as Coefficient[]
   extractions against the basis of derivatives. *)

basisFor[f_] := {
  f,
  dchi[f],
  dchi2[f],
  dtt[f],
  dtt2[f],
  dttchi[f]
};

opMatrix[expr_, basis_] :=
  Table[Coefficient[Expand[expr], b], {b, basis}];

coefPhi00ScalarPsi = opMatrix[Phi00ScalarLM, basisFor[PsiL]];
coefPhi00ScalarPhi = opMatrix[Phi00ScalarLM, basisFor[PhiL]];
coefPsi0ScalarPsi  = opMatrix[Psi0ScalarLM,  basisFor[PsiL]];
coefPsi0ScalarPhi  = opMatrix[Psi0ScalarLM,  basisFor[PhiL]];
coefPhi00BB        = opMatrix[Phi00BLM,      basisFor[BL]];
coefPsi0BB         = opMatrix[Psi0BLM,       basisFor[BL]];

Print["%%COEF_PHI00_SCALAR_PSI_START%%"];
Print[ToString[TeXForm[coefPhi00ScalarPsi]]];
Print["%%COEF_PHI00_SCALAR_PSI_END%%"];

Print["%%COEF_PHI00_SCALAR_PHI_START%%"];
Print[ToString[TeXForm[coefPhi00ScalarPhi]]];
Print["%%COEF_PHI00_SCALAR_PHI_END%%"];

Print["%%COEF_PSI0_SCALAR_PSI_START%%"];
Print[ToString[TeXForm[coefPsi0ScalarPsi]]];
Print["%%COEF_PSI0_SCALAR_PSI_END%%"];

Print["%%COEF_PSI0_SCALAR_PHI_START%%"];
Print[ToString[TeXForm[coefPsi0ScalarPhi]]];
Print["%%COEF_PSI0_SCALAR_PHI_END%%"];

Print["%%COEF_PHI00_B_START%%"];
Print[ToString[TeXForm[coefPhi00BB]]];
Print["%%COEF_PHI00_B_END%%"];

Print["%%COEF_PSI0_B_START%%"];
Print[ToString[TeXForm[coefPsi0BB]]];
Print["%%COEF_PSI0_B_END%%"];

(* --- Angular power spectrum construction (symbolic) ---
   Given two linear operators T_A, T_B acting on the same underlying
   scalar multipole X_{ell m}(tau, chi), the cross-power at two
   light-cone distances (chi, chi') is:

        C_ell^{AB}(chi, chi') =
            T_A(tau, chi, d_tau, d_chi)  T_B(tau', chi', d_tau', d_chi')
            C_ell^{XX}(tau, chi; tau', chi').

   On the observer's past light cone we set tau  = tau_0 - chi,
                                             tau' = tau_0 - chi'.
   We emit the OPERATOR acting symbolically on a covariance kernel
   K[t, c, T, C] (placeholder for C_ell^{XX}(tau, chi; tau', chi')).
*)

(* Build a symbolic covariance kernel as a function of
   (tau_A = tau, chi_A = chi, tau_B = TauP, chi_B = chiP). *)

Clear[kernK];
applyA[op_, K_] := op /. {
  dchi[x_]  :> D[K, chi],  dtt[x_]  :> D[K, tt],
  dchi2[x_] :> D[K, {chi, 2}], dtt2[x_] :> D[K, {tt, 2}],
  dttchi[x_] :> D[K, tt, chi],
  psiLM[tt, chi] -> K, phiLM[tt, chi] -> K, bLM[tt, chi] -> K
};

(* Demo: Phi_00 x Phi_00 cross-spectrum contribution from Psi x Psi. *)
ClearAll[TpAChi, TpAPhi, TpBChi, TpBPhi];

(* Operator form acting on the "A-side" (tau, chi) and "B-side" (tauP, chiP). *)

sideA = Phi00ScalarLM;
sideB = Phi00ScalarLM /. {tt -> tauP, chi -> chiP,
                          psiLM -> psiLMp, phiLM -> phiLMp, bLM -> bLMp};

(* Placeholder: cross-correlation of underlying multipoles is
   < X_{ell m}(tau, chi) X'^*_{ell m}(tauP, chiP) > = KXY[tt, chi, tauP, chiP]. *)

(* Psi x Psi contribution *)
ABkernelPsiPsi =
  sideA sideB /. {
    psiLM[tt, chi] * psiLMp[tauP, chiP] -> KPsiPsi[tt, chi, tauP, chiP],
    Derivative[a_, b_][psiLM][tt, chi] * Derivative[c_, d_][psiLMp][tauP, chiP]
      :> D[KPsiPsi[tt, chi, tauP, chiP], {tt, a}, {chi, b}, {tauP, c}, {chiP, d}]
  };

(* Note: the real expansion is a bilinear in all the scalar fields; the
   above is illustrative.  The user supplies the 3x3 matrix of
   C_ell^{XY}(chi, chi') covariance kernels and combines with the
   transfer operators above. *)

(* --- Emit the six transfer-operator bases more readably --- *)

printOp[label_, coefs_, fieldSymbol_, derivs_] := (
  Print[label];
  Do[
    If[coefs[[i]] =!= 0,
      Print["  ", derivs[[i]], " : ", coefs[[i]]]],
    {i, Length[derivs]}]
);

derivLabels = {"X_{lm}", "d_chi X_{lm}", "d_chi^2 X_{lm}",
               "d_tt X_{lm}", "d_tt^2 X_{lm}", "d_tt d_chi X_{lm}"};

Print["=== Transfer operators (ell, m implicit; L2 = ell(ell+1)) ==="];
printOp["Phi_00[Psi-sector]:", coefPhi00ScalarPsi, "Psi", derivLabels];
printOp["Phi_00[Phi-sector]:", coefPhi00ScalarPhi, "Phi", derivLabels];
printOp["Phi_00[B-sector]:",   coefPhi00BB,        "B",   derivLabels];
printOp["Psi_0[Psi-sector]:",  coefPsi0ScalarPsi,  "Psi", derivLabels];
printOp["Psi_0[Phi-sector]:",  coefPsi0ScalarPhi,  "Phi", derivLabels];
printOp["Psi_0[B-sector]:",    coefPsi0BB,         "B",   derivLabels];

(* --- Tensor sector of Phi_00 and Psi_0 ---

   The tensor perturbation h_{ij}(tau, vec x) is transverse-traceless with 2
   physical polarisation modes.  In harmonic space it is conventionally
   decomposed into gradient (E) and curl (B) modes with respect to the
   celestial sphere.  Acting on an observer-centred shell, a single real
   TT tensor mode admits an expansion in spin-+/-2 harmonics,

        h_{mm}(tau, chi n)  = Sum_{ell m}  h_{ell m}^{(+)}(tau, chi) {}_2 Y_{ell m}(n),
        h_{mn}(tau, chi n)  = Sum_{ell m}  h_{ell m}^{(1)}(tau, chi) {}_1 Y_{ell m}(n),
        h_{nn}(tau, chi n)  = Sum_{ell m}  h_{ell m}^{(0)}(tau, chi)      Y_{ell m}(n),

   where h_{ell m}^{(0,1,2)} are linear combinations of the standard
   E/B tensor-mode harmonic amplitudes at each (tau, chi).

   The driving-field tensor contributions from cosmology.tex, Eqs.
   (Phi00 tensor) and (Psi0 tensor), then become (using D = d_tau - d_chi,
   L2 = ell(ell+1), NP-normalised hat-screenvec = (x+iy)/Sqrt[2]):

        Phi_00^{(1,t)}_{ell m}
          = (E^2 / 4 a^2)  [d_tt^2 + 2 H d_tt + (L2/chi^2) - d_chi^2 - (2/chi) d_chi] h_{ell m}^{(0)},

        Psi_0 ^{(1,t)}_{ell m}
          = (E^2 / 2 a^2) [  D^2 h_{ell m}^{(+)}
                           + 2 D (1/(chi Sqrt[2])) Sqrt[L2 - 2] h_{ell m}^{(1)}
                           + (1/(2 chi^2)) Sqrt[L2 (L2-2)] h_{ell m}^{(0)} ].

   Spin-raising eigenvalues (Goldberg et al., NP-normalised m):
        edth {}_s Y_{l m} = + sqrt((l - s)(l + s + 1)) {}_{s+1} Y_{l m}.
   So spin-0 -> spin-1 brings sqrt(L2), spin-1 -> spin-2 brings sqrt(L2 - 2),
   and the double raise spin-0 -> spin-2 brings sqrt(L2 (L2 - 2)).  The
   coefficient in front of h_{ell m}^{(1)} must therefore be sqrt(L2 - 2),
   NOT sqrt(L2).  (Verified numerically in /tmp/verify_spin_raise.wl.)

   These are printed below symbolically. *)

h0LM = hL0[tt, chi];   (* h_{nn} spin-0 radial multipole *)
h1LM = hL1[tt, chi];   (* h_{mn} spin-1 radial multipole *)
h2LM = hL2[tt, chi];   (* h_{mm} spin-2 radial multipole *)

Phi00TensorLM = -(1/4) (
      dtt2[h0LM] + 2 Hc dtt[h0LM]
    - dchi2[h0LM] - (2/chi) dchi[h0LM] + (L2/chi^2) h0LM
);

Psi0TensorLM = (1/2) (
      Dnull2[h2LM]
    + 2 Dnull[ (1/(chi Sqrt[2])) Sqrt[L2 - 2] h1LM ]
    + (1/(2 chi^2)) Sqrt[L2 (L2 - 2)] h0LM
);

Print["%%TEX_PHI00_TENSOR_LM_START%%"];
Print[ToString[TeXForm[Phi00TensorLM]]];
Print["%%TEX_PHI00_TENSOR_LM_END%%"];

Print["%%TEX_PSI0_TENSOR_LM_START%%"];
Print[ToString[TeXForm[Psi0TensorLM]]];
Print["%%TEX_PSI0_TENSOR_LM_END%%"];

coefPhi00T0 = opMatrix[Phi00TensorLM, basisFor[h0LM]];
coefPsi0T0  = opMatrix[Psi0TensorLM,  basisFor[h0LM]];
coefPsi0T1  = opMatrix[Psi0TensorLM,  basisFor[h1LM]];
coefPsi0T2  = opMatrix[Psi0TensorLM,  basisFor[h2LM]];

printOp["Phi_00[h^{(0)}-sector]:", coefPhi00T0, "h0", derivLabels];
printOp["Psi_0[h^{(0)}-sector]:",  coefPsi0T0,  "h0", derivLabels];
printOp["Psi_0[h^{(1)}-sector]:",  coefPsi0T1,  "h1", derivLabels];
printOp["Psi_0[h^{(2)}-sector]:",  coefPsi0T2,  "h2", derivLabels];

(* ============================================================
   Part 12.  Re-express the tensor sector in the E/B basis.
   ============================================================

   Physical GW modes on each observer-centred shell can be parametrised by
   two parity-distinguished amplitudes per (ell, m):

        h^E_{ell m}(tau, chi)   (parity-even, "gradient"-type)
        h^B_{ell m}(tau, chi)   (parity-odd,  "curl"-type)

   The spin-+s projections of the underlying TT tensor on the sphere are
   related to (h^E, h^B) by the standard total-angular-momentum (TAM)
   decomposition.  Retaining the leading ell-dependent algebraic form (the
   full radial transfer functions are summarised in the appendix), we use

        h^{(2)}_{ell m}(tau, chi) = h^E_{ell m}(tau, chi)
                                    + I h^B_{ell m}(tau, chi),

        h^{(1)}_{ell m}(tau, chi) = alpha1E[ll, chi] h^E_{ell m}(tau, chi)
                                    + I alpha1B[ll, chi] h^B_{ell m}(tau, chi),

        h^{(0)}_{ell m}(tau, chi) = alpha0E[ll, chi] h^E_{ell m}(tau, chi),

   where alpha0E, alpha1E, alpha1B are ell- (and possibly chi-) dependent
   algebraic factors fixed by the 3D TT conditions.  The spin-0 projection
   receives no B contribution by parity.  We keep these factors symbolic
   so that the transfer operators are exact once alpha is supplied. *)

hELM = hE[tt, chi];
hBLM = hB[tt, chi];

(* Conversion factors.  Given explicit form in appendix / conventions. *)
alpha0E[ll_, chi_] := A0E[ll, chi];
alpha1E[ll_, chi_] := A1E[ll, chi];
alpha1B[ll_, chi_] := A1B[ll, chi];

h0InEB = alpha0E[ll, chi] hELM;
h1InEB = alpha1E[ll, chi] hELM + I alpha1B[ll, chi] hBLM;
h2InEB = hELM + I hBLM;

(* Substitute into the tensor transfer operators *)
Phi00TensorEB = (Phi00TensorLM /. {hL0[tt, chi] -> h0InEB,
                                    Derivative[a_, b_][hL0][tt, chi] :>
                                       D[h0InEB, {tt, a}, {chi, b}]});
Psi0TensorEB  = Psi0TensorLM /. {
  hL0[tt, chi] -> h0InEB,
  Derivative[a_, b_][hL0][tt, chi] :> D[h0InEB, {tt, a}, {chi, b}],
  hL1[tt, chi] -> h1InEB,
  Derivative[a_, b_][hL1][tt, chi] :> D[h1InEB, {tt, a}, {chi, b}],
  hL2[tt, chi] -> h2InEB,
  Derivative[a_, b_][hL2][tt, chi] :> D[h2InEB, {tt, a}, {chi, b}]
};

Phi00TensorEB = Expand[Phi00TensorEB];
Psi0TensorEB  = Expand[Psi0TensorEB];

(* Split into E- and B-contributions (Psi_0 is complex; Phi_00 only has E) *)
phi00TensorE = Phi00TensorEB;
psi0TensorE  = Psi0TensorEB /. {hB[tt, chi] -> 0,
                                 Derivative[__][hB][__] -> 0} // Expand;
psi0TensorB  = Psi0TensorEB - psi0TensorE // Expand;

Print["%%TEX_PHI00_TENSOR_EB_START%%"];
Print[ToString[TeXForm[Phi00TensorEB]]];
Print["%%TEX_PHI00_TENSOR_EB_END%%"];

Print["%%TEX_PSI0_TENSOR_E_START%%"];
Print[ToString[TeXForm[psi0TensorE]]];
Print["%%TEX_PSI0_TENSOR_E_END%%"];

Print["%%TEX_PSI0_TENSOR_B_START%%"];
Print[ToString[TeXForm[psi0TensorB]]];
Print["%%TEX_PSI0_TENSOR_B_END%%"];

coefPhi00TensorE = opMatrix[phi00TensorE, basisFor[hELM]];
coefPsi0TensorE  = opMatrix[psi0TensorE,  basisFor[hELM]];
coefPsi0TensorB  = opMatrix[psi0TensorB /. (I x_) :> x, basisFor[hBLM]];

printOp["Phi_00[E-sector]:", coefPhi00TensorE, "E", derivLabels];
printOp["Psi_0[E-sector]:",  coefPsi0TensorE,  "E", derivLabels];
printOp["Psi_0[B-sector, divided by i]:", coefPsi0TensorB, "B", derivLabels];

(* --- Save a machine-readable dump of all transfer operators --- *)
Put[{
   "Phi00_Psi"  -> coefPhi00ScalarPsi,
   "Phi00_Phi"  -> coefPhi00ScalarPhi,
   "Phi00_B"    -> coefPhi00BB,
   "Psi0_Psi"   -> coefPsi0ScalarPsi,
   "Psi0_Phi"   -> coefPsi0ScalarPhi,
   "Psi0_B"     -> coefPsi0BB,
   "Phi00_h0"   -> coefPhi00T0,
   "Psi0_h0"    -> coefPsi0T0,
   "Psi0_h1"    -> coefPsi0T1,
   "Psi0_h2"    -> coefPsi0T2,
   "derivBasis" -> derivLabels
  },
  "scripts/driving_fields_harmonics.m"];
Print["Saved transfer operators to scripts/driving_fields_harmonics.m"];

Print["[done]"];
