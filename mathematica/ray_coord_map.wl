(* ============================================================
   ray_coord_map.wl
   ------------------------------------------------------------
   Symbolic derivation of the past-light-cone coordinate map
     R : (lambda, Omega)  -->  (chi(lambda, Omega), n_hat(lambda, Omega))
   on the perturbed Poisson-gauge FLRW metric, using xAct/xCoba
   coordinate-component calculus.

   Deliverables (printed between %%TEX_*%% markers for pasting into
   sections/cosmology.tex):

     (1) Background radial rate  d chi / d lambda = EE/a = E_0/a^2
         for past-directed photons (EE = local energy, EE = E_0/a).
     (2) Background Jacobian  J_(0) = d(chi, n_hat^A)/d(lambda, Omega^B)
         and its determinant  det J_(0) = EE/a = E_0/a^2.
     (3) Chain rule  d chi/dz = 1/Hconf(z), hence
         chi(z) = int_0^z dz'/Hconf(z').
     (4) First-order null-geodesic RHS  d(delta k^mu)/d lambda
         showing structure of the line-of-sight integrands driving
         delta chi and delta n_hat at O(eps).  The user paper then
         states these *structurally* (paper is Born-level only).

   Conventions (matching affine_to_redshift_general.wl):
     ds^2 = a^2(tau) [-(1 + 2 eps Psi) d tau^2
                      - 2 eps (d_i B) d tau dx^i
                      + ((1 - 2 eps Phi) delta_{ij} + eps h_{ij}) dx^i dx^j]
     u^mu = (1/(a Sqrt[1 + 2 eps Psi]), 0)           (Hubble-flow observer)
     k^mu = EE * (-u^mu + e^mu)                      (past-directed)
     EE   = local photon energy (kept symbolic; equals E_0 at lambda=0,
            redshifted to E_0/a along the background ray).
     E_0  = observer-today energy (constant).

   Run with:
     ~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh \
        scripts/ray_coord_map.wl 300
   ============================================================ *)

$HistoryLength = 0;

(* ---- Phase 1.  Coordinates and perturbation fields ---- *)
coords = {tt, x1, x2, x3};

Psifld = PsiF[tt, x1, x2, x3];
Phifld = PhiF[tt, x1, x2, x3];
Bfld   = BF  [tt, x1, x2, x3];

hh[1, 1] := hF11[tt, x1, x2, x3];
hh[2, 2] := hF22[tt, x1, x2, x3];
hh[3, 3] := hF33[tt, x1, x2, x3];
hh[1, 2] := hF12[tt, x1, x2, x3]; hh[2, 1] := hh[1, 2];
hh[1, 3] := hF13[tt, x1, x2, x3]; hh[3, 1] := hh[1, 3];
hh[2, 3] := hF23[tt, x1, x2, x3]; hh[3, 2] := hh[2, 3];

di[f_, k_Integer] := D[f, coords[[k + 1]]];

(* ---- Phase 2.  Metric in Poisson gauge ---- *)
aa = a[tt];

gM = Table[0, {4}, {4}];
gM[[1, 1]] = -aa^2 (1 + 2 eps Psifld);
Do[ gM[[1, k + 1]] = -aa^2 eps di[Bfld, k];
    gM[[k + 1, 1]] = gM[[1, k + 1]], {k, 3}];
Do[ gM[[i + 1, j + 1]] = aa^2 ((1 - 2 eps Phifld) KroneckerDelta[i, j]
                               + eps hh[i, j]),
    {i, 3}, {j, 3}];

linearize[expr_] := Normal @ Series[expr, {eps, 0, 1}];

g0M    = gM /. eps -> 0;
g0Minv = Inverse[g0M];
dgMet  = gM - g0M;
gMinv  = g0Minv - g0Minv . dgMet . g0Minv;
gMinv  = linearize[gMinv];

Print["[progress] Christoffels ..."];
dgM = Table[D[gM[[bb, nn]], coords[[mm]]], {mm, 4}, {bb, 4}, {nn, 4}];
Gam = Table[
   (1/2) Sum[
     gMinv[[aaI, bb]] (dgM[[mm, bb, nn]] + dgM[[nn, bb, mm]] - dgM[[bb, mm, nn]]),
     {bb, 4}],
   {aaI, 4}, {mm, 4}, {nn, 4}];
Gam = linearize[Gam];

(* ---- Phase 3.  Observer direction Omega = (theta, phi)  ---- *)
(* Parameterise the sky direction on S^2.  The 3-unit vector is        *)
(*    Omega^i = (Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th])             *)

OmegaHat = {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]};

(* Sanity: |Omega| = 1 *)
omegaNorm = Sum[OmegaHat[[ii]]^2, {ii, 3}] // FullSimplify;
Print["[check] |Omega|^2 = ", omegaNorm, "  (expected 1)"];

(* ---- Phase 4.  Observer 4-velocity u^mu and line-of-sight tetrad e^mu ---- *)
(* Following affine_to_redshift_general.wl convention (past-directed        *)
(* photon).  e^mu is the first-order tetrad vector oriented along Omega.    *)

uUp = {1/(aa Sqrt[1 + 2 eps Psifld]), 0, 0, 0};
uUp = linearize[uUp];

uLo = Table[Sum[gM[[mm, nn]] uUp[[nn]], {nn, 4}], {mm, 4}];
uLo = linearize[uLo];

(* Tetrad along Omega: time component carries a B-shift correction;       *)
(* spatial components carry (Phi, h) corrections that normalise n.n = 1.  *)
tetOmega = Module[{vec},
   vec = Table[0, {4}];
   vec[[1]] = -eps Sum[OmegaHat[[ii]] di[Bfld, ii], {ii, 3}]/aa;
   Do[
     vec[[jj + 1]] =
       (1/aa) (OmegaHat[[jj]]
               + eps (Phifld OmegaHat[[jj]]
                      - (1/2) Sum[hh[jj, ii] OmegaHat[[ii]], {ii, 3}])),
     {jj, 3}];
   linearize[vec]
   ];
nUp = tetOmega;

(* Sanity checks on the tetrad *)
unDot = Sum[gM[[mm, nn]] uUp[[mm]] nUp[[nn]], {mm, 4}, {nn, 4}];
unDot = linearize[unDot];
Print["[check] u.n (should be 0) = ", Simplify @ unDot];

nnDot = Sum[gM[[mm, nn]] nUp[[mm]] nUp[[nn]], {mm, 4}, {nn, 4}];
nnDot = linearize[nnDot];
Print["[check] n.n (should be 1) = ", Simplify @ nnDot];

(* ---- Phase 5.  Photon four-momentum  k^mu = EE(-u^mu + n^mu)  ---- *)
kUp = EE (-uUp + nUp);
kUp = linearize[kUp];

kNorm = Sum[gM[[mm, nn]] kUp[[mm]] kUp[[nn]], {mm, 4}, {nn, 4}];
kNorm = linearize[kNorm];
Print["[check] k.k (should be 0) = ", Simplify @ kNorm];

(* ---- Phase 6.  Background radial rate d chi / d lambda ---- *)
(* On the past light cone in flat FLRW, chi = tau_0 - tau (conformal),   *)
(* so d chi/d lambda = - d tau / d lambda = - k^0.                        *)

kBg0 = kUp[[1]] /. eps -> 0 // Simplify;
dChiDlamBg = -kBg0 // Simplify;
Print["[bg] k^0       = ", kBg0,     "  (expected -EE/a)"];
Print["[bg] d chi/dl  = ", dChiDlamBg, "  (expected  EE/a = E_0/a^2)"];

(* ---- Phase 7.  Background Jacobian J_(0) ---- *)
(* Express n_hat in sphere angular coords (th, ph): at background,        *)
(* n_hat(lambda, Omega) = Omega, so d n_hat^A / d Omega^B = delta^A_B.    *)
(* d n_hat^A / d lambda = 0 (parallel ray at zeroth order).               *)
(* d chi / d Omega^A    = 0 (chi is Omega-independent at bg).             *)

JBg = {{dChiDlamBg, 0, 0},
       {0,           1, 0},
       {0,           0, 1}};

detJBg = Det[JBg] // Simplify;
Print["[bg] J_(0) = ", JBg // MatrixForm];
Print["[bg] det J_(0) = ", detJBg, "  (expected EE/a; positive => invertible)"];

(* ---- Phase 8.  chi(z) at background (chain rule) ---- *)
(* dz/d tau = -Hconf (1+z) * a = -Hconf (1+z)  [using a = 1/(1+z)]         *)
(* Substitute 1+z = 1/a so Hconf (1+z) = Hconf/a.                          *)
(* d tau/d lambda = k^0 = -EE/a.                                            *)
(* dz/d lambda = dz/d tau * d tau/d lambda = EE*Hconf*(1+z)/a^2             *)
(* d chi/d z   = (d chi/d lambda) / (dz/d lambda)                            *)

(* Variable names:
     Hconf  = conformal Hubble H := a'/a   (paper's \mathcal{H})
     Hcosm  = cosmic-time Hubble H_c := \dot a/a = \mathcal{H}/a
   The rule below fixes the first definition so every a' in the derivation
   becomes aa * Hconf. Hcosm appears only in the final cosmic-Hubble form. *)
HRules = {Derivative[1][a][tt] -> aa Hconf[tt]};    (* a' = a * Hconf *)
dtaudlamBg = kBg0;                                  (* = -EE/a  *)
dzdtauBg   = -(Derivative[1][a][tt]/aa^2);          (* dz/dtau from 1+z=1/a *)
dzdtauBgH = dzdtauBg /. HRules // Simplify;         (* = -Hconf/a = -Hconf (1+z) *)
dzdlamBg   = dzdtauBgH * dtaudlamBg // Simplify;    (* = EE Hconf / a^2 *)
dchidzBg   = dChiDlamBg / dzdlamBg // Simplify;     (* = a / Hconf = 1 / Hcosm *)

Print["[bg] dz/d tau   = ", dzdtauBgH];
Print["[bg] dz/d lam   = ", dzdlamBg,  "  (expected EE*Hconf/a^2)"];
Print["[bg] d chi/d z (conformal-H form) = ", dchidzBg,  "  (= a/Hconf)"];

(* Re-express in cosmic-Hubble form: Hconf = a * Hcosm, a = 1/(1+z) *)
dchidzBgCosm = dchidzBg /. Hconf[tt] -> aa Hcosm[tt] // FullSimplify;
Print["[bg] d chi/d z (cosmic-H form)    = ", dchidzBgCosm, "  (expected 1/Hcosm(z))"];

(* Sanity: dzdlamBg should equal E0 Hcosm (1+z)^2 when EE = E_0/a, a = 1/(1+z) *)
dzdlamCheck = FullSimplify[(dzdlamBg /. {EE -> E0/aa, Hconf[tt] -> aa Hcosm[tt]}) /. aa -> 1/(1 + zz)];
Print["[bg] dz/d lam in z-only (cosmic-H) form = ", dzdlamCheck,
      "  (expected E0 Hcosm (1+z)^2)"];

(* ---- Phase 9.  First-order null-geodesic RHS (structure) ---- *)
(* Null geodesic:  d k^mu / d lambda = -Gamma^mu_{nu rho} k^nu k^rho.       *)
(* Split Gamma = Gamma_(0) + eps Gamma_(1) and k = k_(0) + eps delta k.    *)
(* At O(eps^0): d k^mu_(0)/d lam = -Gamma^mu_{(0),nu rho} k^nu_(0) k^rho_(0) *)
(* At O(eps^1): d(delta k^mu)/d lam =                                       *)
(*   - Gamma^mu_{(1), nu rho} k^nu_(0) k^rho_(0)                             *)
(*   - 2 Gamma^mu_{(0), nu rho} k^nu_(0) delta k^rho                         *)
(* We emit the FIRST TERM symbolically (the second is routine / involves   *)
(* only background Christoffels).                                           *)

kUp0 = kUp /. eps -> 0;
GamBg = Gam /. eps -> 0 // Simplify;
GamPert = Coefficient[Gam, eps, 1];

(* Background geodesic acceleration (sanity — integrates to standard FLRW) *)
dkBgDlam = Table[
   -Sum[GamBg[[mm, nn, rr]] kUp0[[nn]] kUp0[[rr]], {nn, 4}, {rr, 4}],
   {mm, 4}] // Simplify;
Print["[bg] -Gamma_(0) k k = ", dkBgDlam // Simplify,
      "  (drives the redshift of EE = E_0/a along the ray)"];

(* First-order RHS (from the perturbed Christoffels alone) *)
dkPert1Dlam = Table[
   -Sum[GamPert[[mm, nn, rr]] kUp0[[nn]] kUp0[[rr]], {nn, 4}, {rr, 4}],
   {mm, 4}] // Simplify;

(* Apply Hubble-rule for readability *)
dkPert1DlamH = dkPert1Dlam /. HRules // Simplify;

Print["%%TEX_DELTA_K_RHS_START%%"];
Print[ToString[TeXForm[dkPert1DlamH]]];
Print["%%TEX_DELTA_K_RHS_END%%"];
Print["[O(eps)] -Gamma_(1) k_(0) k_(0)  — line-of-sight integrand"];
Print["         driving the first-order ray deflection and redshift."];
Print["         Time component  -> delta tau / SW-ISW-type integrand;"];
Print["         spatial part    -> delta x^i / transverse gradient of (Phi+Psi)."];

(* ---- Phase 9b.  Scalar-sector source projections ---- *)
(* Specialise to (Psi, Phi) only, set B = 0 and h_{ij} = 0, and split the
   first-order null-geodesic source
       S^mu = -Gamma^mu_{(1),nu rho} k^nu_(0) k^rho_(0)
   into a time component (S^0 -> drives delta tau) and spatial part S^i,
   then project the spatial part onto the radial direction Omega^i and
   the transverse screen triad { x_(0), y_(0) }.  These are the integrands
   of the per-block Jacobian corrections reported in the appendix.          *)

(* Use function-symbol-level substitution so both the bare field and every
   derivative collapse: BF -> (0&) makes Derivative[...][BF][...] = 0 too. *)
scalarSubs = {hF11 -> (0 &), hF22 -> (0 &), hF33 -> (0 &),
              hF12 -> (0 &), hF13 -> (0 &), hF23 -> (0 &),
              BF   -> (0 &)};

dkScalar = (dkPert1DlamH /. scalarSubs) // FullSimplify;

kSrcTime    = dkScalar[[1]] // FullSimplify;
kSrcSpatial = dkScalar[[2 ;; 4]] // FullSimplify;

(* Radial projection: S^i Omega^i   (drives the radial-rate correction)     *)
kSrcRad = Sum[kSrcSpatial[[ii]] OmegaHat[[ii]], {ii, 3}] // FullSimplify;

(* Transverse projection: compute S^i - (S . Omega) Omega^i, the transverse
   part of S^i.  We leave it as a 3-vector with OmegaHat-parameterised
   transverse directions; the appendix presents the invariant form.        *)
kSrcTrans = Table[kSrcSpatial[[ii]] - kSrcRad OmegaHat[[ii]],
                 {ii, 3}] // FullSimplify;

Print["%%TEX_KSRC_TIME_SCALAR_START%%"];
Print[ToString[TeXForm[kSrcTime]]];
Print["%%TEX_KSRC_TIME_SCALAR_END%%"];

Print["%%TEX_KSRC_RAD_SCALAR_START%%"];
Print[ToString[TeXForm[kSrcRad]]];
Print["%%TEX_KSRC_RAD_SCALAR_END%%"];

Print["%%TEX_KSRC_TRANS_SCALAR_START%%"];
Print[ToString[TeXForm[kSrcTrans]]];
Print["%%TEX_KSRC_TRANS_SCALAR_END%%"];

(* As a cross-check, display the transverse source dotted into the flat
   sky-plane basis at theta = 0 (photon along x3): e1 = x_hat_(0),
   e2 = y_hat_(0). The answer should equal -(EE^2/a^2) d_{e_A} (Phi + Psi).  *)
kSrcTransNorth = kSrcTrans /. {th -> 0, ph -> 0} // FullSimplify;
Print["[scalar sector, photon along x3]  S^i_transverse = ", kSrcTransNorth];
Print["  (expected: [-(EE^2/a^2) d_1(Phi+Psi), -(EE^2/a^2) d_2(Phi+Psi), 0])"];

Print[""];
Print["[scalar sector, B = h = 0]"];
Print["  delta tau' RHS   -->   ", kSrcTime];
Print["  delta x^i' radial   --> ", kSrcRad];
Print["  delta x^i' transverse spatial vector = kSrcTrans (full Omega-dependence)."];

(* ---- Phase 10.  Emit TeX snippets for paper ---- *)

Print["%%TEX_DCHIDLAM_BG_START%%"];
Print[ToString[TeXForm[dChiDlamBg]]];
Print["%%TEX_DCHIDLAM_BG_END%%"];

Print["%%TEX_JBG_START%%"];
Print[ToString[TeXForm[JBg]]];
Print["%%TEX_JBG_END%%"];

Print["%%TEX_DETJBG_START%%"];
Print[ToString[TeXForm[detJBg]]];
Print["%%TEX_DETJBG_END%%"];

Print["%%TEX_DCHIDZ_BG_START%%"];
Print[ToString[TeXForm[dchidzBgCosm]]];
Print["%%TEX_DCHIDZ_BG_END%%"];

Print["%%TEX_DZDLAM_BG_START%%"];
Print[ToString[TeXForm[dzdlamCheck]]];
Print["%%TEX_DZDLAM_BG_END%%"];

(* ---- Phase 11.  Archive / regression dump ---- *)
Put[{
   "kSrc_time_scalar"   -> kSrcTime,
   "kSrc_radial_scalar" -> kSrcRad,
   "kSrc_trans_scalar_pole" -> kSrcTransNorth,
   "dchi_dlambda_bg"    -> dChiDlamBg,
   "det_JBg"            -> detJBg,
   "JBg"                -> JBg,
   "dchi_dz_bg"         -> dchidzBg,
   "dzdlamBg_z_only"   -> dzdlamCheck,
   "delta_k_RHS_Gam1"   -> dkPert1DlamH
  },
  "scripts/ray_coord_map.m"];
Print["Saved regression dump to scripts/ray_coord_map.m"];

Print[""];
Print["Done."];
