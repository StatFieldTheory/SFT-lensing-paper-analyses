(* ============================================================
   affine_to_redshift_general.wl
   ------------------------------------------------------------
   Symbolic derivation of  dE/dlambda = k^mu k^nu nabla_nu u_mu
   and  dz/dlambda = (1/E0) dE/dlambda
   on the perturbed Poisson-gauge FLRW metric.

   Conventions (matching cosmology.tex Subsection "From Affine Parameter
   to Redshift" and Subsection "Driving fields"):

     ds^2 = a^2(tau) [ -(1 + 2 eps Psi) d tau^2
                       - 2 eps (d_i B) d tau dx^i
                       + ((1 - 2 eps Phi) delta_ij + eps h_ij) dx^i dx^j ]

     u^mu = ( 1/(a Sqrt[1 + 2 eps Psi]),  0 )       (Hubble-flow observer)

     e^mu = spatial direction of the past-directed photon (unit spacelike,
            orthogonal to u^mu).  In the driving-fields convention this
            is the third tetrad vector along the line of sight:
                e^0 = -eps (d_3 B) / a
                e^i = (1/a) [ delta^{i3}
                              + eps ( Phi delta^{i3} - (1/2) h^{i3} ) ]

     k^mu = EE * ( - u^mu + e^mu )                   (past-directed)
     E  = EE := local photon energy measured by the comoving observer
          the ray has just passed; kept symbolic throughout.
          We do NOT substitute E = E0 / a anywhere: that is a
          background relation, not a definition.

     E0 = E(lambda = 0): constant photon energy measured by today's
          observer; appears only in  1+z = E/E0.

   The calculation proceeds to first order in the perturbation bookkeeper
   eps.  Only u_mu and the Christoffels enter nabla_nu u_mu with first-
   order corrections; k^mu itself is kept at zeroth order in the eps
   expansion (Born-like treatment of the four-momentum).  The coupled
   first-order correction to k^mu is computed separately in
   scripts/affine_to_redshift_perturbed.wl.

   Run with:
     ~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh \
        scripts/affine_to_redshift_general.wl 300
   ============================================================ *)

$HistoryLength = 0;

(* ---- Phase 1. Coordinates and perturbation fields ---- *)
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

linearize[expr_] := Normal @ Series[expr, {eps, 0, 1}];

(* ---- Phase 3. Inverse metric and Christoffels ---- *)
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

(* ---- Phase 4. Observer four-velocity u^mu, u_mu ---- *)
uUp = {1/(aa Sqrt[1 + 2 eps Psifld]), 0, 0, 0};
uUp = linearize[uUp];

uLo = Table[Sum[gM[[mm, nn]] uUp[[nn]], {nn, 4}], {mm, 4}];
uLo = linearize[uLo];

uNorm = Sum[gM[[mm, nn]] uUp[[mm]] uUp[[nn]], {mm, 4}, {nn, 4}];
uNorm = linearize[uNorm];
Print["[check] u.u (should be -1) = ", Simplify @ uNorm];

(* ---- Phase 5. Line-of-sight tetrad vector e^mu ---- *)
(* Following driving_fields_poisson.wl convention (past-directed photon   *)
(* along +x3).  e^mu is spatial (u.e = 0), unit spacelike (e.e = 1),     *)
(* expressed as the third spatial tetrad vector.                          *)
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
nUp = tet[3];

(* Sanity checks on the tetrad. *)
unDot = Sum[gM[[mm, nn]] uUp[[mm]] nUp[[nn]], {mm, 4}, {nn, 4}];
unDot = linearize[unDot];
Print["[check] u.n (should be 0) = ", Simplify @ unDot];

nnDot = Sum[gM[[mm, nn]] nUp[[mm]] nUp[[nn]], {mm, 4}, {nn, 4}];
nnDot = linearize[nnDot];
Print["[check] n.n (should be 1) = ", Simplify @ nnDot];

(* ---- Phase 6. Photon four-momentum k^mu = EE ( -u^mu + e^mu ) ---- *)
(* EE is the LOCAL photon energy measured by the comoving observer       *)
(* the ray has just passed.  It is kept as a free symbol: no relation    *)
(* EE = E0/a is imposed.                                                  *)
kUp = EE (-uUp + nUp);
kUp = linearize[kUp];

(* Null check *)
kNorm = Sum[gM[[mm, nn]] kUp[[mm]] kUp[[nn]], {mm, 4}, {nn, 4}];
kNorm = linearize[kNorm];
Print["[check] k.k (should be 0) = ", Simplify @ kNorm];

(* u.k at observer (lambda=0, where perturbations also hold) *)
ukDot = Sum[gM[[mm, nn]] uUp[[mm]] kUp[[nn]], {mm, 4}, {nn, 4}];
ukDot = linearize[ukDot];
Print["[check] u.k (should be +EE under past-directed E = u_mu k^mu) = ",
      Simplify @ ukDot];

(* ---- Phase 7. Compute k^mu k^nu nabla_nu u_mu ---- *)
nabluLowered = Table[
   D[uLo[[mm]], coords[[nn]]] - Sum[Gam[[lam, nn, mm]] uLo[[lam]], {lam, 4}],
   {mm, 4}, {nn, 4}];
nabluLowered = linearize[nabluLowered];

dEdlam = Sum[kUp[[mm]] kUp[[nn]] nabluLowered[[mm, nn]], {mm, 4}, {nn, 4}];
dEdlam = linearize[dEdlam] // Expand;

(* Split into background (eps^0) and first-order (eps^1) *)
dEdlamBg   = dEdlam /. eps -> 0                 // Simplify;
dEdlamPert = Coefficient[dEdlam, eps, 1]        // Simplify;

(* Substitute a' -> a * Hcf (conformal Hubble) for clarity *)
HubbleRules = {Derivative[1][a][tt] -> aa Hcf[tt]};
dEdlamBgH   = dEdlamBg   /. HubbleRules // Simplify;
dEdlamPertH = dEdlamPert /. HubbleRules // Simplify;

Print["%%TEX_dEdlamBg_START%%"];
Print[ToString[TeXForm[dEdlamBgH]]];
Print["%%TEX_dEdlamBg_END%%"];

Print["%%TEX_dEdlamPert_START%%"];
Print[ToString[TeXForm[dEdlamPertH]]];
Print["%%TEX_dEdlamPert_END%%"];

(* dz/dlambda = dE/dlambda / E0 *)
dzdlamBgH   = dEdlamBgH  / E0 // Simplify;
dzdlamPertH = dEdlamPertH / E0 // Simplify;

Print["%%TEX_dzdlamBg_START%%"];
Print[ToString[TeXForm[dzdlamBgH]]];
Print["%%TEX_dzdlamBg_END%%"];

Print["%%TEX_dzdlamPert_START%%"];
Print[ToString[TeXForm[dzdlamPertH]]];
Print["%%TEX_dzdlamPert_END%%"];

(* ---- Phase 8. Cross-check in cosmic-time coordinates ---- *)
(* As an independent sanity check, redo the background calculation       *)
(* in cosmic time  ds^2 = -dt^2 + a^2 delta_ij dx^i dx^j  to confirm     *)
(* that the resulting (1+z)-power agrees with the conformal result       *)
(* via the identity  H_conformal = a * H_cosmic.                          *)
Print["[progress] Cosmic-time cross-check ..."];
Module[{gMct, g0ctInv, dgMct, GamCt, ut, ul, nct, kct, nab, ratect, nullct, dE},
  (* cosmic-time background metric *)
  gMct = DiagonalMatrix[{-1, aa^2, aa^2, aa^2}];
  g0ctInv = Inverse[gMct];
  dgMct = Table[D[gMct[[bb, nn]], coords[[mm]]], {mm, 4}, {bb, 4}, {nn, 4}];
  GamCt = Table[
     (1/2) Sum[g0ctInv[[aaI, bb]] (dgMct[[mm, bb, nn]] + dgMct[[nn, bb, mm]] - dgMct[[bb, mm, nn]]),
               {bb, 4}],
     {aaI, 4}, {mm, 4}, {nn, 4}];
  ut = {1, 0, 0, 0};                       (* u^mu in cosmic time *)
  ul = {-1, 0, 0, 0};                      (* u_mu *)
  nct = {0, 0, 0, 1/aa};                   (* n^mu, unit spacelike along x3 *)
  kct = EE (-ut + nct);                    (* k^mu *)
  nullct = Sum[gMct[[mm, nn]] kct[[mm]] kct[[nn]], {mm, 4}, {nn, 4}];
  Print["[check cosmic-time] k.k (should be 0) = ", Simplify @ nullct];
  nab = Table[
     D[ul[[mm]], coords[[nn]]] - Sum[GamCt[[lam, nn, mm]] ul[[lam]], {lam, 4}],
     {mm, 4}, {nn, 4}];
  dE = Sum[kct[[mm]] kct[[nn]] nab[[mm, nn]], {mm, 4}, {nn, 4}] // Simplify;
  ratect = dE /. {Derivative[1][a][tt] -> aa Hcosm[tt]} // Simplify;
  Print["[check cosmic-time] dE/dlambda = ", Simplify @ ratect];
  Print["  (expected: EE^2 * Hcosm[tt], independent of a)"];
];

(* ---- Phase 9. Input-form archive ---- *)
Print[""];
Print["%%INPUT_ARCHIVE_START%%"];
Print["dE/dlambda (bg)   = ", ToString[InputForm[dEdlamBgH]]];
Print["dE/dlambda (pert) = ", ToString[InputForm[dEdlamPertH]]];
Print["dz/dlambda (bg)   = ", ToString[InputForm[dzdlamBgH]]];
Print["dz/dlambda (pert) = ", ToString[InputForm[dzdlamPertH]]];
Print["%%INPUT_ARCHIVE_END%%"];

Print[""];
Print["Done."];
