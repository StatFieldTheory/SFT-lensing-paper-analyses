(* ============================================================
   affine_to_redshift_perturbed.wl
   ------------------------------------------------------------
   Perturbed affine-to-redshift mapping on the Poisson-gauge FLRW
   metric, keeping z = z_bg = 1/a(tau) - 1 as the background redshift.

   We compute the first-order correction to k^0 = d tau / d lambda
   along a past-directed null geodesic moving in the +x3 direction,
   starting from the observer at lambda = 0.  The fractional correction
   to dlambda/dz_bg is then

        delta(d lambda/dz_bg) / (d lambda/dz_bg)_bg  =  - delta k^0 / bar k^0.

   Conventions (matching cosmology.tex Subsection "Driving fields"):
     ds^2 = a^2(tau) [ -(1 + 2 eps Psi) d tau^2
                       - 2 eps (d_i B) d tau dx^i
                       + ((1 - 2 eps Phi) delta_ij + eps h_ij) dx^i dx^j ]
     k^mu = Ef (-u^mu + e^mu)        (past-directed)
     bar k^mu = (-E0/a^2, 0, 0, +E0/a^2) on the background ray
     Observer at (tt = tt0, x^i = 0); a_0 = a[tt0] = 1.

   Strategy:
     1. Build perturbed metric and first-order Christoffels Gamma^0_{munu}
        (direct coordinate computation, same style as driving_fields_poisson.wl).
     2. Null condition 0 = g_{munu} k^mu k^nu at first order:
          fixes an algebraic relation between delta k^0, delta k^3
          and the metric perturbations.
     3. Geodesic equation for k^0 at first order:
          d(delta k^0)/dlambda
             = - delta Gamma^0_{munu} bar k^mu bar k^nu
               - 2 bar Gamma^0_{munu} bar k^mu delta k^nu
        Convert to conformal time via d/dlambda = bar k^0 d/dtau
        (along the ray parameterised by tt alone, since transverse
        positions x^{1,2} are zero at zeroth order and metric perturbations
        along the ray depend only on tt and x^3 = tt0 - tt).
     4. Output  delta k^0 / bar k^0 = dk0Frac  as an explicit expression
        (integral along the ray in general; the integrand is what we emit).
     5. Present both the time-component form and the fractional form.

   Run with:
     ~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh \
        scripts/affine_to_redshift_perturbed.wl 300
   ============================================================ *)

$HistoryLength = 0;

(* ---- Phase 1. Coordinates and fields ---- *)
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

(* Background Christoffels (for the 2 bar Gamma bar k dk term) *)
GamBar = Gam /. eps -> 0;

(* ---- Phase 4. Background photon along x3 axis, past-directed ---- *)
(* bar k^0 = -E0/a^2 (past-directed: d tau / d lambda < 0)            *)
(* bar k^3 = +E0/a^2 (ray moves in +x3 as lambda grows into past)      *)
kBarUp = {-E0/aa^2, 0, 0, E0/aa^2};

(* ---- Phase 5. Perturbed k^mu ---- *)
(* Four unknowns: delta k^mu (mu = 0,1,2,3), functions of (tt, xi).    *)
(* For the past-directed ray with n^i = delta^{i3}, transverse           *)
(* components delta k^{1,2} describe the deflection and do NOT appear   *)
(* in the time-component geodesic equation at leading order (because   *)
(* bar Gamma^0_{0,1}, bar Gamma^0_{0,2} vanish on FLRW background      *)
(* and because bar k^{1,2} = 0). So for delta k^0 we only need         *)
(* delta k^0 and delta k^3 coupled.                                     *)
dk = {dk0[tt], 0, 0, dk3[tt]}; (* transverse = 0 at leading order along the ray *)

(* ---- Phase 6. Null condition at first order ---- *)
(* Build k^mu = bar k^mu + eps dk^mu and expand g_{munu} k^mu k^nu       *)
(* through linear order in eps.                                          *)
kFullUp = kBarUp + eps dk;
nullFull = Sum[gM[[mm, nn]] kFullUp[[mm]] kFullUp[[nn]], {mm, 4}, {nn, 4}];
nullFull = Series[nullFull, {eps, 0, 1}] // Normal;
nullOrder1 = Coefficient[nullFull, eps, 1] // Expand;

(* Solve for dk3[tt] in terms of dk0[tt] + metric perturbations *)
dk3Sol = Solve[nullOrder1 == 0, dk3[tt]][[1, 1, 2]] // Simplify;
Print["[check] dk3 expression length = ", LeafCount[dk3Sol]];

(* ---- Phase 7. Geodesic equation for k^0 at first order ---- *)
(* k^mu nabla_mu k^0 = 0:                                              *)
(*   (k^mu d_mu) k^0 + Gamma^0_{mu nu} k^mu k^nu = 0                    *)
(* At zeroth order this is satisfied automatically.                    *)
(* At first order:                                                      *)
(*   d(delta k^0)/dlambda + delta Gamma^0_{mn} bar k^m bar k^n          *)
(*      + 2 bar Gamma^0_{mn} bar k^m delta k^n = 0                      *)
(* With d/dlambda = bar k^0 d/dtt (at zeroth order, along x3 = x3_0 + k^3 lambda *)
(* but since metric perturbations evaluated ON the ray depend on        *)
(* tt and x3, we need to account for both derivatives. Specifically,    *)
(* on the background ray (past-directed, +x3): dtt/dlambda = bar k^0,    *)
(* dx3/dlambda = bar k^3 = -bar k^0. Hence                               *)
(*   d F(tt,x3)/d lambda = bar k^0 (d F/d tt) + bar k^3 (d F/d x3)      *)
(*                        = bar k^0 (F_{,tt} - F_{,x3}) ≡ bar k^0 DOP[F] *)

(* The user's "D" operator from Subsec 1:  D = d_tau - n^i d_i.         *)
(* Here n^i d_i = d_{x3} (photon along +x3). So:                         *)
DOP[f_] := D[f, tt] - D[f, x3];

(* Compute delta Gamma^0_{mn} bar k^m bar k^n *)
GamPert = Coefficient[Gam, eps, 1]; (* first-order part of Gamma *)

deltaGammaTerm = Sum[GamPert[[1, mm, nn]] kBarUp[[mm]] kBarUp[[nn]],
                     {mm, 4}, {nn, 4}] // Expand;

(* Compute 2 bar Gamma^0_{mn} bar k^m delta k^n *)
(* Recall: bar k^{1,2} = 0; dk^{1,2} = 0 at leading order.                *)
barGammaDk = 2 Sum[GamBar[[1, mm, nn]] kBarUp[[mm]] dk[[nn]],
                   {mm, 4}, {nn, 4}] // Expand;

(* Geodesic equation for k^0:
      d/dlambda (bar k^0 + delta k^0) = -(deltaGammaTerm + barGammaDk)
   where d/dlambda acts on functions of (tt,x3) as bar k^0 DOP[]
   (at zeroth order kinematics), and d(bar k^0)/dlambda satisfies
   the zeroth-order geodesic which we've verified above.              *)

(* The first-order ODE:
      d(delta k^0)/dlambda = - deltaGammaTerm - barGammaDk
   Substitute dk^3 from null condition, substitute d/dlambda:          *)
(* Eliminate dk3 in barGammaDk via nullOrder1 *)
barGammaDkSub = barGammaDk /. dk3[tt] -> dk3Sol;

(* The equation becomes:
      bar k^0 DOP[delta k^0](tt, x3)  =  - RHS *)
(* where RHS = deltaGammaTerm + barGammaDkSub. Along the ray, we          *)
(* treat delta k^0(tt, x3) restricted to (tt, x3_0 + (tt0 - tt)).           *)
(* We can integrate:                                                       *)
(*      delta k^0(tt, x3)_along_ray = - int [RHS / bar k^0] d lambda'.    *)
(* With dlambda = dtt / bar k^0 along the ray,                             *)
(*      = - int RHS / (bar k^0)^2 dtt.                                     *)

RHS = -(deltaGammaTerm + barGammaDkSub) // Expand;

(* Fractional correction to k^0:
   delta k^0 / bar k^0 = (1/bar k^0) int [ RHS / bar k^0 ] d lambda'
                      along the past ray.
   Equivalently, d(delta k^0 / bar k^0)/d lambda
      = (1/bar k^0) d(delta k^0)/d lambda - (d bar k^0/ d lambda)/(bar k^0)^2 * delta k^0.
   But since the zeroth-order geodesic gives a specific d(bar k^0)/d lambda,
   the simplest output is the time-component delta k^0 directly.           *)

Print["[check] RHS leaf count = ", LeafCount[RHS]];

(* ---- Phase 8. Emit the integrand for delta k^0 along the ray ---- *)
(* d(delta k^0)/dlambda = RHS, so                                       *)
(* delta k^0(lambda) = int_0^lambda RHS dlambda' (Born integral)         *)

Print[""];
Print["%%TEX_RHS_full_START%%"];
Print[ToString[TeXForm[RHS]]];
Print["%%TEX_RHS_full_END%%"];

(* Simplify by isolating scalar, vector, tensor sectors *)
RHSscalar = RHS /. {
   BF[tt,x1,x2,x3] -> 0,
   hF11[tt,x1,x2,x3] -> 0, hF22[tt,x1,x2,x3] -> 0, hF33[tt,x1,x2,x3] -> 0,
   hF12[tt,x1,x2,x3] -> 0, hF13[tt,x1,x2,x3] -> 0, hF23[tt,x1,x2,x3] -> 0,
   Derivative[___][BF][__] -> 0,
   Derivative[___][hF11][__] -> 0, Derivative[___][hF22][__] -> 0, Derivative[___][hF33][__] -> 0,
   Derivative[___][hF12][__] -> 0, Derivative[___][hF13][__] -> 0, Derivative[___][hF23][__] -> 0
} // Simplify;

RHSvector = RHS /. {
   PsiF[tt,x1,x2,x3] -> 0,
   PhiF[tt,x1,x2,x3] -> 0,
   hF11[tt,x1,x2,x3] -> 0, hF22[tt,x1,x2,x3] -> 0, hF33[tt,x1,x2,x3] -> 0,
   hF12[tt,x1,x2,x3] -> 0, hF13[tt,x1,x2,x3] -> 0, hF23[tt,x1,x2,x3] -> 0,
   Derivative[___][PsiF][__] -> 0, Derivative[___][PhiF][__] -> 0,
   Derivative[___][hF11][__] -> 0, Derivative[___][hF22][__] -> 0, Derivative[___][hF33][__] -> 0,
   Derivative[___][hF12][__] -> 0, Derivative[___][hF13][__] -> 0, Derivative[___][hF23][__] -> 0
} // Simplify;

RHStensor = RHS /. {
   PsiF[tt,x1,x2,x3] -> 0,
   PhiF[tt,x1,x2,x3] -> 0,
   BF[tt,x1,x2,x3] -> 0,
   Derivative[___][PsiF][__] -> 0, Derivative[___][PhiF][__] -> 0,
   Derivative[___][BF][__] -> 0
} // Simplify;

Print["%%TEX_RHS_scalar_START%%"];
Print[ToString[TeXForm[RHSscalar]]];
Print["%%TEX_RHS_scalar_END%%"];

Print["%%TEX_RHS_vector_START%%"];
Print[ToString[TeXForm[RHSvector]]];
Print["%%TEX_RHS_vector_END%%"];

Print["%%TEX_RHS_tensor_START%%"];
Print[ToString[TeXForm[RHStensor]]];
Print["%%TEX_RHS_tensor_END%%"];

(* ---- Phase 9. Closed-form delta k^0 / bar k^0 for scalar sector ---- *)
(* On the ray, we evaluate scalars at x1=x2=0 (ray passes through origin)  *)
(* and note that the past-directed ray at x3 = tt0 - tt satisfies           *)
(*   d(Phi[tt, x3]) / d lambda = bar k^0 DOP[Phi]                         *)
(* So the Born-limit integrated correction is                               *)
(*   delta k^0(lambda) = int_0^lambda RHS_scalar d lambda'                 *)
(*                    = -int_0^tt (RHS_scalar / bar k^0) dtt'  (past-directed) *)

(* Express RHS_scalar /( bar k^0 ) symbolically; this is the integrand       *)
(* in terms of conformal time.                                              *)
scalarIntegrandConfTime = RHSscalar / kBarUp[[1]] // Simplify;
Print["%%TEX_scalar_integrand_confTime_START%%"];
Print[ToString[TeXForm[scalarIntegrandConfTime]]];
Print["%%TEX_scalar_integrand_confTime_END%%"];

(* ---- Phase 10. Input-form archive ---- *)
Print[""];
Print["%%INPUT_ARCHIVE_START%%"];
Print["nullOrder1   = ", ToString[InputForm[nullOrder1]]];
Print["dk3Sol        = ", ToString[InputForm[dk3Sol]]];
Print["deltaGammaTerm= ", ToString[InputForm[deltaGammaTerm]]];
Print["barGammaDk    = ", ToString[InputForm[barGammaDkSub]]];
Print["RHS           = ", ToString[InputForm[RHS]]];
Print["RHSscalar     = ", ToString[InputForm[RHSscalar]]];
Print["RHSvector     = ", ToString[InputForm[RHSvector]]];
Print["RHStensor     = ", ToString[InputForm[RHStensor]]];
Print["%%INPUT_ARCHIVE_END%%"];

Print[""];
Print["Done."];
