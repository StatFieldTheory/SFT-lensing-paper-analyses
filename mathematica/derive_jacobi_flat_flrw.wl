(* ============================================================
   derive_jacobi_flat_flrw.wl
   ------------------------------------------------------------
   Symbolic verification that D(lambda) := a(lambda) chi(lambda)
   solves the Jacobi geodesic-deviation equation
       D''(lambda) = bar Phi_00(lambda) * D(lambda),
       D(0) = 0,  D'(0) = 1
   on a flat FLRW background.

   Inputs we take as GIVEN (cited; not re-derived here):

     (P1) FLRW metric is conformal-time,
              ds^2 = a^2(eta) [-d eta^2 + dx_i dx^i].
          (sections/cosmology.tex Eq. ``poisson gauge'' with
          perturbations zeroed.)

     (P2) Past-directed null geodesic, E_0 = 1.  Then on the radial line
              d eta / d lambda = -1 / a^2,   d chi / d lambda = +1 / a^2.
          (sections/cosmology.tex Eq. ``dchi dlambda bg'' and
          scripts/ray_coord_map.wl Phases 6, 8.)

     (P3) bar Phi_00 in the FLRW limit:
              bar Phi_00(lambda) = -(calH^2(eta) - calH'(eta)) / a^4(eta),
          with calH := a' / a (prime = d/d eta).
          (sections/cosmology.tex Eq. ``driving bg'' and
          sections/sachs_dynamics.tex Eq. ``Phi 00 from Ricci'';
          numerical implementation in
          scripts/ccl/driving_bg.calH2_minus_calHprime_lcdm.)

     (P4) Friedmann pair on flat FLRW with matter, radiation, and Lambda:
              3 calH^2          = 8 pi G a^2 (rho_m + rho_r + rho_L),
              calH^2 - calH'    = 4 pi G a^2 (rho_m + (4/3) rho_r).
          (Standard textbook; consistent with cosmology.tex.)

   What this script PROVES symbolically:

     (T1) D = a chi satisfies D'' = bar Phi_00 D, derived only from
          the chain rule (P1)+(P2) and the form of bar Phi_00 (P3).
          The proof is independent of the matter content; it works
          for any flat-FLRW a(eta).

     (T2) D(0) = 0 and D'(0) = 1, using a(lambda=0) = 1,
          chi(lambda=0) = 0, and d chi/d lambda(lambda=0) = 1
          (E_0 = 1 with a=1 at the observer).

     (T3) The Friedmann pair (P4) implies the fractional identity
              calH^2 - calH' = calH^2 [(3/2) Omega_m + 2 Omega_r],
          which is what scripts/ccl/driving_bg.calH2_minus_calHprime_lcdm
          implements.

   Run with:
       wolframscript -file scripts/derive_jacobi_flat_flrw.wl
   ============================================================ *)

ClearAll["Global`*"];

Print["======================================================="];
Print["  Jacobi flat-FLRW: verify D(lambda) = a(lambda) chi(lambda) "];
Print["======================================================="];

(* --------------------------------------------------------------
   Phase 1.  Chain rule from eta to lambda.
   -------------------------------------------------------------- *)

calH[eta_]      := a'[eta] / a[eta];
calHprime[eta_] := D[calH[etaIn], etaIn] /. etaIn -> eta;

(* Substitution rules to express derivatives of a in terms of (calH,
   calH').  These come from differentiating calH = a'/a once. *)
substAprime  = a'[eta]  -> a[eta] calH[eta];
substA2prime = a''[eta] -> a[eta] (calH[eta]^2 + calHprime[eta]);

(* From (P2): d eta/d lambda = -1/a^2,  d chi/d lambda = +1/a^2. *)
detaDlam = -1 / a[eta]^2;
dchiDlam =  1 / a[eta]^2;

(* da/dlambda = a'(eta) * d eta/d lambda = -a'(eta)/a^2 = -calH/a. *)
daDlam = a'[eta] * detaDlam // Simplify;
daDlamHform = daDlam /. substAprime // Simplify;

(* d2 a / d lambda^2 = d/d lambda (da/d lambda)
                     = (d eta/d lambda) * d/d eta (da/d lambda). *)
d2aDlam = D[a'[eta] * detaDlam, eta] * detaDlam // Simplify;
d2aDlamHform = d2aDlam /. substA2prime /. substAprime // Simplify;

(* d2 chi / d lambda^2 = d/d lambda (1/a^2) = -2 / a^3 * (da/d lambda)
                       = -2/a^3 * (-calH/a) = 2 calH / a^4. *)
d2chiDlam = D[dchiDlam, eta] * detaDlam /. substAprime // Simplify;

Print["Phase 1:  da/dlambda      = ", daDlamHform];
Print["          d2a/dlambda^2   = ", d2aDlamHform];
Print["          d chi/dlambda   = ", dchiDlam];
Print["          d2chi/dlambda^2 = ", d2chiDlam];

(* --------------------------------------------------------------
   Phase 2.  Compute D'(lambda) and D''(lambda) for D := a chi.
   -------------------------------------------------------------- *)

DP  = daDlamHform * chi + a[eta] * dchiDlam // Simplify;
DPP = d2aDlamHform * chi
    + 2 daDlamHform * dchiDlam
    + a[eta] * d2chiDlam // Simplify;

Print["\nPhase 2:  D'(lambda)  = ", DP];
Print["          D''(lambda) = ", DPP];

(* --------------------------------------------------------------
   Phase 3.  Compare D'' to bar Phi_00 D.  bar Phi_00 from (P3).
   -------------------------------------------------------------- *)

barPhi00 = -(calH[eta]^2 - calHprime[eta]) / a[eta]^4;
RHS = barPhi00 * (a[eta] chi) // Simplify;

T1Residual = Simplify[DPP - RHS];

Print["\nPhase 3:  bar Phi_00 D = ", RHS];
Print["          T1 residual D'' - bar Phi_00 D = ", T1Residual];
Print["          (PASS iff residual is 0)"];

(* --------------------------------------------------------------
   Phase 4.  Initial conditions D(0) = 0, D'(0) = 1.
   At the observer (lambda = 0): a = 1, chi = 0, EE = 1 (E_0 = 1).
   So d chi/d lambda |_{lam=0} = 1/a^2 = 1.
   -------------------------------------------------------------- *)

ICRule = {a[eta] -> 1, chi -> 0};
DAtZero  = (a[eta] chi) /. ICRule;
DPAtZero = DP /. ICRule;

Print["\nPhase 4:  D(0)  = a(0) chi(0) = ", DAtZero,  "   (expected 0)"];
Print["          D'(0) = ", DPAtZero, "   (expected 1)"];

(* --------------------------------------------------------------
   Phase 5.  Friedmann fractional form (T3).
   This is the algebraic identity used by
   scripts/ccl/driving_bg.calH2_minus_calHprime_lcdm.
   -------------------------------------------------------------- *)

Print["\nPhase 5:  Friedmann pair => fractional form  ..."];

eq1 = 3 cH2 == 8 Pi G a[eta]^2 (rhoM + rhoR + rhoL);
eq2 = cH2 - cHp == 4 Pi G a[eta]^2 (rhoM + (4/3) rhoR);
sol = Solve[{eq1, eq2}, {cH2, cHp}][[1]];

cH2val = cH2 /. sol;
cHpval = cHp /. sol;

(* rho_crit(a) = 3 calH^2 / (8 pi G a^2). *)
rhoCrit = 3 cH2val / (8 Pi G a[eta]^2) // Simplify;
OmegaM = rhoM / rhoCrit // Simplify;
OmegaR = rhoR / rhoCrit // Simplify;

cH2MinusCHp = cH2val - cHpval // Simplify;
expectedFrac = cH2val * ((3/2) OmegaM + 2 OmegaR) // Simplify;

T3Residual = Simplify[cH2MinusCHp - expectedFrac];

Print["          calH^2 - calH'                              = ", cH2MinusCHp];
Print["          calH^2 [(3/2) Omega_m + 2 Omega_r]          = ", expectedFrac];
Print["          T3 residual (LHS - RHS)                     = ", T3Residual];
Print["          (PASS iff residual is 0)"];

(* --------------------------------------------------------------
   Phase 6.  theta^(sa) integral identity and R-propagator (T4, T5).
   These are the identities the Python pipeline uses to bypass any
   numerical integration of theta^(sa); they are the user-facing
   reason we never integrate theta^(sa) and instead read D from
   theta_sa.npz.

     theta^(sa)(lambda) := D'(lambda) / D(lambda)
       so theta^(sa) is the lambda-derivative of ln D(lambda).

   (T4) Integral identity
     int_{lambda'}^{lambda} theta^(sa)(tau) d tau
        = int_{lambda'}^{lambda} (d/d tau) ln D(tau) d tau
        = ln D(lambda) - ln D(lambda').

   (T5) Response propagator
     R(lambda, lambda') = exp(-2 int_{lambda'}^{lambda} theta^(sa) d tau)
                        = exp(-2 (ln D(lambda) - ln D(lambda')))
                        = (D(lambda') / D(lambda))^2.

   sections/path_int.tex Eq. ``explicit resp op'' is the source for R.
   Mathematica check is straightforward: differentiate ln D(tau) and
   integrate it back.
   -------------------------------------------------------------- *)

Print["\nPhase 6:  theta^(sa) integral and R-propagator  ..."];

(* T4 via the differential form: verify theta^(sa)(tau) = d/dtau ln D(tau).
   Integration from lp to l then follows by the fundamental theorem of
   calculus, with no further symbolic work. *)

DTest[tau_]      := dFun[tau];
thetaSaTest[t_]  := DTest'[t] / DTest[t];

T4Diff = thetaSaTest[tau] - D[Log[DTest[tau]], tau] // Simplify;
Print["  T4: theta^(sa)(tau) - d/dtau ln D(tau) = ", T4Diff,
      "    (PASS iff 0)"];
Print["      => by FTC,  int_{lp}^{l} theta^(sa) d tau"];
Print["                   = ln D(l) - ln D(lp)."];
T4Residual = T4Diff;

(* T5: R(l, lp) = exp(-2 integral) = (D(lp)/D(l))^2.
   Take the FTC result as given (T4 above), then simplify
   exp(-2 (ln D(l) - ln D(lp))) directly. *)
RFromIntegral = Exp[-2 (Log[DTest[l]] - Log[DTest[lp]])];
RRatioForm    = (DTest[lp] / DTest[l])^2;
T5Residual = FullSimplify[RFromIntegral - RRatioForm,
   Assumptions -> {DTest[l] > 0, DTest[lp] > 0}];
Print["  T5: exp(-2 (ln D(l) - ln D(lp))) - (D(lp)/D(l))^2 = ",
      T5Residual, "    (PASS iff 0)"];

(* --------------------------------------------------------------
   Final summary.
   -------------------------------------------------------------- *)

t1Pass = (T1Residual === 0);
t2Pass = (DAtZero === 0) && (DPAtZero === 1);
t3Pass = (T3Residual === 0);
t4Pass = (T4Residual === 0);
t5Pass = (T5Residual === 0);

Print["\n======================================================="];
Print["Summary"];
Print["  T1: D'' = bar Phi_00 D            -> ", If[t1Pass, "PASS", "FAIL"]];
Print["  T2: IC D(0)=0, D'(0)=1            -> ", If[t2Pass, "PASS", "FAIL"]];
Print["  T3: calH^2 - calH' fractional ID  -> ", If[t3Pass, "PASS", "FAIL"]];
Print["  T4: int theta^(sa) = ln D(l)-ln D(lp) -> ", If[t4Pass, "PASS", "FAIL"]];
Print["  T5: R = (D(lp)/D(l))^2            -> ", If[t5Pass, "PASS", "FAIL"]];
If[t1Pass && t2Pass && t3Pass && t4Pass && t5Pass,
   Print["\n[PASS]  All five identities verified."];
,
   Print["\n[FAIL]  Some assertion did not collapse to zero; inspect above."];
];
Print["======================================================="];
