(* ::Package:: *)

(* Machine verification of the redshift-dependence statements in
   notes/theory_redshift_dependence.md.

   Run:  wolframscript -file verify_redshift_scalings.wl

   Five check groups, each printing PASS/FAIL:

   V1  power-law soft factor S_inf and its consistency with the
       projection-slice identity (3D xi -> transverse projection);
   V2  power-law hard factor H;
   V3  per-shell tree scaling exponents in (1+z), chi, gamma, assembled
       from the source-read prefactors (spec section 2), certified here;
   V4  low-redshift fold: closed-form window, the exact O0 and FK fold
       integrals for power-law shells, the chi_s exponents of both, and
       the gamma-independence of their ratio;
   V5  numeric cross-check of the V4 closed forms against NIntegrate.

   Conventions: P(k) = A k^n at z = 0; all lengths in one unit; gamma the
   flat-sky external angle; chi the shell distance; chi_s the source.
   The fold expressions transcribe fk_analytic.py eq. (FK) and
   driver_stats.order0_mc verbatim (window G, response u^4 factors). *)

report[label_, ok_] := Print[If[TrueQ[ok], "PASS  ", "FAIL  "], label];

(* ---------------------------------------------------------------- V1 *)
(* Soft factor for P = A k^n at shell chi:
   S_inf(gamma; chi) = INT_0^inf dv v J0(v gamma) A (v/chi)^n           *)

sInf = Assuming[gamma > 0 && chi > 0 && amp > 0 && -2 < n < -1/2,
   Integrate[v BesselJ[0, v gamma] amp (v/chi)^n, {v, 0, Infinity}]];

cS = Simplify[sInf/(amp chi^-n gamma^(-(n + 2)))];
report["V1a  S_inf = A chi^-n gamma^-(n+2) * c_S(n),  c_S = " <>
   ToString[InputForm[cS]],
  FreeQ[cS, gamma] && FreeQ[cS, chi] && FreeQ[cS, amp]];

(* c_S against the known Hankel closed form 2^(1+n) Gamma(1+n/2)/Gamma(-n/2);
   the two differ by the recurrence Gamma(1-n/2) = (-n/2) Gamma(-n/2), which
   plain Simplify does not apply, so FullSimplify the ratio instead. *)
report["V1b  c_S matches 2^(1+n) Gamma(1+n/2)/Gamma(-n/2)",
  FullSimplify[cS/(2^(1 + n) Gamma[1 + n/2]/Gamma[-n/2]),
    -2 < n < -1/2] === 1];

(* Projection-slice consistency:
   xi3D(r) = 1/(2 pi^2) INT dk k^2 P(k) sin(kr)/(kr);
   Sigma(rp) = INT_-inf^inf dDelta xi3D(sqrt(rp^2 + Delta^2));
   identity:  S_inf = 2 pi chi^2 Sigma(chi gamma).                       *)

xi3D = Assuming[r > 0 && amp > 0 && -3 < n < -1,
   Integrate[amp k^(2 + n) Sin[k r]/(k r), {k, 0, Infinity}]/(2 Pi^2)];

sigmaProj = Assuming[rp > 0 && amp > 0 && -2 < n < -1,
   Integrate[xi3D /. r -> Sqrt[rp^2 + Delta^2], {Delta, -Infinity,
     Infinity}]];

lhsOverRhs = Assuming[gamma > 0 && chi > 0 && amp > 0 && -2 < n < -1,
   FullSimplify[sInf/(2 Pi chi^2 (sigmaProj /. rp -> chi gamma))]];
report["V1c  projection-slice identity S_inf = 2 pi chi^2 Sigma(chi gamma)",
  lhsOverRhs === 1];

(* ---------------------------------------------------------------- V2 *)
(* Hard factor for P = A k^n:
   H(chi; L) = INT_lc^L du u (17/7 - n/2) A (u/chi)^n                    *)

hHard = Assuming[el > lc > 0 && chi > 0 && amp > 0 && -2 < n < 0,
   Integrate[u (17/7 - n/2) amp (u/chi)^n, {u, lc, el}]];
hPred = amp chi^-n (17/7 - n/2) (el^(2 + n) - lc^(2 + n))/(2 + n);
report["V2   H = A chi^-n (17/7 - n/2)(L^(2+n) - lc^(2+n))/(2+n)",
  Simplify[hHard - hPred] === 0];

(* ---------------------------------------------------------------- V3 *)
(* Per-shell tree scaling exponents. Inputs read from the source
   (input_ready_b_to_zeta_spec.md section 2, KAPPA3:3537-3626):
     per Phi00 leg:  A(a) (1+z)^4  with A(a) = -(3/2) Om H0^2 (1+z)
     three legs:     R_TTT = [A(a) (1+z)^4]^3
     radial measure: jac_lambda = (d lambda / d chi)^2 = (1+z)^-4
     geometry:       chi^-4
     growth:         D(z)^4  (B evaluated at z=0, then growth^4)
   plus the factorized z=0 quadrature (2 pi)^2 S_inf H from V1/V2.       *)

perLeg = -(3/2) om h0^2 (1 + z) (1 + z)^4;
prefShell = perLeg^3 (1 + z)^-4 chi^-4 growthD^4;
zetaShell = prefShell (2 Pi)^2 (amp chi^-n gamma^(-(n + 2)) cS0) *
    (amp chi^-n cH0 el^(2 + n));

report["V3a  explicit (1+z) exponent (before growth) = 11",
  Exponent[Simplify[zetaShell /. growthD -> 1], 1 + z] === 11];
report["V3b  chi exponent = -4 - 2n",
  Simplify[Exponent[zetaShell, chi] - (-4 - 2 n)] === 0];
report["V3c  gamma exponent = -(2 + n)",
  Simplify[Exponent[zetaShell, gamma] - (-(2 + n))] === 0];
report["V3d  growth exponent = 4 (tree: B(z) = D^4 B(0), scale-free)",
  Exponent[zetaShell, growthD] === 4];

(* ---------------------------------------------------------------- V4 *)
(* Low-redshift fold. Background: D_J = a chi -> lambda; u(l) = l/ls.
   Window (driver_stats eq. 3):
     G(l) = INT_l^ls (D_J(l)/D_J(t))^2 dt -> l^2 (1/l - 1/ls)
   FK fold (fk_analytic eq. FK):
     xi_FK = INT_0^ls dl G(l) u(l)^4 Hresp(l) Z(l; gamma),
     Hresp(l) = INT_l^ls G(t) u(t)^-4 dt
   O0 fold (driver_stats order0_mc):
     xi_O0 = INT_0^ls dl G(l)^2 C2(l; gamma)
   Power-law shells (from V1-V3 at z ~ 0):
     Z(l)  = kz l^(-4 - 2 n) gamma^(-(n + 2))     (three-point vertex)
     C2(l) = k2 l^(-2 - n)   gamma^(-(n + 2))     (two-point density)   *)

gWin = Assuming[ls > l > 0, Integrate[(l/t)^2, {t, l, ls}]];
report["V4a  window G(l) = l (1 - l/ls)",
  Simplify[gWin - l (1 - l/ls)] === 0];

uOf[l_] := l/ls;

hResp = Assuming[ls > l > 0,
   Integrate[(t (1 - t/ls)) uOf[t]^-4, {t, l, ls}]];

(* Convergence domains: the O0 integrand ~ l^{-n} at the origin needs
   n < 1; the FK integrand reduces to l^{-1-2n}(ls-l)^3/(2 ls^3), so it
   needs n < 0. Both closed forms therefore hold on -2 < n < 0, which
   covers BOTH sides of the sign-flip point n = -1 of the ratio exponent;
   the soft-factor closed form of V1 (inherited by the anchor's
   interpretation) is the binding constraint -2 < n < -1/2. *)
xiO0 = Assuming[ls > 0 && gamma > 0 && -2 < n < 0,
   Integrate[(l (1 - l/ls))^2 k2 l^(-2 - n) gamma^(-(n + 2)),
    {l, 0, ls}]];

xiFK = Assuming[ls > 0 && gamma > 0 && -2 < n < 0,
   Integrate[(l (1 - l/ls)) uOf[l]^4 hResp kz l^(-4 - 2 n) *
     gamma^(-(n + 2)), {l, 0, ls}]];

expO0 = Exponent[xiO0, ls];
expFK = Exponent[xiFK, ls];
ratioFKO0 = FullSimplify[xiFK/xiO0];

Print["      xi_O0  = ", InputForm[Simplify[xiO0]]];
Print["      xi_FK  = ", InputForm[Simplify[xiFK]]];
Print["      chi_s exponents:  O0 -> ", InputForm[expO0],
  "   FK -> ", InputForm[expFK],
  "   ratio -> ", InputForm[Simplify[expFK - expO0]]];

report["V4b  both folds converge for -2 < n < 0 (closed forms exist)",
  FreeQ[xiO0, Integrate] && FreeQ[xiFK, Integrate]];
report["V4c  FK/O0 fold ratio is gamma-independent (power law)",
  FreeQ[Assuming[gamma > 0 && ls > 0,
     FullSimplify[PowerExpand[ratioFKO0]]], gamma]];
report["V4d  FK/O0 chi_s exponent = -(n + 1)  [derived, not assumed]",
  Simplify[expFK - expO0 - (-(n + 1))] === 0];

(* ---------------------------------------------------------------- V5 *)
(* Numeric cross-check of the closed forms on BOTH sides of n = -1
   (the sign-flip point of the ratio exponent): n = -3/5 and n = -3/2. *)

nNum = -3/5; lsNum = 1;
o0Closed = xiO0 /. {n -> nNum, ls -> lsNum, k2 -> 1, gamma -> 1};
fkClosed = xiFK /. {n -> nNum, ls -> lsNum, kz -> 1, gamma -> 1};
o0Num = NIntegrate[(l (1 - l))^2 l^(-2 - nNum), {l, 0, 1},
   WorkingPrecision -> 20];
hRespN[l_?NumericQ] := NIntegrate[t (1 - t) t^-4, {t, l, 1},
   WorkingPrecision -> 20];
fkNum = NIntegrate[l (1 - l) l^4 hRespN[l] l^(-4 - 2 nNum), {l, 0, 1},
   WorkingPrecision -> 15, PrecisionGoal -> 8];
report["V5a  O0 closed form vs NIntegrate  (rel err < 1e-8)",
  Abs[o0Closed/o0Num - 1] < 10^-8];
report["V5b  FK closed form vs NIntegrate  (rel err < 1e-6)",
  Abs[fkClosed/fkNum - 1] < 10^-6];

(* n = -3/2: the n < -1 side. Exact rational values of the closed forms
   (independent sympy audit 2026-08-26: xi_FK = 1/120, xi_O0 = 16/315
   at ls = 1, gamma = 1, unit prefactors) checked against NIntegrate.  *)
nN2 = -3/2;
o0C2 = xiO0 /. {n -> nN2, ls -> 1, k2 -> 1, gamma -> 1};
fkC2 = xiFK /. {n -> nN2, ls -> 1, kz -> 1, gamma -> 1};
o0N2 = NIntegrate[(l (1 - l))^2 l^(-2 - nN2), {l, 0, 1},
   WorkingPrecision -> 20];
fkN2 = NIntegrate[l (1 - l) l^4 hRespN[l] l^(-4 - 2 nN2), {l, 0, 1},
   WorkingPrecision -> 15, PrecisionGoal -> 8];
report["V5c  n = -3/2 (n < -1 side): closed forms = {16/315, 1/120} " <>
   "and match NIntegrate",
  Simplify[o0C2 - 16/315] === 0 && Simplify[fkC2 - 1/120] === 0 &&
   Abs[o0C2/o0N2 - 1] < 10^-8 && Abs[fkC2/fkN2 - 1] < 10^-6];

Print["\nDone."];
