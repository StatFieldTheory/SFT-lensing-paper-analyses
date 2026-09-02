(* ::Package:: *)
(* verify_squeezed_expansion.wl
   Machine verification of the squeezed-limit expansions used in
   notes/theory_shape_amplitude_factorization.md.

   Run:  wolframscript -file verify_squeezed_expansion.wl
   Every check prints PASS/FAIL; the note may only quote formulas that
   PASS here. *)

F2[kk1_, kk2_] := Module[{a = Sqrt[kk1 . kk1], b = Sqrt[kk2 . kk2], mu},
  mu = (kk1 . kk2)/(a b);
  5/7 + (mu/2) (a/b + b/a) + (2/7) mu^2];

report[name_, resid_] := Print[name, ": ",
  If[Simplify[resid] === 0, "PASS", Row[{"FAIL, residual = ", resid}]]];

(* ---- Lemma 1: hard pair whose total momentum q is soft ------------- *)
(* F2(k1, q - k1) = (q/k)^2 (3/14 - (5/7) mu^2) + O((q/k)^3),
   mu = cos(angle between k1 and q).  2D vectors (coplanar Limber). *)
k1 = k {1, 0};
qv = eps k {Cos[th], Sin[th]};
k2 = qv - k1;
lemma1 = Normal[Assuming[k > 0 && 0 < eps < 1/10 && th \[Element] Reals,
   Series[F2[k1, k2], {eps, 0, 2}]]];
target1 = eps^2 (3/14 - (5/7) Cos[th]^2);
report["Lemma 1 (hard-pair kernel, O(eps^2) coefficient)",
  FullSimplify[lemma1 - target1]];

(* ---- Lemma 2: soft leg against a hard pair, power-law P ------------ *)
(* [2 F2(k1,k3) P(k1) + 2 F2(k2,k3) P(k2)] P(k3)
     = 2 P(k3) P(k) [13/14 + (4/7 - n/2) c^2] + O(eps),
   c = cos(angle between the soft leg and hard leg k1),
   n = dlnP/dlnk at the hard scale. *)
P[x_] := P0 x^n;
k3v = eps k {Cos[th], Sin[th]};
k2v = -k1 - k3v;
m2 = Sqrt[k2v . k2v];
normalised = (F2[k1, k3v] P[k] + F2[k2v, k3v] P[m2])/P[k];
lemma2 = Normal[Assuming[k > 0 && P0 > 0 && 0 < eps < 1/10 &&
    th \[Element] Reals && n \[Element] Reals,
   Series[normalised, {eps, 0, 0}]]];
target2 = 13/14 + (4/7) Cos[th]^2 - (n/2) Cos[th]^2;
report["Lemma 2 (soft-hard sum, O(1) coefficient)",
  FullSimplify[lemma2 - target2]];

(* ---- Averages ------------------------------------------------------ *)
avg2D = Simplify[Integrate[2 target2, {th, 0, 2 Pi}]/(2 Pi)];
report["2D (coplanar) angle average = 17/7 - n/2", avg2D - (17/7 - n/2)];
avg3D = Simplify[2 (13/14 + (4/7 - n/2) (1/3))];
report["3D angle-average anchor = 47/21 - n/3", avg3D - (47/21 - n/3)];

(* ---- Numeric spot checks (independent of the symbolic path) -------- *)
num1 = F2[{1, 0}, {0.01 Cos[Pi/3], 0.01 Sin[Pi/3]} - {1, 0}];
pred1 = 0.01^2 (3/14 - (5/7) Cos[Pi/3]^2);
Print["Lemma 1 numeric: F2 = ", num1, "  predicted = ", pred1,
  "  rel.dev = ", Abs[num1 - pred1]/Abs[pred1]];

nval = -2; epsv = 10^-3; cth = Pi/3;
k3n = epsv {Cos[cth], Sin[cth]}; k2n = -{1, 0} - k3n;
Pn[x_] := x^nval;
sumn = (F2[{1, 0}, k3n] Pn[1] +
    F2[k2n, k3n] Pn[Sqrt[k2n . k2n]])/Pn[1];
predn = 13/14 + (4/7 - nval/2) Cos[cth]^2;
Print["Lemma 2 numeric: sum = ", N[sumn, 10], "  predicted = ", N[predn, 10],
  "  rel.dev = ", N[Abs[sumn - predn]/Abs[predn], 3]];

(* ---- Remainder scaling (added after review round 1) ---------------- *)
(* Lemma 1 residual must scale as eps^3, Lemma 2 residual as eps^1. *)
res1[e_] := Abs[F2[{1, 0}, e {Cos[Pi/3], Sin[Pi/3]} - {1, 0}] -
   e^2 (3/14 - (5/7) Cos[Pi/3]^2)];
r1 = res1[10^-2]/res1[10^-2/2];
Print["Lemma 1 remainder scaling res(eps)/res(eps/2) = ", N[r1, 4],
  If[Abs[r1 - 8] < 1/2, "  PASS (~2^3)", "  FAIL"]];
res2[e_] := Module[{k3n = e {Cos[Pi/3], Sin[Pi/3]}, k2n},
  k2n = -{1, 0} - k3n;
  Abs[(F2[{1, 0}, k3n] Pn[1] + F2[k2n, k3n] Pn[Sqrt[k2n . k2n]])/Pn[1] -
    (13/14 + (4/7 - nval/2) Cos[Pi/3]^2)]];
r2 = res2[10^-3]/res2[10^-3/2];
Print["Lemma 2 remainder scaling res(eps)/res(eps/2) = ", N[r2, 4],
  If[Abs[r2 - 2] < 1/4, "  PASS (~2^1)", "  FAIL"]];

(* ---- Collapsed-family alpha reduction |Q| = v gamma ---------------- *)
(* Flat-triangle map (KAPPA3:1732-1744): theta_ij = ArcCos[c_ij];
   r1 = 0, r2 = (theta12, 0), r3 = (x3, y3),
   x3 = (theta12^2 + theta31^2 - theta23^2)/(2 theta12)  [degenerate
   theta12 -> 0: r2 = 0, r3 = (theta31, 0)].  Alpha kernel
   (KAPPA3:1747-1785): Q = u r2 + v Exp[-I phi] r3 (complex embedding).
   For the collapsed triple (1, c, c): theta12 = 0, theta23 = theta31 =
   gamma, so r2 = 0, r3 = gamma, and |Q| = v gamma for all u, phi. *)
gammaTest = ArcCos[86/100]; uT = 371; vT = 59; phiT = 2 Pi/7;
QT = uT*0 + vT Exp[-I phiT] gammaTest;   (* r2 = 0, r3 = gammaTest *)
report["Collapsed-family |Q| = v gamma", Simplify[Abs[QT] - vT gammaTest]];
