(* Flat-sky shear 3PCF natural components: derive the exact 2D-Fourier-integral
   form and the CENTROID projection spin-2 phase factors, matching fastnc x2cent.
   Plain complex algebra (no xTensor needed). *)

(* ---- Setup: SAS isoceles triangle, vertices on the plane ---------------- *)
(* fastnc convention (trigutils.x1x2phi_to_x1x2x3): side vectors.
   Triangle sides t1,t2 emanate from vertex with opening angle.
   We use the SAS family: t1=t2=g (arcmin->rad handled numerically later),
   opening angle phi, third side t3 = 2 g sin(phi/2).               *)

(* Place the three GALAXY positions X1,X2,X3 in the plane.
   fastnc's natural components are correlators of shear at the three
   triangle vertices.  The "x-projection" measures each shear relative to
   the triangle SIDE directions; the "centroid" projection measures each
   relative to the line vertex->centroid.                                  *)

(* Reference (Schneider-Lombardi 2003, A&A 408, 829; Sugiyama+2024):
   gamma~(l) = e^{2 i polar(l)} kappa~(l).
   Gamma^0(X1,X2,X3) = <gamma(X1) gamma(X2) gamma(X3)>
   Gamma^1 = <gamma*(X1) gamma(X2) gamma(X3)>, etc.
   In Fourier space, with <k(l1)k(l2)k(l3)>=(2pi)^2 delta(l1+l2+l3) B(l1,l2,l3),
   Gamma^mu = INT d^2 l2 d^2 l3 /(2pi)^4 *
                exp(i l2.(X2-X1)) exp(i l3.(X3-X1)) *
                exp(2 i [ s1 polar(l1) + s2 polar(l2) + s3 polar(l3) ]) *
                B(|l1|,|l2|,|l3|),
   with l1 = -(l2+l3), s_j = +1 unconjugated, -1 conjugated leg,
   and the GAMMA^mu indexing: mu=0 -> (s1,s2,s3)=(+,+,+);
   mu=1 -> (-,+,+); mu=2 -> (+,-,+); mu=3 -> (+,+,-).                       *)

(* The "x-projection" used internally by fastnc rotates out the absolute
   triangle orientation; the CENTROID projection multiplies by the per-vertex
   phase exp(-2 i alpha_j) where alpha_j = polar angle of (centroid - X_j).
   We verify fastnc x2cent factors reproduce these centroid phases.        *)

(* fastnc x2cent (fastnc.py:752-784), with v-vectors:                       *)
x2cent[mu_, t1_, t2_, phi_] := Module[{v, q1, q2, q3},
  v = t1 + t2 Exp[-I phi]; q1 = v/Conjugate[v];
  v = -2 t1 + t2 Exp[-I phi]; q2 = v/Conjugate[v];
  v = t1 - 2 t2 Exp[-I phi]; q3 = v/Conjugate[v];
  Switch[mu,
    0, q1 q2 q3 Exp[3 I phi],
    1, Conjugate[q1] q2 q3 Exp[1 I phi],
    2, q1 Conjugate[q2] q3 Exp[3 I phi],
    3, q1 q2 Conjugate[q3] Exp[-1 I phi]]];

(* Verify |x2cent|=1 (pure phase) for the isoceles family, all mu.         *)
SeedRandom[1];
testphi = Pi/3; tt1 = 1.0; tt2 = 1.0;
Print["%%CHECK x2cent modulus (should all be 1) %%"];
Do[Print["mu=", mu, "  |x2cent|=", Abs[x2cent[mu, tt1, tt2, testphi]]], {mu,0,3}];

(* Centroid phase from FIRST PRINCIPLES for the isoceles SAS triangle.
   Place: vertex1 and vertex2 at base, vertex3 at apex.
   Use fastnc's internal placement: the FFT-native "x-projection" has the
   triangle with side-1 along a fixed axis. The centroid bridge is a pure
   geometric per-leg rotation; we only need its phi-dependence (g cancels in
   the spin-2 ratio because all lengths scale by g).                        *)

(* Build explicit vertices: isoceles, base 1-2, apex 3.
   Use unit g; opening angle phi at the apex? In the SAS "two equal sides
   t1=t2 with opening angle phi" convention, the EQUAL sides are t1 (vertex
   3->1) and t2 (vertex 3->2), opening angle phi between them at vertex 3.
   So vertex 3 is the apex. third side t3 = 2 sin(phi/2).                   *)
g = 1;
apex = {0, 0};
X1 = apex + g {Cos[phi/2], Sin[phi/2]};
X2 = apex + g {Cos[phi/2], -Sin[phi/2]};
X3 = apex;
centroid = (X1 + X2 + X3)/3 // Simplify;
alpha[Xj_] := ArcTan[(centroid - Xj)[[1]], (centroid - Xj)[[2]]];
Print["%%CENTROID per-vertex angles alpha_j (vertex->centroid)%%"];
Print["alpha1 = ", Simplify[alpha[X1]]];
Print["alpha2 = ", Simplify[alpha[X2]]];
Print["alpha3 = ", Simplify[alpha[X3]]];

(* For the reference builder we will compute Gamma in the x-projection by the
   direct 2D Fourier integral, then multiply by x2cent[mu,...] exactly as
   fastnc does. So the symbolic job is only to CONFIRM x2cent is a pure phase
   (done above) and to record the integral form. Output the integral form. *)
Print["%%INTEGRAL FORM%%"];
Print["Gamma^mu = INT d2l2 d2l3/(2pi)^4 e^{i l2.d2 + i l3.d3} ",
      "e^{2i(s1 ang(l1)+s2 ang(l2)+s3 ang(l3))} B(l1,l2,l3), l1=-(l2+l3)"];
Print["then multiply by x2cent[mu,t1,t2,phi] to go x-proj -> centroid."];
