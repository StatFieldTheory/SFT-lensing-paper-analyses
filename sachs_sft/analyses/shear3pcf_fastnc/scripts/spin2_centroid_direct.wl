(* Direct centroid-projection natural-component construction.
   Confirm the per-vertex centroid reference angles alpha_j and the
   displacement vectors d2=X2-X1, d3=X3-X1 for fastnc's SAS isoceles family.
   We need the algebraic centroid angles to build the numeric reference. *)

(* SAS isoceles: equal sides t1=t2=g meet at the APEX (vertex 3) with opening
   angle phi.  Base vertices = 1 and 2, third side t3 = 2 g sin(phi/2).
   Place apex at origin, symmetric about +x axis.                          *)
g = 1;
X3 = {0, 0};                                  (* apex *)
X1 = g {Cos[phi/2],  Sin[phi/2]};             (* base vertex 1 *)
X2 = g {Cos[phi/2], -Sin[phi/2]};             (* base vertex 2 *)
cen = (X1 + X2 + X3)/3 // Simplify;

ang[v_] := ArcTan[v[[1]], v[[2]]];
a1 = Simplify[ang[cen - X1]];
a2 = Simplify[ang[cen - X2]];
a3 = Simplify[ang[cen - X3]];
Print["centroid = ", cen];
Print["alpha1 (X1->cen) = ", a1, "  = ", N[a1 /. phi -> Pi/3]];
Print["alpha2 (X2->cen) = ", a2, "  = ", N[a2 /. phi -> Pi/3]];
Print["alpha3 (X3->cen) = ", a3, "  = ", N[a3 /. phi -> Pi/3]];
Print["d2 = X2-X1 = ", Simplify[X2-X1]];
Print["d3 = X3-X1 = ", Simplify[X3-X1]];
(* side lengths sanity *)
Print["|X1-X3| = ", Simplify[Norm[X1-X3]], "  |X2-X3| = ", Simplify[Norm[X2-X3]],
      "  |X1-X2| = ", Simplify[Norm[X1-X2]], " (should be 2 sin(phi/2))"];
