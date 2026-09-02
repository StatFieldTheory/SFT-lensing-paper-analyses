(* Exact toy fold, evaluated numerically at high precision.
   W(p)=Exp[-p^2/2], J(l)=int_l^inf W, so G(d)=int W(p)J(p+d)dp = Pi Erfc[d/2].
   The fold is  A(sigma) = int f(x) G(sigma x) dx  with f the density of the
   SHIFT that the regulator induces in the kink argument.                     *)
S[x_] := Piecewise[{{Exp[x]/2, x < 0}}, 1 - Exp[-x]/2];
G[d_] := Pi Erfc[d/2];
fMax[x_] := Exp[-Abs[x]] S[x];                  (* Max[U,V]  : both legs smeared *)
fSm[x_] := (1 + 2 Abs[x]) Exp[-2 Abs[x]]/2;     (* (U+V)/2   : SMOOTH control    *)
A3[s_] := NIntegrate[fMax[x] G[s x], {x, -Infinity, 0, Infinity},
                     WorkingPrecision -> 40];
A1[s_] := G[0]/2 + NIntegrate[Exp[-x]/2 G[s x], {x, 0, Infinity},
                              WorkingPrecision -> 40];
Asm[s_] := NIntegrate[fSm[x] G[s x], {x, -Infinity, 0, Infinity},
                      WorkingPrecision -> 40];
lead = N[Pi, 40];
Print["predicted linear coefficients  (deficit / (A(0) sigma)):"];
Print["   both legs smeared  E[Max[U,V]]/Sqrt[Pi] = 3/(4 Sqrt[Pi]) = ",
      N[3/(4 Sqrt[Pi]), 12]];
Print["   one  leg  smeared  E[Max[U,0]]/Sqrt[Pi] = 1/(2 Sqrt[Pi]) = ",
      N[1/(2 Sqrt[Pi]), 12]];
Print["   smooth control                                          = 0"];
Print[""];
Print[StringForm["`1`  `2`  `3`  `4`", PaddedForm["sigma", 10],
      PaddedForm["kink2/sig", 14], PaddedForm["kink1/sig", 14],
      PaddedForm["smooth/sig", 14]]];
Do[Print[StringForm["`1`  `2`  `3`  `4`",
   PaddedForm[N[sg, 8], {10, 7}],
   PaddedForm[N[(lead - A3[sg])/lead/sg, 10], {14, 10}],
   PaddedForm[N[(lead - A1[sg])/lead/sg, 10], {14, 10}],
   PaddedForm[N[(lead - Asm[sg])/lead/sg, 10], {14, 10}]]],
 {sg, {1/2, 1/4, 1/8, 1/16, 1/32, 1/64}}];
Print[""];
Print["fitted exponent d log(deficit) / d log(sigma), successive pairs:"];
Do[Module[{d1 = (lead - A3[sg])/lead, d2 = (lead - A3[sg/2])/lead,
          e1 = (lead - Asm[sg])/lead, e2 = (lead - Asm[sg/2])/lead},
   Print["   sigma ", N[sg, 6], " -> ", N[sg/2, 6],
         "   kink2 exponent = ", N[Log[d1/d2]/Log[2], 8],
         "   smooth exponent = ", N[Log[e1/e2]/Log[2], 8]]],
 {sg, {1/2, 1/4, 1/8, 1/16, 1/32}}];
