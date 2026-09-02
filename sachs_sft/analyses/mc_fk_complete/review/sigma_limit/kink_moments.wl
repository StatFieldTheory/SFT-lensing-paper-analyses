(* ==========================================================================
   Why the FK estimator's regulator error is O(sigma_lambda), not O(sigma^2).

   After the exact Wick identity (verified to 1e-16 in identity_check.py) the
   estimator's expectation IS the FK diagram of the colored driving field:

     <k2_A k1_B> = sum_j dlam W_j F_Abc sum_{m,n<=j-1} U[j,m] U[j,n]
                   sum_p dlam W_p <f_b[m] f_c[n] f_B[p]>_c .

   As a function of the three drive positions the kernel is

     T(m,n,p) = W(p) D(m)^2 D(n)^2 J(Max[m,n]),   J(l) = int_l^lf W/D^4,
     J'(l) = -W(l)/D(l)^4,   H = D^4 J,

   the Max coming from the PRODUCT of the two retarded response legs
   theta(m<j) theta(n<j).  The white-noise limit collapses the induced cumulant
   to zeta(l) delta delta and gives  sum dlam W H F:zeta -- the local FK
   integral.  At finite sigma each delta is the AR(1) kernel k(u)=Exp[-|u|]/2 of
   width sigma, so the collapse point is displaced by sigma * (a shift variable),
   and the leading error is sigma * E[shift] * J'.

     Q on the far observable leg  -> BOTH F-legs smeared -> shift = Max[U,V]
     Q on an F-vertex input leg   -> ONE  F-leg  smeared -> shift = Max[U,0]

   A SMOOTH kernel argument would give E[shift]=0 by the two-sided exponential's
   parity, hence O(sigma^2).  The Max is a kink, and its mean is not zero.
   ========================================================================== *)

k[u_] := Exp[-Abs[u]]/2;

Print["=== 1. AR(1) smearing kernel: normalisation and moments ==="];
Print["  int k     = ", Integrate[k[u], {u, -Infinity, Infinity}]];
Print["  int u k   = ", Integrate[u k[u], {u, -Infinity, Infinity}],
      "   <- ZERO: any SMOOTH direction costs only O(sigma^2)"];
Print["  int u^2 k = ", Integrate[u^2 k[u], {u, -Infinity, Infinity}]];

Print["=== 2. the two kink moments (the coefficients of the O(sigma) term) ==="];
muMax = Integrate[k[u] k[v] Max[u, v], {u, -Infinity, Infinity},
                  {v, -Infinity, Infinity}];
muHalf = Integrate[k[u] Max[u, 0], {u, -Infinity, Infinity}];
Print["  E[Max[U,V]] = ", muMax, "   -> T1 deficit = (3/4) sigma * sum dlam W^2 F:X3"];
Print["  E[Max[U,0]] = ", muHalf, "   -> T2 deficit = (1/2) sigma * sum dlam W^2 F:(X1+X2)"];
Print["  E[Max[U,V]^2] = ",
      Integrate[k[u] k[v] Max[u, v]^2, {u, -Infinity, Infinity},
                {v, -Infinity, Infinity}],
      " , E[Max[U,0]^2] = ", Integrate[k[u] Max[u, 0]^2, {u, -Infinity, Infinity}],
      "   <- finite: the next term is sigma^2, with NO sigma^2 Log sigma"];

Print["=== 3. the same constants in smoothed-step language ==="];
S[x_] := Piecewise[{{Exp[x]/2, x < 0}}, 1 - Exp[-x]/2];   (* int_-inf^x k *)
i1 = Integrate[S[x] - UnitStep[x], {x, -Infinity, Infinity}];
i2 = Integrate[S[x]^2 - UnitStep[x], {x, -Infinity, Infinity}];
i3 = Integrate[UnitStep[x] (S[x] - 1), {x, -Infinity, Infinity}];
Print["  int (S - theta)       = ", i1,
      "   <- ONE smoothed causal step alone is free at O(sigma)"];
Print["  int (S^2 - theta)     = ", i2,
      "   <- theta is idempotent, S is NOT: this is the whole effect"];
Print["  int theta (S - theta) = ", i3];
Print["  -int(S^2-theta) == E[Max[U,V]] : ", Simplify[-i2 - muMax] === 0,
      " ;  -int theta(S-1) == E[Max[U,0]] : ", Simplify[-i3 - muHalf] === 0];

Print["=== 4. universality: O(sigma) cannot be tuned away by a nicer kernel ==="];
Print["  For ANY symmetric regulator kernel with iid draws U,V:"];
Print["    E[Max[U,V]] = E[|U-V|]/2 > 0  unless U = V almost surely, and"];
Print["    E[Max[U,0]] = E[U+]        > 0  unless U = 0 almost surely."];
Print["  Both are strictly positive for every non-degenerate kernel, so the"];
Print["  leading regulator error is O(sigma) with a DEFICIT sign for every"];
Print["  choice of colored noise; only the two constants change."];
Print["  identity Max[u,v] == (u+v+|u-v|)/2 : ",
      Simplify[Max[u, v] - (u + v + Abs[u - v])/2,
               Assumptions -> {u \[Element] Reals, v \[Element] Reals}] === 0];
