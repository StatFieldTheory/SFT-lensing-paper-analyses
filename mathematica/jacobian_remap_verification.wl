(* ==================================================================
   jacobian_remap_verification.wl

   Numerical verification of the harmonic-space identities used in
   Appendix "First-Order Jacobian-Kernel Remapping of Angular
   Statistics" of the STF-lensing paper.

   The appendix claims:
   (i)   Spherical-harmonic orthogonality.
   (ii)  The "conjugated Gaunt integral"
            G[l,m,L,M,l',m'] ≡ ∫ Y*_{l,m} Y_{L,M} Y_{l',m'} dΩ
         equals (−1)^m √[(2l+1)(2L+1)(2l'+1)/(4π)]
                × 3j(l,L,l'; 0,0,0) × 3j(l,L,l'; −m,M,m').
   (iii) The gradient-gradient kernel
            F[l,m,L,M,l',m'] ≡ ∫ Y*_{l,m} (∇^A Y_{L,M})(∇_A Y_{l',m'}) dΩ
         equals ½[L(L+1) + l'(l'+1) − l(l+1)] × G[l,m,L,M,l',m'].

   Identity (iii) is the central claim of Appendix D: the transverse-
   Jacobian kernel (gradient-gradient coupling) reduces to the Gaunt
   kernel × a Casimir-difference factor. Combined with (ii), this
   gives the Wigner-3j-factored form of the CMB-lensing convolution
   used in Eq. (eq: appendix dCl ang).

   Reproduction:
     wolframscript -file scripts/jacobian_remap_verification.wl

   Run time: ~30 s on a laptop (numerical integration of a dozen
   spherical-harmonic integrals at low ℓ).
   ================================================================== *)

Print["==========================================================="];
Print[" STF-lensing: Jacobian-Kernel Remapping Verification"];
Print["==========================================================="];

ClearAll[Gaunt3j, GauntNum, GradGradNum, GradGradExpected, fmt];

(* ---------- Wigner-3j form of the conjugated Gaunt integral ---------- *)
(* ∫ Y*_{l,m} Y_{L,M} Y_{l',m'} dΩ, nonzero only when m = M + m'. *)
Gaunt3j[l_Integer, m_Integer, L_Integer, MM_Integer,
        lp_Integer, mp_Integer] :=
  (-1)^m * Sqrt[(2 l + 1) (2 L + 1) (2 lp + 1) / (4 Pi)] *
    ThreeJSymbol[{l, 0}, {L, 0}, {lp, 0}] *
    ThreeJSymbol[{l, -m}, {L, MM}, {lp, mp}];

(* ---------- Gaunt integral by direct numerical integration ---------- *)
GauntNum[l_Integer, m_Integer, L_Integer, MM_Integer,
         lp_Integer, mp_Integer] :=
  NIntegrate[
    Conjugate[SphericalHarmonicY[l, m, th, ph]] *
      SphericalHarmonicY[L, MM, th, ph] *
      SphericalHarmonicY[lp, mp, th, ph] * Sin[th],
    {th, 0, Pi}, {ph, 0, 2 Pi},
    WorkingPrecision -> 20, AccuracyGoal -> 10, PrecisionGoal -> 10];

(* ---------- Gradient-gradient integral by direct integration ---------- *)
(* On the unit 2-sphere with ds² = dθ² + sin²θ dφ², the dot product
     (∇^A f)(∇_A g) = ∂_θ f ∂_θ g + (1/sin²θ) ∂_φ f ∂_φ g
   and the volume element is sinθ dθ dφ. *)
GradGradNum[l_Integer, m_Integer, L_Integer, MM_Integer,
            lp_Integer, mp_Integer] :=
  NIntegrate[
    Conjugate[SphericalHarmonicY[l, m, th, ph]] *
      (D[SphericalHarmonicY[L, MM, th, ph], th] *
         D[SphericalHarmonicY[lp, mp, th, ph], th]
       + (1 / Sin[th]^2) *
         D[SphericalHarmonicY[L, MM, th, ph], ph] *
         D[SphericalHarmonicY[lp, mp, th, ph], ph]
      ) * Sin[th],
    {th, 0, Pi}, {ph, 0, 2 Pi},
    WorkingPrecision -> 20, AccuracyGoal -> 8, PrecisionGoal -> 8];

(* ---------- Gradient-gradient identity prediction ---------- *)
(* F = ½[L(L+1) + l'(l'+1) − l(l+1)] * G (conjugated Gaunt). *)
GradGradExpected[l_Integer, m_Integer, L_Integer, MM_Integer,
                 lp_Integer, mp_Integer] :=
  (1/2) * (L (L + 1) + lp (lp + 1) - l (l + 1)) *
    Gaunt3j[l, m, L, MM, lp, mp];

fmt[x_] := ToString@NumberForm[N[Chop[x, 10^-12], 6], {6, 4}];

(* -------------------------------------------------------------------
    TEST 1: Orthogonality of spherical harmonics
   ------------------------------------------------------------------- *)
Print["\n[Test 1] Spherical-harmonic orthogonality"];
Print["         ∫ Y*_{l,m} Y_{l',m'} dΩ = δ_{l,l'} δ_{m,m'}"];
Module[{cases, l, m, lp, mp, val, expected, err, pass, anyFail = False},
  cases = {
    {2, 1, 2, 1}, {2, 1, 3, 1}, {2, 1, 2, -1},
    {3, 0, 3, 0}, {4, 2, 4, 2}, {4, 2, 2, 2}};
  Do[
    {l, m, lp, mp} = c;
    val = NIntegrate[
      Conjugate[SphericalHarmonicY[l, m, th, ph]] *
        SphericalHarmonicY[lp, mp, th, ph] * Sin[th],
      {th, 0, Pi}, {ph, 0, 2 Pi}];
    expected = If[l == lp && m == mp, 1, 0];
    err = Abs[val - expected];
    pass = err < 10^-8;
    If[!pass, anyFail = True];
    Print["   (l,m,l',m') = ", c,
      "   num = ", fmt[val],
      "   expected = ", expected,
      "   |err| = ", fmt[err],
      If[pass, "   PASS", "   FAIL"]],
    {c, cases}];
  Print["   ", If[anyFail, "*** FAILURES in Test 1 ***", "All Test 1 checks passed."]]];

(* -------------------------------------------------------------------
    TEST 2: Gaunt integral vs Wigner 3j decomposition
   ------------------------------------------------------------------- *)
Print["\n[Test 2] Gaunt integral (conjugated): numerical vs Wigner 3j"];
Print["         ∫ Y*_{l,m} Y_{L,M} Y_{l',m'} dΩ = 3j-factored G"];
Module[{cases, l, m, L, MM, lp, mp, val, expected, err, pass, anyFail = False},
  (* All rows satisfy m = M + m' (else both sides are zero). *)
  cases = {
    {2, 1, 1, 0, 3, 1},
    {2, 0, 2, 0, 2, 0},
    {3, 2, 2, -1, 3, 3},
    {4, 0, 2, 0, 4, 0},
    {3, 1, 2, 0, 3, 1},
    {2, -1, 2, -1, 2, 0}};
  Do[
    {l, m, L, MM, lp, mp} = c;
    val = GauntNum[l, m, L, MM, lp, mp];
    expected = Gaunt3j[l, m, L, MM, lp, mp];
    err = Abs[val - expected];
    pass = err < 10^-6;
    If[!pass, anyFail = True];
    Print["   (l,m,L,M,l',m') = ", c,
      "   num = ", fmt[val],
      "   3j = ", fmt[expected],
      "   |err| = ", fmt[err],
      If[pass, "   PASS", "   FAIL"]],
    {c, cases}];
  Print["   ", If[anyFail, "*** FAILURES in Test 2 ***", "All Test 2 checks passed."]]];

(* -------------------------------------------------------------------
    TEST 3: Gradient-gradient kernel identity
       F = ½[L(L+1) + l'(l'+1) − l(l+1)] * G
   ------------------------------------------------------------------- *)
Print["\n[Test 3] Gradient-gradient kernel identity"];
Print["         F ≡ ∫ Y*_{l,m} (∇Y_{L,M}).(∇Y_{l',m'}) dΩ"];
Print["           = ½[L(L+1) + l'(l'+1) − l(l+1)] × G"];
Module[{cases, l, m, L, MM, lp, mp, val, expected, err, pass, anyFail = False},
  cases = {
    {2, 1, 1, 0, 3, 1},
    {2, 0, 2, 0, 2, 0},
    {3, 2, 2, -1, 3, 3},
    {4, 0, 2, 0, 4, 0},
    {3, 1, 2, 0, 3, 1},
    {3, 0, 1, 0, 2, 0},
    {2, -1, 2, -1, 2, 0}};
  Do[
    {l, m, L, MM, lp, mp} = c;
    val = GradGradNum[l, m, L, MM, lp, mp];
    expected = GradGradExpected[l, m, L, MM, lp, mp];
    err = Abs[val - expected];
    pass = err < 10^-4;
    If[!pass, anyFail = True];
    Print["   (l,m,L,M,l',m') = ", c,
      "   num = ", fmt[val],
      "   identity = ", fmt[expected],
      "   |err| = ", fmt[err],
      If[pass, "   PASS", "   FAIL"]],
    {c, cases}];
  Print["   ", If[anyFail, "*** FAILURES in Test 3 ***", "All Test 3 checks passed."]]];

(* -------------------------------------------------------------------
    TEST 4: Trivial edge cases (L = 0 ⇒ F = 0; l' = 0 ⇒ F = 0)
   ------------------------------------------------------------------- *)
Print["\n[Test 4] Trivial edge cases"];
Print["         L = 0  ⇒  ∇Y_{0,0} = 0  ⇒  F = 0 (and identity prefactor = l'(l'+1)−l(l+1))"];
Print["         l' = 0 ⇒  ∇Y_{0,0} = 0  ⇒  F = 0"];
Module[{pairs, l, m, val, pass, anyFail = False},
  (* L = 0 case: l = l' forces triangle inequality; else Gaunt = 0. *)
  pairs = {{2, 1}, {3, 0}, {4, 2}};
  Do[
    {l, m} = p;
    val = GradGradNum[l, m, 0, 0, l, m];  (* l' = l for non-zero Gaunt *)
    pass = Abs[val] < 10^-6;
    If[!pass, anyFail = True];
    Print["   F(l=", l, ",m=", m, ",L=0,M=0,l'=", l, ",m'=", m, ") = ",
      fmt[val], If[pass, "   PASS (≈0)", "   FAIL"]],
    {p, pairs}];
  Do[
    {l, m} = p;
    val = GradGradNum[l, m, l, m, 0, 0];  (* l'=0, m'=0; need M = m *)
    pass = Abs[val] < 10^-6;
    If[!pass, anyFail = True];
    Print["   F(l=", l, ",m=", m, ",L=", l, ",M=", m, ",l'=0,m'=0) = ",
      fmt[val], If[pass, "   PASS (≈0)", "   FAIL"]],
    {p, pairs}];
  Print["   ", If[anyFail, "*** FAILURES in Test 4 ***", "All Test 4 checks passed."]]];

Print["\n==========================================================="];
Print[" Verification complete."];
Print[" If all tests PASS, the identity"];
Print["   F[l,m,L,M,l',m'] = ½[L(L+1) + l'(l'+1) − l(l+1)] G[l,m,L,M,l',m']"];
Print[" used in Appendix D is numerically consistent to O(10^-4) or"];
Print[" better across the tested (l,m,L,M,l',m') tuples."];
Print["==========================================================="];
