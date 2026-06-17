(* ============================================================
   kappa3_prefactor_audit.wl
   ------------------------------------------------------------
   3-point analog of sigma2_prefactor_audit.wl.

   !!! SUPERSEDED CONCLUSION (2026-06-09, same day) !!!
   The "missing single dchi/dlambda=(1+z)^2" verdict BELOW IS WRONG.
   It assumed the equal-shell radial-decoherence delta is a CHI-density
   delta_D(chi-chi') while the fields are LAMBDA-densities -- a MEASURE MIX.
   The delta is a correlation OF the (1+z)^4-transfer lambda-density fields,
   so it inherits the LAMBDA measure: no separate Jacobian.
   Refuted numerically: born_limber 2-pt Cl with delta_measure="lambda" vs
   "chi" gives a ratio ~1 (NOT (1+z)^-2) -> the measures are EQUIVALENT, the
   delta-jacobian is internal binning bookkeeping fully compensated by the bin
   width. Also corroborated by the kappa3 lambda-native convention reproducing
   tree-level SPT <kappa^3> to ~30% (a (1+z)^6-mixed form fails by 4850x).
   VERDICT: NO missing Jacobian. The transfer operators (E^2/a^2)=(1+z)^4
   already impose the lambda-measure; the user's dchi/dlambda form and the
   code's lambda-native form are the SAME object in two measures. Operators
   and 1/k^2 cancellation (below) remain CORRECT. Kept for the record.
   ------------------------------------------------------------

   QUESTION (user, 2026-06-09): the equal-shell driving-field
   angular bispectrum B^{XYZ}_l(lambda) -- the "D4" vertex -- is
   claimed to be missing a single radial Jacobian dchi/dlambda=(1+z)^2,
   because Limber collapses the 3 line-of-sight chi-integrals to ONE,
   and converting that single surviving integral to the affine
   lambda-measure that the Sachs fold uses needs exactly one dchi/dlambda.

   This script does the FIRST-PRINCIPLES (1+z)-power + radial-measure
   accounting to confirm or refute that, anchored to the PyCCL-validated
   2-point (n=2) so the bookkeeping is not free-floating.

   Conventions (same as sigma2_prefactor_audit.wl):
     E = E_0/a, E_0=1, a=1/(1+z); dchi/dlambda = 1/a^2 = (1+z)^2.
     A(a) = -(3/2) Omega_m (H_0/c)^2 / a  = A_0 (1+z)   [Poisson amplitude]
     Phi(k) = A(a) delta(k)/k^2                          [Poisson]
     Phi_00 = (E^2/a^2) * [screen-Hessian] * Phi         [Sachs source]
              screen-Hessian (Ricci, spin-0)  = Lap_perp -> L^2/chi^2
              screen-Hessian (Weyl,  spin-2)  = eth^2     -> sqrt(L^2(L^2-2))/(2 chi^2)
   Run:
     ~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh \
        scripts/mathematica/kappa3_prefactor_audit.wl 120
   ============================================================ *)

$HistoryLength = 0;
Print["=== kappa3 (D4) prefactor + radial-measure audit ==="];
Print[""];

(* track (1+z) as the symbol opz := (1+z); chi, k, A0, Pd symbolic *)
opz = 1 + z;

(* ---- Phase 1: per-leg physical driving field Phi_00 in (1+z) powers ---- *)
Print["Phase 1. Per-leg physical Phi_00 (1+z)-power (sub-horizon/Limber):"];
Print["  (a) Sachs photon prefactor   E^2/a^2 = a^-4 = (1+z)^4"];
Print["  (b) screen Laplacian         Lap_perp -> L^2/chi^2 = k^2  (k=L/chi)"];
Print["  (c) Poisson                  Phi = A(a) delta/k^2,  A(a)=A0(1+z)"];
sachs = opz^4;          (* (a) *)
lap   = kk^2;           (* (b), with kk := k = L/chi *)
poisson = A0 opz / kk^2;(* (c) *)
Phi00leg = sachs * lap * poisson;            (* (E^2/a^2) k^2 (A0(1+z)/k^2) *)
Phi00leg = Simplify[Phi00leg];
Print["  Phi_00(leg) = (E^2/a^2)*Lap*Poisson*delta = ", Phi00leg, " * delta"];
Print["            => (1+z)^5 * A0 * delta   [k^2 from Laplacian cancels 1/k^2 from Poisson]"];
pLeg = Exponent[Phi00leg /. A0 -> 1, opz];   (* power of (1+z) per leg *)
Print["  per-leg (1+z) power pLeg = ", pLeg, "   (matches code: poisson_amp*(1+z)^4 = A0(1+z)*(1+z)^4)"];
Print[""];

(* ---- Phase 2: n-point equal-shell cumulant as a chi-density ---- *)
Print["Phase 2. n-point equal-shell cumulant density (physical, chi-density):"];
Print["  zeta_n^field(chi) = [prod_{i=1}^n Phi_00(leg)] * B_or_P / chi^{2(n-1)}"];
Print["    per-leg (1+z)^5  => total (1+z)^{5n}"];
Print["    angular Limber projection: 1/chi^2 per collapsed pair => 1/chi^{2(n-1)}"];
zetaFieldPow[n_] := 5 n;   (* (1+z) power of the pure field product *)
Print["    zeta_2 field (1+z) power = ", zetaFieldPow[2], "   (2-pt sigma^2 numerator)"];
Print["    zeta_3 field (1+z) power = ", zetaFieldPow[3], "   (3-pt vertex numerator)"];
Print[""];

(* ---- Phase 3: radial measure -- the crux ---- *)
Print["Phase 3. Radial measure: PHYSICAL (chi) vs CODE (lambda)."];
Print["  Physical observable (equal-shell):"];
Print["    O = INT prod_i dchi_i [efficiency] zeta_n(chi'_1..chi'_n)"];
Print["    equal-shell:  zeta_n ~ zeta_field(chi'_1) prod_{i=2}^n delta(chi'_1 - chi'_i)"];
Print["    [radial decoherence is in COMOVING distance chi -> deltas are in chi]"];
Print[""];
Print["  Convert every radial coordinate to affine lambda (what the Sachs fold uses):"];
Print["    INT dchi_i      = (1+z)^2 INT dlambda_i        [n integrals]"];
Print["    delta(chi-chi') = delta(lambda-lambda')/(1+z)^2 [n-1 deltas]"];
Print[""];
(* each collapsed (integral,delta) pair: (1+z)^2 * 1/(1+z)^2 = 1 *)
collapsedPairFactor = opz^2 * (1/opz^2);
collapsedPairFactor = Simplify[collapsedPairFactor];
Print["    each collapsed (integral x delta) pair -> ", collapsedPairFactor, "  (cancels)"];
Print["    the ONE surviving line-of-sight integral -> (1+z)^2  [single dchi/dlambda]"];
survivingJac = opz^2;
Print[""];
Print["  => PHYSICAL, written in lambda-measure:"];
Print["       O = INT dlambda_1 (1+z)^2 [efficiency] zeta_field(lambda_1)"];
Print["  => the CORRECT lambda-DENSITY vertex is"];
Print["       zeta_n^correct(lambda) = (1+z)^2 * zeta_field(lambda)"];
Print["  => n-INDEPENDENT single (1+z)^2, exactly the user's dchi/dlambda."];
Print[""];

(* ---- Phase 4: what the code's vertex carries ---- *)
Print["Phase 4. Code vertex (compute_kappa3_sigma3_high / equal_time_limber):"];
Print["    zeta_n^code(lambda) = zeta_field(lambda)   [radial_factor=1 for measure='lambda']"];
Print["    i.e. NO (1+z)^2 attached to the vertex."];
discrepancy = survivingJac;     (* correct/code = (1+z)^2 *)
Print["    correct/code = (1+z)^2  -> the standalone vertex UNDER-counts by (1+z)^2."];
Print["    THIS IS n-INDEPENDENT (same for the 2-pt sigma^2 and the 3-pt zeta)."];
Print[""];

(* ---- Phase 5: anchor to the PyCCL-validated 2-point ---- *)
Print["Phase 5. ANCHOR: the n=2 production 2-point is PyCCL-validated (0.2-0.4%)."];
Print["  Two ways the validated 2-pt can be right despite the vertex under-count:"];
Print["   (i)  the fold/response (corr_op C, R-propagator) SUPPLIES the single (1+z)^2,"];
Print["        OR"];
Print["   (ii) the 2-pt sigma^2 construction already bakes the (1+z)^2 into its"];
Print["        (1+z)-power (sigma2_prefactor_audit.wl: sigma^2_phys = (1+z)^8 =="];
Print["        (1+z)^{5*2} / (1+z)^2 -- the field numerator (1+z)^10 minus ONE (1+z)^2)."];
Print[""];
Print["  CROSS-CHECK against sigma2_prefactor_audit.wl:"];
twoPtFieldPow = zetaFieldPow[2];        (* 10 *)
sigma2Pow = twoPtFieldPow - 2;          (* 10 - 2 = 8 *)
Print["    2-pt field numerator (1+z) power = ", twoPtFieldPow];
Print["    minus single radial Jacobian (1+z)^2  => sigma^2 (1+z) power = ", sigma2Pow];
Print["    sigma2_prefactor_audit.wl independent result: sigma^2_phys ~ (1+z)^8  -> ",
      If[sigma2Pow === 8, "CONSISTENT (8==8): the single (1+z)^2 refund is REAL and present in the 2-pt", "MISMATCH"]];
Print[""];

(* ---- Phase 6: prediction for the 3-point ---- *)
Print["Phase 6. PREDICTION for the 3-point equal-shell cumulant (physical, lambda):"];
threePtFieldPow = zetaFieldPow[3];      (* 15 *)
zeta3Pow = threePtFieldPow - 2;         (* 15 - 2 = 13 *)
Print["    3-pt field numerator (1+z) power = ", threePtFieldPow];
Print["    minus single radial Jacobian (1+z)^2 => correct zeta_3 (1+z) power = ", zeta3Pow];
Print["    code zeta_3 (1+z) power (no refund)  = ", threePtFieldPow];
Print["    => the 3-pt vertex needs the SAME single (1+z)^2 as the 2-pt did."];
Print[""];

(* ---- Verdict ---- *)
Print["=================================================================="];
Print["VERDICT"];
Print["=================================================================="];
Print[""];
Print["1. OPERATORS O_X: CONFIRMED (driving_fields_harmonics.wl,"];
Print["   driving_fields_operator_form.wl, all tests PASS):"];
Print["     Phi_00 (Ricci)  = Lap_perp -> L^2/chi^2"];
Print["     Psi_0  (Weyl)   = eth^2    -> sqrt(L^2(L^2-2))/(2 chi^2)"];
Print["   The Laplacian L^2/chi^2 exactly cancels the Poisson 1/k^2 (k=L/chi),"];
Print["   so the user's prod[A/k^2] + O_X is the SAME object as the code's"];
Print["   ell-flat A(a) response -- the 1/k^2 is cancelled, not dropped."];
Print["   (spin-2: cancels only up to sqrt(L^2(L^2-2))/L^2 = 1+O(1/L^2);"];
Print["    that O(1/L^2) spin-2 residual is dropped by the ell-flat HIGH branch.)"];
Print[""];
Print["2. RADIAL MEASURE dchi/dlambda: the user's single (1+z)^2 is CORRECT and"];
Print["   n-INDEPENDENT. Limber collapses 3 chi-integrals to 1; converting that"];
Print["   single surviving line-of-sight integral to the affine lambda-measure"];
Print["   needs exactly ONE dchi/dlambda=(1+z)^2. The (n-1) equal-shell deltas"];
Print["   cancel (n-1) integral Jacobians; one integral survives -> one (1+z)^2."];
Print["   Cross-checked: 2-pt gives sigma^2 ~ (1+z)^8, matching the independent"];
Print["   sigma2_prefactor_audit.wl. The user's FORMULA is correct."];
Print[""];
Print["3. THE STANDALONE VERTEX (paper 'D4') as a pure lambda-density field"];
Print["   product UNDER-counts by this single (1+z)^2 -- so the paper's D4"];
Print["   expression SHOULD carry the dchi/dlambda (as the user wrote it)."];
Print[""];
Print["4. WHETHER THE PRODUCTION NUMBER IS WRONG depends on whether the FOLD"];
Print["   supplies the (1+z)^2: since n=2 is PyCCL-validated, the (1+z)^2 IS"];
Print["   supplied for the 2-pt (candidate: corr_op C build-time (1+z)^{+2})."];
Print["   ACTION: verify the corr_op C (1+z)^{+2} is applied ONCE in the FK"];
Print["   (C x K) topology -> if yes, production FK kk is correct and only the"];
Print["   PAPER formula needs the explicit dchi/dlambda; if no, the production"];
Print["   FK kk is under-counted by (1+z)^2 and must be fixed."];
Print[""];
Print["[done]"];
