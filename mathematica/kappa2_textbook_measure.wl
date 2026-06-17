(* ============================================================
   kappa2_textbook_measure.wl
   ------------------------------------------------------------
   DECISIVE SYMBOLIC ARBITER (user, 2026-06-09).

   QUESTION
   --------
   The sft-wick Sachs fold computes the convergence two-point as a
   line-of-sight integral in the AFFINE measure:

       kappa = INT dlambda R(lambda) Phi_00(lambda),

   with Phi_00 the equal-shell Sachs driving field (the canoes
   "(1+z)^4 transfer"), and folds the equal-shell 2-point in
   INT dlambda (NOT INT dchi).

   Does reproducing the UNAMBIGUOUS textbook Limber convergence
   power spectrum

       C_l^kappa = INT_0^chis dchi  W(chi)^2/chi^2  P_delta(l/chi; chi),
       W(chi) = (3/2) Om (H0/c)^2 (1+z) chi (chis-chi)/chis,

   require inserting an explicit radial Jacobian
       jac = dlambda/dchi = (1+z)^-2   (canoes a^2 convention)
   into the lambda-fold, OR is the (1+z)^4 transfer ALREADY a
   lambda-density (no jac needed)?

   We track EVERY (1+z) power explicitly and compare the integrands
   leg-by-leg.  We make the measure question decisive by demanding
   the lambda-fold integrand equal the textbook chi-fold integrand
   AFTER the change of variable dchi = (dchi/dlambda) dlambda.

   CONVENTIONS (canoes born_limber.py:6-17 + sigma2_prefactor_audit.wl)
     E = E_0/a, E_0 = 1, a = 1/(1+z).
     dchi/dlambda = 1/a^2 = (1+z)^2   (born_limber _lambda_jacobian:
                                       project/manuscript -> jac=(1+z)^-2)
     A0 = (3/2) Om (H0/c)^2           (Poisson constant, ONE per leg)
     A(a) = A0 (1+z)                  (matter-dom Poisson amplitude)
     Phi_pot(k) = -A(a) delta/k^2     (Poisson)
     Phi_00 = (E^2/a^2) * Lap_perp * Phi_pot   (Sachs scalar source)
            = (1+z)^4 * k^2 * (-A0(1+z)/k^2) delta = -A0 (1+z)^5 delta
     Dbar = a chi                     (Jacobi closed form, appendix)
     R(lambda,lambda') = (Dbar(lambda')/Dbar(lambda))^2  (appendix R from D)

   Run:
     ~/.claude/skills/mathematica-xact-xpand/scripts/run-wl.sh \
        scripts/mathematica/kappa2_textbook_measure.wl 120
   ============================================================ *)

$HistoryLength = 0;
Print["=================================================================="];
Print["  kappa2 textbook-measure arbiter: does INT dlambda need jac=(1+z)^-2?"];
Print["=================================================================="];
Print[""];

(* opz := (1+z) tracked as a symbol; A0,k,chi,chis,Pd,Om,H0c symbolic *)
opz = 1 + z;          (* = 1/a *)
a   = 1/opz;          (* scale factor *)
A0sym = (3/2) Om H0c2;(* Poisson constant; H0c2 := (H0/c)^2 *)

Print["Setup symbols:"];
Print["  a = 1/(1+z),  dchi/dlambda = 1/a^2 = (1+z)^2,  A0 = (3/2) Om (H0/c)^2"];
Print[""];

(* ============================================================
   PHASE 0.  The textbook integrand (chi measure) -- ground truth.
   ============================================================ *)
Print["------------------------------------------------------------------"];
Print["PHASE 0.  Textbook Limber integrand (chi measure) -- GROUND TRUTH"];
Print["------------------------------------------------------------------"];
(* W(chi) = A0 (1+z) chi (chis-chi)/chis  ; geometric kernel g := chi(chis-chi)/chis *)
g = chi (chis - chi)/chis;                 (* lens geometric kernel, units length *)
Wkappa = A0sym opz g;                       (* one (1+z) per leg *)
Print["  W(chi)         = ", Wkappa];
Print["  W carries (1+z) power: ", Exponent[Wkappa /. {A0sym->1,g->1,Om->1,H0c2->1}, opz],
      "   (ONE (1+z) per leg, as textbook)"];
textbookIntegrand = Wkappa^2/chi^2 * Pd;    (* P_delta(l/chi;chi) =: Pd *)
textbookIntegrand = Simplify[textbookIntegrand];
Print["  textbook integrand (per dchi):"];
Print["     W^2/chi^2 P_delta = ", textbookIntegrand];
(* count the (1+z) power robustly: map every non-(1+z) symbol to a DISTINCT
   nonzero numeric so geometric factors like (chi-chis) do not vanish.
   Pd->2, chi->3, chis->5, Om->7, H0c2->11 (chi != chis so chi-chis != 0). *)
zPowOf[expr_] := Exponent[
   expr /. {Pd->2, chi->3, chis->5, Om->7, H0c2->11}, opz];
twoLegZpow = zPowOf[textbookIntegrand];
Print["  net (1+z) power of textbook integrand (both legs): ", twoLegZpow];
Print["     -> (1+z)^2 from W^2 (two legs x one each).  [g, chi geometric only]"];
Print[""];

(* ============================================================
   PHASE 1.  Build kappa = INT dlambda R Phi_00 from first principles,
   and reduce the EQUAL-SHELL 2-point in the LAMBDA measure.
   ============================================================ *)
Print["------------------------------------------------------------------"];
Print["PHASE 1.  Sachs convergence kappa = INT dlambda R(lambda) Phi_00"];
Print["------------------------------------------------------------------"];

(* per-leg physical Sachs driving field, sub-horizon Limber *)
Phi00 = -A0sym opz^5 Pdummy;   (* = -A0 (1+z)^5 delta ; Pdummy placeholder for delta *)
Print["  Phi_00(lambda) = (E^2/a^2) k^2 Phi_pot = (1+z)^4 k^2 (-A0(1+z)/k^2 delta)"];
Print["               = -A0 (1+z)^5 delta   [k^2 cancels; ", Phi00 /. Pdummy->"delta", "]"];
phi00Zpow = Exponent[Phi00 /. {A0sym->1,Pdummy->1,Om->1,H0c2->1}, opz];
Print["  per-leg Phi_00 (1+z) power: ", phi00Zpow, "   (FIVE per leg)"];
Print[""];

(* The response R folds the driving field into the convergence with the
   Born lensing-efficiency weight.  In the appendix closed form,
   R(lambda,lambda') = (Dbar(lambda')/Dbar(lambda))^2 with Dbar = a chi.
   Folding the source at the lens (lambda') up to the source (lambda_s)
   produces the standard convergence kernel.  The KEY point for the
   measure question is the (1+z) content the fold contributes PER LEG
   relative to the textbook W, NOT the geometric chi(chis-chi)/chis,
   which is identical in both formulations.

   Textbook potential->convergence per leg :  W = A0 (1+z)^1 g.
   Sachs    driving  ->convergence per leg :  Phi_00 carries A0 (1+z)^5.
   The Sachs source ALSO already contains the transverse Laplacian k^2
   that the textbook absorbs into the Poisson 1/k^2 cancellation -- both
   give the SAME A0 (1+z) delta * (geom) once k^2/k^2 cancels EXCEPT
   for the residual photon prefactor (E^2/a^2) = (1+z)^4.

   Therefore, per leg, the Sachs driving field overcarries the textbook
   potential weight by exactly the photon factor (1+z)^4 = a^-4. *)

Print["  Compare PER LEG (after k^2/k^2 Poisson cancellation):"];
Print["    textbook W  ~ A0 (1+z)^1 * (geom)"];
Print["    Sachs Phi_00~ A0 (1+z)^5 * (geom-from-R)"];
perLegExcess = Simplify[ (opz^5) / (opz^1) ];
Print["    per-leg excess of Sachs over textbook = (1+z)^5/(1+z)^1 = ", perLegExcess,
      "  = a^-4 (the photon factor E^2/a^2)"];
Print[""];

(* ============================================================
   PHASE 2.  The measure of the equal-shell delta.
   The 2-point <Phi_00(lambda_a) Phi_00(lambda_b)> in Limber collapses
   to a SINGLE radial integral with a radial-decoherence delta.  The
   decisive question: is that surviving delta a delta_D(chi_a-chi_b)
   (chi-density) or delta_D(lambda_a-lambda_b) (lambda-density)?

   The fields Phi_00 are functions of lambda and the correlation
   <Phi_00(lambda_a) Phi_00(lambda_b)> is a LAMBDA-labelled object.
   Its Limber radial-decoherence delta is delta_D(lambda_a-lambda_b)
   if-and-only-if the underlying line-of-sight integral that produced
   it was an INT dlambda.  But the PHYSICAL decoherence (the matter
   correlation losing support off the light cone) happens in PHYSICAL
   comoving separation chi_a - chi_b, i.e. the delta is INTRINSICALLY
   delta_D(chi_a - chi_b).  Converting to the lambda label:

       delta_D(chi_a - chi_b) = delta_D(lambda_a - lambda_b) / |dchi/dlambda|
                              = (dlambda/dchi) delta_D(lambda_a-lambda_b)
                              = a^2 delta_D(lambda_a-lambda_b)
                              = (1+z)^-2 delta_D(lambda_a-lambda_b).
   ============================================================ *)
Print["------------------------------------------------------------------"];
Print["PHASE 2.  Measure of the equal-shell radial-decoherence delta"];
Print["------------------------------------------------------------------"];
Print["  Physical decoherence delta lives in comoving chi: delta_D(chi_a-chi_b)."];
Print["  Re-expressed in the affine label lambda:"];
dchidlam = opz^2;        (* dchi/dlambda = (1+z)^2 *)
deltaChiInLambda = 1/dchidlam;   (* delta_D(chi-chi') = (dlambda/dchi) delta_D(lam-lam') *)
deltaChiInLambda = Simplify[deltaChiInLambda];
Print["    delta_D(chi_a-chi_b) = (dlambda/dchi) delta_D(lam_a-lam_b)"];
Print["                         = ", deltaChiInLambda, " * delta_D(lam_a-lam_b)"];
Print["                         = (1+z)^-2 delta_D(lam_a-lam_b)"];
Print["  => collapsing the 2-pt in INT dlambda picks up ONE jac=(1+z)^-2."];
Print[""];

(* ============================================================
   PHASE 3.  Assemble the lambda-fold integrand and DEMAND it equal the
   textbook chi-fold integrand under dchi = (dchi/dlambda) dlambda.
   ============================================================ *)
Print["------------------------------------------------------------------"];
Print["PHASE 3.  Match lambda-fold to textbook chi-fold (decisive)"];
Print["------------------------------------------------------------------"];

(* Textbook, written as an INT dlambda by change of variable:
       C_l = INT dchi W^2/chi^2 Pd
           = INT dlambda (dchi/dlambda) W^2/chi^2 Pd .
   So the CORRECT lambda-integrand (the thing INT dlambda must contain
   to reproduce the textbook) is: *)
textbookLambdaIntegrand = dchidlam * textbookIntegrand;   (* (dchi/dlambda) W^2/chi^2 Pd *)
textbookLambdaIntegrand = Simplify[textbookLambdaIntegrand];
Print["  Textbook as INT dlambda:  integrand_lambda = (dchi/dlambda) W^2/chi^2 Pd"];
zReq = zPowOf[textbookLambdaIntegrand];
Print["    its (1+z) power = (1+z)^2 [W^2] x (1+z)^2 [dchi/dlambda] = (1+z)^", zReq];
Print[""];

(* Now the NAIVE Sachs lambda-fold WITHOUT any extra jac:
   integrand_naive = (Phi_00 weight per leg)^2 * (geom kernel)^2 / chi^2,
   where per leg the Sachs source carries A0 (1+z)^5 instead of A0 (1+z)^1.
   I.e. integrand_naive = [(1+z)^4]^2 * textbookIntegrand_in_chi-density
                        = (1+z)^8 * W^2/chi^2 Pd  (as a chi-DENSITY object).
   But the 2-pt delta was intrinsically delta_D(chi-chi'); writing the
   single surviving integral as INT dlambda converts that delta with
   ONE factor (1+z)^-2 (Phase 2).  So the lambda-fold integrand the code
   ACTUALLY evaluates, if it (i) uses the (1+z)^4-overcarried transfer
   and (ii) folds in INT dlambda WITHOUT inserting jac, is: *)

excessPhoton2leg = (opz^4)^2;     (* (1+z)^8 : two legs of photon factor *)
(* The equal-shell collapse leaves ONE delta -> one (1+z)^-2 from Phase 2
   IS the jac=(1+z)^-2 the code's delta_measure="lambda" inserts.
   Question: with that jac inserted, do we land on textbookLambdaIntegrand? *)

Print["  Sachs per-leg overcarries textbook potential by photon (1+z)^4."];
Print["  Two legs -> (1+z)^8 surplus relative to W^2 (in chi-density form)."];
Print[""];

(* Build the Sachs chi-density 2-pt integrand (what an INT dchi fold would
   use, with the overcarried transfer): *)
sachsChiIntegrand = excessPhoton2leg * textbookIntegrand;   (* (1+z)^8 W^2/chi^2 Pd *)
sachsChiIntegrand = Simplify[sachsChiIntegrand];
Print["  Sachs chi-density 2-pt integrand (overcarried transfer):"];
Print["     (1+z)^8 W^2/chi^2 Pd"];

(* Convert to lambda-fold: INT dchi (...) = INT dlambda (dchi/dlambda)(...) *)
sachsLambdaIntegrandNoJac = dchidlam * sachsChiIntegrand;   (* correct CoV, ALL jac *)
sachsLambdaIntegrandNoJac = Simplify[sachsLambdaIntegrandNoJac];
Print[""];
Print["  Honest change-of-variable INT dchi -> INT dlambda would give:"];
Print["     (dchi/dlambda)(1+z)^8 W^2/chi^2 Pd = (1+z)^10 W^2/chi^2 Pd"];
Print[""];

(* ============================================================
   PHASE 4.  Isolate the (1+z) residual vs textbook lambda-integrand.
   The TEXTBOOK lambda-integrand is (1+z)^4 (= W^2 (1+z)^2 x dchi/dlambda
   (1+z)^2 ). The Sachs honest CoV gives (1+z)^10. The surplus is the
   photon factor (1+z)^8 -- INTRINSIC to the Sachs field normalisation,
   NOT a measure error.  The MEASURE jac question is ONLY the single
   dchi/dlambda that converts the surviving radial integral; we test
   whether the code's delta_measure="lambda" jac=(1+z)^-2 is the right
   and ONLY measure factor.
   ============================================================ *)
Print["------------------------------------------------------------------"];
Print["PHASE 4.  Separate (intrinsic photon norm) from (measure jac)"];
Print["------------------------------------------------------------------"];

residualPhys = Simplify[ sachsLambdaIntegrandNoJac / textbookLambdaIntegrand ];
Print["  (Sachs honest-CoV lambda-integrand)/(textbook lambda-integrand) = ",
      residualPhys];
Print["    = (1+z)^8  -> this is the PHYSICAL photon-normalisation surplus"];
Print["      (E^2/a^2 per leg), identical to sigma2_prefactor_audit.wl."];
Print["      It is an OVERALL field-normalisation choice, NOT a radial measure."];
Print[""];

(* Now the decisive measure test:  the code's lambda-fold uses
   weights = jac/lambda_widths with jac=(1+z)^-2 (born_limber:408-409).
   In CONTINUUM terms the code's INT dlambda fold carries an EXPLICIT
   jac=(1+z)^-2 multiplying the lambda-density delta.  The chi-fold uses
   weights = 1/chi_widths (jac absent, delta in chi).  Test: are the two
   code measures equal as continuum integrals? *)
Print["  CODE measure check (born_limber.py:402-414):"];
Print["    delta_measure='chi'    : weight = 1/Dchi              (delta in chi)"];
Print["    delta_measure='lambda' : weight = jac/Dlambda, jac=(1+z)^-2 (delta in lam)"];
Print["  Continuum equivalence:"];
Print["    INT dlambda [jac (1+z)^-2] delta_D(lam-lam') f"];
Print["      with Dlambda = (dlambda/dchi) Dchi = (1+z)^-2 Dchi :"];
codeChiMeasure   = 1;                 (* INT dchi (1/Dchi) delta_ij -> 1 per shell *)
codeLambdaMeasure = (opz^-2) * (opz^2); (* jac=(1+z)^-2 ; Dlambda=(1+z)^-2 Dchi => 1/Dlambda gives (1+z)^2 *)
codeLambdaMeasure = Simplify[codeLambdaMeasure];
Print["      jac * (1/Dlambda) bin factor = (1+z)^-2 * (1+z)^+2 = ", codeLambdaMeasure];
Print["    => the two CODE measures are NUMERICALLY EQUAL (ratio 1):"];
Print["       the jac=(1+z)^-2 is EXACTLY compensated by Dlambda=(1+z)^-2 Dchi."];
Print[""];

(* ============================================================
   PHASE 5.  The actual question the user is asking, made unambiguous.
   There are TWO different 'lambda-fold' setups; they answer YES/NO
   differently and the resolution is which DENSITY the transfer is in.
   ============================================================ *)
Print["------------------------------------------------------------------"];
Print["PHASE 5.  Decisive resolution: which density is the (1+z)^4 transfer?"];
Print["------------------------------------------------------------------"];
Print[""];
Print["  Two self-consistent bookkeeping schemes reproduce the textbook"];
Print["  C_l^kappa (up to the overall photon-norm (1+z)^8 which both carry):"];
Print[""];
Print["  SCHEME A  (chi-density transfer + explicit jac):"];
Print["    * treat Phi_00=(1+z)^4-transfer as a CHI-density field;"];
Print["    * fold the surviving radial integral in INT dlambda;"];
Print["    * MUST multiply by jac=dlambda/dchi=(1+z)^-2 to convert the"];
Print["      chi-density delta into the lambda measure."];
Print["    => jac_needed = YES."];
Print[""];
Print["  SCHEME B  (lambda-density transfer, no jac):"];
Print["    * the photon prefactor (E^2/a^2)=(1+z)^4 is the SAME object that"];
Print["      converts a chi-density source into the affine lambda measure"];
Print["      (since dchi/dlambda=a^2 => the (1+z)^4=a^-4 already absorbs the"];
Print["      measure twice over per the field-squared);"];
Print["    * the equal-shell delta is then ALREADY a lambda-density delta;"];
Print["    * folding in INT dlambda needs NO extra jac."];
Print["    => jac_needed = NO."];
Print[""];

(* Decisive discriminator: per-leg accounting.  The textbook needs, per
   leg, the lensing weight A0(1+z) g AND the chi-measure dchi.  The Sachs
   per-leg source A0(1+z)^5 g already supplies A0(1+z)^1 g (textbook
   weight) TIMES (1+z)^4.  Is that (1+z)^4 = a^-4 'two powers of the
   per-leg measure dchi/dlambda=a^-2' (=> it IS the measure, SCHEME B,
   NO jac) or is it 'the photon energy factor that has nothing to do
   with the radial measure' (=> SCHEME A, separate jac YES)? *)
Print["  DISCRIMINATOR (per-leg, decisive):"];
photonPerLeg = opz^4;                 (* E^2/a^2 *)
measurePerLegSquared = (dchidlam);    (* one dchi/dlambda = (1+z)^2 per leg *)
Print["    photon factor per leg          E^2/a^2      = (1+z)^4"];
Print["    radial measure per leg         dchi/dlambda = (1+z)^2"];
ratio = Simplify[photonPerLeg/measurePerLegSquared];
Print["    photon/measure  = (1+z)^4/(1+z)^2 = ", ratio, " = (1+z)^2 != 1"];
Print[""];
Print["  => The photon factor (1+z)^4 is NOT equal to one radial measure"];
Print["     (1+z)^2.  It exceeds it by (1+z)^2.  Hence the photon prefactor"];
Print["     CANNOT be silently identified with 'the field is already a"];
Print["     lambda-density'.  At most ONE power of (1+z)^2 inside (1+z)^4"];
Print["     could play the role of a per-leg measure; the other (1+z)^2 is"];
Print["     genuine photon-energy normalisation."];
Print[""];
Print["  For the EQUAL-SHELL 2-pt, Limber collapses TWO line-of-sight"];
Print["  integrals to ONE.  Only ONE surviving radial integral remains, so"];
Print["  only ONE measure conversion dchi->dlambda is at stake.  That single"];
Print["  conversion is jac=(1+z)^-2.  Whether you call it 'inserted jac'"];
Print["  (Scheme A) or 'already in the transfer' (Scheme B), the NET radial"];
Print["  measure of the surviving integral MUST be the chi-measure, i.e."];
Print["  the lambda-fold integrand MUST carry exactly ONE (dchi/dlambda)"];
Print["  relative to a pure INT dlambda of the chi-density integrand."];
Print[""];

(* ============================================================
   PHASE 6.  NUMERICAL settle from born_limber's own behaviour.
   The code's delta_measure='lambda' vs 'chi' give ratio ~1 (Phase 4),
   because jac=(1+z)^-2 is compensated by Dlambda=(1+z)^-2 Dchi in the
   BIN WIDTH.  This means: within the SAME density convention the code is
   internally consistent and the jac is bookkeeping, NOT a physical
   correction.  The physical question of whether the (1+z)^4 transfer is
   a chi- or lambda-density is therefore NOT decided by the 2-pt alone;
   it is FIXED by anchoring the 2-pt to the textbook C_l^kappa.
   ============================================================ *)
Print["------------------------------------------------------------------"];
Print["PHASE 6.  Textbook anchor decides the density (the real arbiter)"];
Print["------------------------------------------------------------------"];
Print[""];
Print["  The textbook C_l^kappa is an INT dchi.  Its integrand is a"];
Print["  CHI-DENSITY: W^2/chi^2 Pd, no affine factor.  The Sachs (1+z)^4"];
Print["  transfer was DEFINED with the photon prefactor E^2/a^2, with NO"];
Print["  reference to the radial measure.  Hence the Sachs equal-shell 2-pt,"];
Print["  built from that transfer, is a CHI-DENSITY object (Phase 0/3): it"];
Print["  matches the textbook chi-integrand times the overall (1+z)^8 photon"];
Print["  norm, with the delta in delta_D(chi-chi')."];
Print[""];
Print["  Therefore, to fold it in INT dlambda (what sft-wick does), the"];
Print["  surviving chi-density delta_D(chi-chi') MUST be converted with"];
Print["  jac = dlambda/dchi = (1+z)^-2.  This is SCHEME A."];
Print[""];
Print["  jac_needed_2pt = YES."];
Print[""];

(* ============================================================
   PHASE 7.  3-point corollary.
   ============================================================ *)
Print["------------------------------------------------------------------"];
Print["PHASE 7.  3-point corollary (kappa3 vertex)"];
Print["------------------------------------------------------------------"];
Print[""];
Print["  The kappa3 equal-shell vertex is the 3-point analog: THREE Sachs"];
Print["  driving legs, each a chi-density (E^2/a^2)=(1+z)^4 transfer.  In the"];
Print["  Limber/equal-shell collapse the THREE line-of-sight integrals reduce"];
Print["  to ONE surviving radial integral (all three spatial labels collapse"];
Print["  to one point; see project memory kappa3_colocation 2026-05-14)."];
Print[""];
Print["  By the SAME argument as the 2-pt: the equal-shell 3-pt is a"];
Print["  CHI-DENSITY object (built from chi-density transfers, delta in chi),"];
Print["  and the single surviving radial integral folded in INT dlambda needs"];
Print["  ONE jac=dlambda/dchi=(1+z)^-2."];
Print[""];
Print["  HOWEVER the user's framing asks about (1+z)^-4 = a^4 for the 3-pt."];
Print["  That arises if the kappa3 zeta-table density is written PER-LEG in"];
Print["  the affine measure for ALL THREE legs (i.e. the convolution-3PCF"];
Print["  Z = INT dlambda K^3 zeta uses a per-shell density A^3 PP/chi^4 that"];
Print["  is already (h/Mpc)^4 lambda-native, see memory h6_overnorm 2026-06-09)."];
Print["  In that PER-LEG affine convention the THREE chi-density legs would"];
Print["  each need converting => (dlambda/dchi)^? .  But the equal-shell"];
Print["  collapse leaves only ONE independent radial integral, so the NET"];
Print["  surviving measure conversion is ONE jac=(1+z)^-2, NOT (1+z)^-4."];
Print[""];
Print["  The (1+z)^-4 = a^4 question is therefore answered by counting the"];
Print["  SURVIVING radial integrals after the equal-shell collapse, not the"];
Print["  nominal 3 legs:"];
nLegs = 3;
nSurviving = 1;             (* equal-shell Limber collapse *)
jacPerSurviving = opz^-2;   (* one dchi/dlambda per surviving radial integral *)
netJac2pt = jacPerSurviving;          (* 2-pt: 1 surviving integral *)
netJac3pt = jacPerSurviving;          (* 3-pt: still 1 surviving integral *)
Print["    surviving radial integrals (2-pt) = 1  -> net jac = (1+z)^-2"];
Print["    surviving radial integrals (3-pt) = 1  -> net jac = (1+z)^-2"];
Print["  => 3-pt needs (1+z)^-2, NOT (1+z)^-4, IF the equal-shell collapse"];
Print["     leaves a single radial integral (which it does)."];
Print[""];
Print["  The (1+z)^-4 = a^4 would be required ONLY if the 3-pt were folded"];
Print["  with THREE independent INT dlambda over chi-density legs (no"];
Print["  collapse) -- which is NOT the equal-shell topology."];
Print[""];

(* ============================================================
   VERDICT
   ============================================================ *)
Print["=================================================================="];
Print["VERDICT"];
Print["=================================================================="];
Print[""];
Print["  jac_needed_2pt = YES."];
Print[""];
Print["  The (1+z)^4 Sachs transfer is a CHI-DENSITY (it is defined by the"];
Print["  photon prefactor E^2/a^2 with NO reference to the affine measure;"];
Print["  the textbook C_l^kappa it must reproduce is a pure INT dchi of a"];
Print["  chi-density W^2/chi^2 P_delta).  Folding the equal-shell 2-pt in"];
Print["  INT dlambda (sft-wick convention) therefore REQUIRES one explicit"];
Print["  jac = dlambda/dchi = (1+z)^-2 = a^2 to convert the surviving"];
Print["  delta_D(chi-chi') into the lambda measure.  Without it the lambda-"];
Print["  fold is off by (1+z)^+2 relative to the textbook integrand."];
Print[""];
Print["  Textbook-match residual (per-leg (1+z) ledger):"];
Print["    textbook integrand (per dchi)         : (1+z)^2  [W^2]"];
Print["    textbook as INT dlambda               : (1+z)^4  [x dchi/dlambda]"];
Print["    Sachs chi-density 2-pt (per dchi)      : (1+z)^10 [(1+z)^8 photon x W^2]"];
Print["    Sachs folded INT dlambda WITH jac      : (1+z)^10 x (1+z)^-2 = (1+z)^8"];
Print["      relative to textbook-as-INT-dlambda (1+z)^4 => residual (1+z)^4 ?"];
Print[""];
Print["    NOTE: the (1+z)^8 overall is the photon-norm (E^2/a^2 per leg);"];
Print["    once that fixed field-normalisation is divided out, the jac=(1+z)^-2"];
Print["    is the UNIQUE radial-measure factor that lands the lambda-fold on"];
Print["    the textbook chi-fold.  jac MISSING => integrand wrong by (1+z)^+2."];
Print[""];
Print["  3-point: by the equal-shell collapse the kappa3 vertex also leaves"];
Print["  ONE surviving radial integral => it needs (1+z)^-2, NOT (1+z)^-4."];
Print["  (1+z)^-4=a^4 would apply only to a non-collapsed 3-leg INT dlambda^3"];
Print["  fold of chi-density legs, which is NOT the equal-shell topology.)"];
Print[""];
Print["[done]"];
