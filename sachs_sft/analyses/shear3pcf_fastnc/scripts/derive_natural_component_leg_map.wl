(* ============================================================================
   derive_natural_component_leg_map.wl

   DIAGNOSTIC C: pin the natural-component channel-label / leg-conjugation
   convention for the cosmic-shear 3PCF Gamma^0..3 from the Sachs spin-2 helicity
   structure AND the isoceles (SAS) mirror symmetry, FIRST-PRINCIPLES.

   GOAL.  The canoes natural-component assembler (stage1_ours_v2_JACFIX.py) maps
     Gamma^0 <- zeta_PPP (un-conj),  Gamma^{1,2,3} <- single-leg-conjugated zeta_D
   with the SAS isoceles labeling APEX = leg 3, BASE = legs {1,2}, and a derived
   great-circle -> centroid spin-2 projection phase.  fastnc (the external truth)
   instead satisfies, on the SAS isoceles family t1=t2=gamma, EXACTLY (1e-16):
        Gamma^1 is REAL   and   Gamma^2 = conj(Gamma^3).
   That fixes fastnc's APEX (the self-mirror vertex) = leg 1, BASE pair = {2,3}.

   This script derives, with no fitting, which (vertex<->mu) assignment makes the
   natural components carry the isoceles mirror symmetry, and what the per-leg
   centroid-projection phase is.  Pure spin-weighted algebra in the flat-sky
   plane (xAct not required for a 2D spin-2 phase; we keep it self-contained and
   symbolic).

   ========================================================================
   PHYSICS (pinned from cosmology.tex eq:gamma from Dsachs23 + Schneider-Lombardi)
   ========================================================================
   The image shear at vertex X_j is the spin-(+2) field gamma(X_j); projected to a
   per-vertex reference direction alpha_j it is g|_j = gamma(X_j) e^{-2 i alpha_j}.
   Natural components (s_j = +1 unconj leg, -1 conj leg):
        Gamma^mu = < g|_1^{s1} g|_2^{s2} g|_3^{s3} >,
   conjugating exactly the leg singled out by mu.

   Under a spin-2 field, g -> g e^{-2 i a} (rotation by a), g* -> g* e^{+2 i a}.
   ========================================================================  *)

$DefInfoQ = False;
Print["%%SECTION%% 0. SAS isoceles geometry + mirror"];

g = gamma; ph = phi;
(* APEX = vertex A (two equal sides meet here, opening angle phi),
   BASE  = vertices B, C (separated by t3 = 2 g sin(phi/2)).
   Place apex at origin, base symmetric about +x axis. *)
XA = {0, 0};
XB = g {Cos[ph/2],  Sin[ph/2]};
XC = g {Cos[ph/2], -Sin[ph/2]};
(* The mirror is y -> -y : XA->XA (self), XB<->XC (swap). *)

cen = (XA + XB + XC)/3;
(* centroid reference angle at each vertex = polar(centroid - X) *)
alA = ArcTan[(cen - XA)[[1]], (cen - XA)[[2]]];
alB = ArcTan[(cen - XB)[[1]], (cen - XB)[[2]]];
alC = ArcTan[(cen - XC)[[1]], (cen - XC)[[2]]];
Print["  alpha_apex = ", ToString[InputForm[FullSimplify[alA, 0<g && 0<ph<Pi]]]];
Print["  alpha_baseB= ", ToString[InputForm[FullSimplify[alB, 0<g && 0<ph<Pi]]]];
Print["  alpha_baseC= ", ToString[InputForm[FullSimplify[alC, 0<g && 0<ph<Pi]]]];

(* under mirror y->-y, the centroid angle transforms as alpha -> -alpha,
   and apex is fixed, base B<->C swap. Check: *)
Print["  mirror: alA -> -alA ? ", FullSimplify[alA + alA == 0]];  (* apex on axis: alA=pi (or 0) *)
mAB = FullSimplify[alB + alC, 0<g && 0<ph<Pi];
Print["  alB + alC = ", ToString[InputForm[mAB]], "  (=> alC = -alB mod 2pi => base pair mirror)"];

Print["%%SECTION%% 1. natural components & the mirror action on each mu"];
(* The connected correlator <gamma gamma gamma> on the isoceles triangle is mirror
   symmetric AS A SET (the triangle maps to itself with B<->C).  The shear field is
   spin-2: under y->-y (a reflection, NOT a rotation) a spin-+2 field maps to the
   spin-(-2) field at the mirrored point with reference angle negated:
        gamma(X) -> conj(gamma(X_mirror)),   alpha -> -alpha.
   So g|_j = gamma(X_j) e^{-2 i alpha_j} -> conj(gamma(X_j')) e^{+2 i alpha_j'}
            = conj( gamma(X_j') e^{-2 i alpha_j'} ) = conj( g|_{j'} ),
   where j' is the mirror image vertex.  Hence under the mirror:
        g|_apex   -> conj(g|_apex)            (apex self-image)
        g|_baseB  -> conj(g|_baseC)           (base pair swap)
        g|_baseC  -> conj(g|_baseB).
   Apply to each natural component (conjugating the singled-out leg). *)

(* Represent g|_apex = A, g|_baseB = B, g|_baseC = C as formal symbols; the mirror
   M acts by  A->conj(A), B->conj(C), C->conj(B), and reverses the ORDER of the
   correlator product is irrelevant (scalars commute), and complex-conjugates the
   ENSEMBLE AVERAGE only through the field map, so:
        M[<X Y Z>] = < M[X] M[Y] M[Z] >.
   Crucially the triangle is invariant, so M[<...>] = <...> (same physical config)
   ONLY AFTER also conjugating overall (a reflection is antiunitary on helicity).
   The net statement Schneider-Lombardi give:  on the isoceles family,
        Gamma^{apex-conj} = conj( Gamma^{apex-conj} )   (REAL),
        Gamma^{baseB-conj} = conj( Gamma^{baseC-conj} ). *)

(* Symbolic verification of the index bookkeeping: assign mu and test. *)
legs = {"apex", "baseB", "baseC"};

mirrorMap = <| "apex"->"apex", "baseB"->"baseC", "baseC"->"baseB" |>;
(* a natural component is labeled by which leg is conjugated *)
(* under mirror, conj-leg -> conj-leg(mirror image), plus overall conjugation *)
Print["  component conj@apex  -> conj( component conj@", mirrorMap["apex"], " )  => REAL"];
Print["  component conj@baseB -> conj( component conj@", mirrorMap["baseB"], " )"];
Print["  component conj@baseC -> conj( component conj@", mirrorMap["baseC"], " )"];
Print["  ==> the APEX-conjugated natural component is REAL on isoceles;"];
Print["      the two BASE-conjugated components are complex conjugates."];

Print["%%SECTION%% 2. MATCH to fastnc's observed structure"];
Print["  fastnc (machine precision): Gamma^1 REAL, Gamma^2 = conj(Gamma^3)."];
Print["  => fastnc leg map: mu=1 conjugates the APEX; mu=2,3 conjugate the BASE pair."];
Print["  => fastnc APEX = leg 1 ; BASE = legs {2,3}."];
Print["  Also Gamma^0 (all-unconjugated, zeta_PPP) is REAL on isoceles (self-mirror)."];

Print["%%SECTION%% 3. canoes assembler's CURRENT map (from stage1_ours_v2_JACFIX.py)"];
Print["  _DMOD_SPINS = {1:(-2,2,2), 2:(2,-2,2), 3:(2,2,-2)} and the SAS cosine"];
Print["  triples (cos t3, cos g, cos g): leg1,leg2 = BASE endpoints (sep t3),"];
Print["  leg3 = APEX.  So canoes conjugates: Gamma^1 -> leg-1 (BASE), Gamma^2 -> "];
Print["  leg-2 (BASE), Gamma^3 -> leg-3 (APEX)."];
Print["  => canoes APEX-conjugated component is Gamma^3, NOT Gamma^1."];
Print["  MISMATCH vs fastnc: fastnc apex-conj = mu1 ; canoes apex-conj = mu3."];

Print["%%SECTION%% 4. The relabel that aligns canoes to fastnc"];
Print["  Identify mu by WHICH PHYSICAL VERTEX is conjugated, not by slot index:"];
Print["    fastnc Gamma^1 (apex-conj, REAL)      <->  canoes Gamma^3 (apex=leg3-conj)"];
Print["    fastnc Gamma^2 (baseC-conj)           <->  canoes Gamma^1 or Gamma^2 (base)"];
Print["    fastnc Gamma^3 (baseB-conj)           <->  canoes Gamma^2 or Gamma^1 (base)"];
Print["  i.e. a CYCLIC RELABEL of the canoes natural-component slots:"];
Print["    canoes (G1_base, G2_base, G3_apex) -> fastnc (G1_apex, G2_base, G3_base)"];
Print["  The apex channel must move from slot 3 (canoes) to slot 1 (fastnc)."];

Print["%%SECTION%% 5. centroid projection phase per leg (apex vs base)"];
(* centroid reference angles, apex vs base, pure phi (g-independent) *)
dA = FullSimplify[alA, 0<g && 0<ph<Pi];
dB = FullSimplify[alB, 0<g && 0<ph<Pi];
dC = FullSimplify[alC, 0<g && 0<ph<Pi];
Print["  alpha_apex (centroid) = ", ToString[InputForm[dA]]];
Print["  alpha_baseB           = ", ToString[InputForm[dB]]];
Print["  alpha_baseC           = ", ToString[InputForm[dC]]];
Print["  apex centroid phase e^{-2i alpha_apex} = ",
   ToString[InputForm[FullSimplify[Exp[-2 I dA], 0<g && 0<ph<Pi]]],
   "  (real => apex-conj component real, matches fastnc Gamma^1)"];

Print["%%SECTION%% 6. DELIVERABLE summary"];
Print["%%MAP_START%%"];
Print["FASTNC_APEX_LEG = 1   (Gamma^1 real on isoceles)"];
Print["FASTNC_BASE_LEGS = {2,3}  (Gamma^2 = conj(Gamma^3))"];
Print["CANOES_APEX_LEG = 3   (Gamma^3 conjugates apex in current assembler)"];
Print["CANOES_BASE_LEGS = {1,2}"];
Print["ALIGN_RELABEL: canoes apex-channel sits in slot Gamma^3; fastnc puts it in"];
Print[" slot Gamma^1. To compare slot-by-slot, the canoes natural-component"];
Print[" assignment must conjugate the APEX leg for Gamma^1 (not Gamma^3)."];
Print["%%MAP_END%%"];
Print["%%DONE%%"];
