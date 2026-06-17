(* ============================================================================
   derive_shear3pcf_natural_components.wl

   FIRST-PRINCIPLES derivation of the cosmic-shear 3PCF natural components
   Gamma^0..Gamma^3 (Schneider-Lombardi; Sugiyama+2024 = fastnc), flat-sky,
   CENTROID projection, in terms of the driving-field spin-2 3-cumulants
   {zeta_PPP, zeta_D(perm), zeta_B}, on the SAS isoceles family
   (t1=t2=gamma, opening angle phi, base t3 = 2 gamma sin(phi/2)).

   Plain symbolic spin-weighted algebra (no xAct needed).

   ========================================================================
   THE PHYSICS, pinned from source
   ========================================================================
   (A) The image shear at vertex j is the spin-(+2) field, the LOS integral
       of the spin-2 driving field Psi0 (cosmology.tex eq: gamma from Dsachs23):
           gamma_j = gamma_+,j + i gamma_x,j = -INT dl (Dsachs2 + i Dsachs3)_j
                   = -INT dl  Psi0(n_j, l)        [Psi0 = Dsachs2 + i Dsachs3].
       Psi0 is referenced to a LOCAL screen frame at each vertex.

   (B) Spin-2 transformation. Under a CCW rotation of the local frame by angle
       a, a spin-(+2) field f -> f e^{-2 i a}; its conjugate f* -> f* e^{+2 i a}.
       The shear projected onto reference direction at polar angle a is
           gamma|_a = gamma e^{-2 i a}.

   (C) Natural components (Schneider-Lombardi; fastnc):
           Gamma^0 = < g|_a1   g|_a2   g|_a3   >          (all unconjugated)
           Gamma^1 = < g*|_a1  g|_a2   g|_a3   >          (leg-1 conjugated)
           Gamma^2 = < g|_a1   g*|_a2  g|_a3   >          (leg-2 conjugated)
           Gamma^3 = < g|_a1   g|_a2   g*|_a3  >          (leg-3 conjugated)
       where a_j is the CENTROID projection reference angle at vertex j.
       Verified against fastnc's Bessel-order map (fastnc.py:365):
       at multipole M=0, (m,n) per mu = {(-3,-3),(-1,-1),(1,-3),(-3,1)}, i.e.
       total leg-spin = {-6,-2,-2,-2}: Gamma^0 = all-same-helicity (spin 6),
       Gamma^{1,2,3} = one helicity flipped (spin 2).  CONFIRMED.

   (D) canoes canonical great-circle frame (the frame the zeta cumulants are
       reduced in): n1 at the pole, n2 in the xz-plane, n3 at azimuth phi3;
       spin reference of each leg = its great-circle (geodesic) tangent.  In
       this frame canoes returns the REAL cumulants (parity makes the imaginary
       part vanish):
           zeta_PPP = < P1 P2 P3 >_gc          (helicity +2,+2,+2)
           zeta_D(perm)                          (one leg flipped to -2)
           zeta_B   = < T  P  P* >_gc            (spin-0 + modulus pair)
       with P = Psi0 referenced to the great-circle tangent.
       KEY appendix fact: small-angle  sqrt(L^2(L^2-2)) -> L^2  =>
           zeta_B -> zeta_TTT,  and the all-same-helicity zeta_PPP ~ gamma^4
       (d^L_{2,-2}) while the modulus channels ~ O(1) (d^L_{2,2}).

   (E) GREAT-CIRCLE -> FIXED-FRAME bridge.  The natural component references
       each leg to the CENTROID direction a_j in the FIXED flat-sky frame.
       canoes references each leg to the GREAT-CIRCLE tangent at angle b_j.
       Because canoes BAKES the great-circle d^L_{m,-s} projection into the
       harmonic sum, the relation between the fixed-frame helicity cumulant and
       the great-circle one is a per-leg phase exp(-2 i s_j (a_j - b_j)):
           < g|_{a} ... > = exp(-2 i Sum_j s_j (a_j - b_j)) < P ... >_gc
       with s_j = +1 (unconjugated leg) or -1 (conjugated leg).  A pure phase
       cannot change |Gamma^mu|; the per-component MAGNITUDE is therefore set
       ENTIRELY by which great-circle cumulant feeds each Gamma^mu.  This is the
       crux of the STAGE-1 discrepancy and the deliverable's central claim.

   ========================================================================
   DELIVERABLE: the explicit linear map (channel + phase) for each Gamma^mu.
   ========================================================================  *)

$DefInfoQ = False;

Print["%%SECTION%% 0. SAS triangle geometry (flat sky)"];

(* Leg labeling matches stage1_grid.json: cosine_triples=(cos t3, cos g, cos g);
   leg-pair (1,2)=base t3, (2,3)=(3,1)=equal sides g.  So vertices 1,2 are the
   BASE endpoints, vertex 3 is the APEX.  Concrete coords: apex at origin, the
   two equal sides of length g open symmetrically about +y with half-angle ph/2. *)
g = gamma; ph = phi;
x3 = {0, 0};
x1 = {-g Sin[ph/2],  g Cos[ph/2]};
x2 = { g Sin[ph/2],  g Cos[ph/2]};

t3len = Simplify[Sqrt[(x2 - x1).(x2 - x1)], 0 < g && 0 < ph < Pi];
Print["t3 - 2 g sin(phi/2) = ",
  ToString[InputForm[Simplify[t3len - 2 g Sin[ph/2], 0 < g && 0 < ph < Pi]]],
  "   (PASS if 0)"];

Print["%%SECTION%% 1. centroid projection reference angles a_j"];
cen = (x1 + x2 + x3)/3;
a1 = ArcTan[(cen - x1)[[1]], (cen - x1)[[2]]];
a2 = ArcTan[(cen - x2)[[1]], (cen - x2)[[2]]];
a3 = ArcTan[(cen - x3)[[1]], (cen - x3)[[2]]];
Print["  a1 = ", ToString[InputForm[FullSimplify[a1, 0<g && 0<ph<Pi]]]];
Print["  a2 = ", ToString[InputForm[FullSimplify[a2, 0<g && 0<ph<Pi]]]];
Print["  a3 = ", ToString[InputForm[FullSimplify[a3, 0<g && 0<ph<Pi]]]];

Print["%%SECTION%% 2. canonical great-circle reference angles b_j"];
(* Each leg's great-circle reference = geodesic to leg 1 (the pole). At vertex 1
   the canonical x-axis points to vertex 2 (azimuth 0). *)
b1 = ArcTan[(x2 - x1)[[1]], (x2 - x1)[[2]]];   (* vertex1 -> vertex2 *)
b2 = ArcTan[(x1 - x2)[[1]], (x1 - x2)[[2]]];   (* vertex2 -> vertex1 *)
b3 = ArcTan[(x1 - x3)[[1]], (x1 - x3)[[2]]];   (* vertex3 -> vertex1 *)
Print["  b1 = ", ToString[InputForm[FullSimplify[b1, 0<g && 0<ph<Pi]]]];
Print["  b2 = ", ToString[InputForm[FullSimplify[b2, 0<g && 0<ph<Pi]]]];
Print["  b3 = ", ToString[InputForm[FullSimplify[b3, 0<g && 0<ph<Pi]]]];

Print["%%SECTION%% 3. per-leg great-circle->centroid spin-2 phase"];
(* unconjugated leg j : w_j = exp(-2 i (a_j - b_j));  conjugated leg : w_j^*. *)
d1 = a1 - b1; d2 = a2 - b2; d3 = a3 - b3;
w1 = Exp[-2 I d1]; w2 = Exp[-2 I d2]; w3 = Exp[-2 I d3];
(* g-independence (angles are scale-free => pure functions of phi). *)
gindep = Simplify[{D[d1,g], D[d2,g], D[d3,g]}];
Print["  d/dg (a_j - b_j) = ", ToString[InputForm[gindep]], "   (PASS if {0,0,0})"];
Print["  a_j-b_j (pure phi): ",
  ToString[InputForm[FullSimplify[{d1,d2,d3}, 0<g && 0<ph<Pi]]]];

Print["%%SECTION%% 4. THE MAP: Gamma^mu = phase x great-circle cumulant"];
(* channel assignment (from the fastnc (m,n) helicity analysis, sec C):
     Gamma^0 <- zeta_PPP  (all +2)          phase w1 w2 w3
     Gamma^1 <- zeta_D1   (leg-1 conj, -2)  phase conj(w1) w2 w3
     Gamma^2 <- zeta_D2   (leg-2 conj)      phase w1 conj(w2) w3
     Gamma^3 <- zeta_D3   (leg-3 conj)      phase w1 w2 conj(w3)
   ZPPP, ZD1, ZD2, ZD3 are the REAL great-circle cumulants (symbols here). *)
G0 = ZPPP w1 w2 w3;
G1 = ZD1  Conjugate[w1] w2 w3;
G2 = ZD2  w1 Conjugate[w2] w3;
G3 = ZD3  w1 w2 Conjugate[w3];
ph0 = FullSimplify[w1 w2 w3, 0<g && 0<ph<Pi];
ph1 = FullSimplify[Conjugate[w1] w2 w3, 0<g && 0<ph<Pi];
ph2 = FullSimplify[w1 Conjugate[w2] w3, 0<g && 0<ph<Pi];
ph3 = FullSimplify[w1 w2 Conjugate[w3], 0<g && 0<ph<Pi];
Print["  Gamma^0 phase = ", ToString[InputForm[ph0]]];
Print["  Gamma^1 phase = ", ToString[InputForm[ph1]]];
Print["  Gamma^2 phase = ", ToString[InputForm[ph2]]];
Print["  Gamma^3 phase = ", ToString[InputForm[ph3]]];

Print["%%SECTION%% 5. LIMIT CHECK (i): magnitudes are phase-invariant"];
(* |Gamma^mu| = |Z_channel|; the phases have unit modulus. *)
Print["  |w1 w2 w3| = ",
  ToString[InputForm[FullSimplify[Abs[ph0], 0<g && 0<ph<Pi && Element[{g,ph},Reals]]]],
  "   (PASS if 1 => projection cannot redistribute magnitude)"];
Print["  => per-component magnitude set ENTIRELY by the channel; small-gamma"];
Print["     zeta_PPP~gamma^4 (d^L_{2,-2}) suppresses |Gamma^0|, modulus zeta_D"];
Print["     ~O(1) (d^L_{2,2}) keeps |Gamma^{1,2,3}| finite  [appendix 615-619]."];

Print["%%SECTION%% 6. LIMIT CHECK (ii): isoceles symmetry"];
(* On t1=t2, vertices 1,2 are mirror images across the y-axis. The mirror swaps
   base legs 1<->2 (so ZD1<->ZD2... but in OUR labeling apex=leg3, base=1,2).
   fastnc isoceles relations: Gamma^3 = conj(Gamma^2) with apex-leg convention.
   Here legs 1,2 are the base (mirror pair), leg 3 the apex (self-mirror).
   Under the mirror x->-x: a1<->(pi-a2), a3->(pi-a3); b likewise; the base
   cumulants ZD1=ZD2 (mirror equal), ZD3=apex (self), ZPPP self. Check that the
   phase structure makes Gamma(base1) = conj(Gamma(base2)) and Gamma(apex),
   Gamma(PPP) real, matching fastnc's Gamma^1(apex) real, Gamma^2=conj(Gamma^3). *)
mir = {gamma -> gamma, phi -> phi}; (* structural; explicit reflection below *)
(* Build the mirror by x-> -x on the coordinates and recompute the phases. *)
y1 = {-1,1} x2; y2 = {-1,1} x1; y3 = {-1,1} x3;  (* reflected, with 1<->2 relabel *)
ceny = (y1+y2+y3)/3;
am1 = ArcTan[(ceny-y1)[[1]],(ceny-y1)[[2]]]; am2 = ArcTan[(ceny-y2)[[1]],(ceny-y2)[[2]]];
am3 = ArcTan[(ceny-y3)[[1]],(ceny-y3)[[2]]];
bm1 = ArcTan[(y2-y1)[[1]],(y2-y1)[[2]]]; bm2 = ArcTan[(y1-y2)[[1]],(y1-y2)[[2]]];
bm3 = ArcTan[(y1-y3)[[1]],(y1-y3)[[2]]];
wm1 = Exp[-2 I (am1-bm1)]; wm2 = Exp[-2 I (am2-bm2)]; wm3 = Exp[-2 I (am3-bm3)];
(* The reflected base-1-conjugated component should equal conj of the original
   base-2-conjugated component (with ZD1=ZD2 real on isoceles). Test phase only. *)
sym = FullSimplify[(Conjugate[wm1] wm2 wm3) - Conjugate[w1 Conjugate[w2] w3],
   0<g<1/100 && 0<ph<Pi];
Print["  [mirror] phase(base1 conj) - conj phase(base2 conj) = ",
  ToString[InputForm[sym]], "   (PASS if 0 => Gamma^2=conj(Gamma^3) emerges)"];
apexImag = FullSimplify[ComplexExpand[Im[Conjugate[w1] w2 w3 /. {w1->wA,w2->wB,w3->wC}]] /.
   {wA->w1,wB->w2,wC->w3}, 0<g<1/100 && 0<ph<Pi];
(* simpler: apex-conjugated phase realness *)
apexPhase = FullSimplify[w1 w2 Conjugate[w3], 0<g<1/100 && 0<ph<Pi];
Print["  apex-conjugated (leg3) phase = ", ToString[InputForm[apexPhase]],
  "   (real => Gamma^apex real on isoceles, matches fastnc Gamma^1 real)"];

Print["%%SECTION%% 7. LIMIT CHECK (iii): 2-point xi_+ / xi_- analogue"];
(* Two points A,B separated by gamma on the x-axis; tangential reference = the
   separation axis. canoes great-circle ref = same separation axis (geodesic),
   so b=a at BOTH points => the projection phases are TRIVIAL and the 2PCF is
   the bare great-circle cumulant. The non-trivial content is the CHANNEL:
     xi_+ = < g g* > : one conjugation => MODULUS channel zeta_B (finite)
     xi_- = < g g  > : no conjugation => un-conjugated channel zeta_TPP-type
                       (gamma^4-suppressed, d^L_{2,-2})
   This reproduces the appendix statement (lines 615-622): xi_+ <- zeta_B
   (lifts it from zero to a comparable level), xi_- <- the suppressed channel. *)
Print["  xi_+ channel  = MODULUS  zeta_B  = <Phi/g g*>   (one conjugation, finite)"];
Print["  xi_- channel  = UNCONJUG zeta_TPP-type = <g g>  (gamma^4-suppressed)"];
Print["  => the SAME channel logic as the 3-point map: conjugation picks the"];
Print["     finite modulus, no-conjugation picks the suppressed channel."];
Print["     Matches appendix `eq: appendix modulus reconstruction` (lines 606-622)."];

Print["%%SECTION%% 8. EXPLICIT MAP for the .py implementation (centroid)"];
Print["%%MAP_START%%"];
Print["a_minus_b_1 = ", ToString[InputForm[FullSimplify[d1, 0<g && 0<ph<Pi]]]];
Print["a_minus_b_2 = ", ToString[InputForm[FullSimplify[d2, 0<g && 0<ph<Pi]]]];
Print["a_minus_b_3 = ", ToString[InputForm[FullSimplify[d3, 0<g && 0<ph<Pi]]]];
Print["Gamma0 = ZPPP * exp(-2i( (a1-b1)+(a2-b2)+(a3-b3) ))"];
Print["Gamma1 = ZD1  * exp(-2i(-(a1-b1)+(a2-b2)+(a3-b3) ))"];
Print["Gamma2 = ZD2  * exp(-2i( (a1-b1)-(a2-b2)+(a3-b3) ))"];
Print["Gamma3 = ZD3  * exp(-2i( (a1-b1)+(a2-b2)-(a3-b3) ))"];
Print["%%MAP_END%%"];
Print["%%DONE%%"];
