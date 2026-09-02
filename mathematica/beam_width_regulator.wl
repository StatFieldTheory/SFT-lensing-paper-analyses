(* ::Package:: *)
(* beam_width_regulator.wl

   Can a finite beam width supply the ultraviolet regulator that the FK
   term needs?

   WHY THIS IS ASKED. The FK contribution samples the driving-field
   three-point cumulant in the COLLAPSED configuration: two of the three
   legs sit at the same point on the ray. A cumulant with coincident legs
   is a variance-like object, and after the Limber radial collapse the hard
   pair contributes as a two-dimensional transverse integral
   INT dk k P(k), whose integrand goes as k^(1+n_eff). That converges only
   once n_eff < -2. The linear spectrum reaches -2 at k = 0.21 h/Mpc, so
   the tree-level FK term has a finite limit. A halofit spectrum, flattened
   by its one-halo term, only reaches -2 at k = 8.4 h/Mpc and sits at
   n_eff = -1.15 at k = 1, i.e. essentially AT the peak of the integrand,
   so every decade contributes comparably and the integral does not
   converge on any scale the model controls.

   The paper's formalism carries no smoothing scale: a search of sections/
   finds no beam width, no coarse-graining, and ell_max never appears in
   the text at all. So the quantity being computed is a coincident-point
   moment of a mathematically point-like beam, and that quantity has no
   finite value for a nonlinear input.

   The physical candidate for a regulator is that a real beam has finite
   cross-section, so the driving fields entering the Sachs system are
   averages over that cross-section rather than point values.

   WHAT THIS SCRIPT ESTABLISHES, and what it does not.

   It establishes the ALGEBRA of screen-plane averaging: the correction to
   the point value, the Fourier-space window, the fact that an axisymmetric
   window does not mix spin weights, and that the windowed transverse
   integral converges for any power-law spectrum. It also gives the
   numerical map between a beam width and an effective ell_max.

   It does NOT establish that beam averaging is the physically correct
   definition of the observable, nor which scale to use. That is a choice
   about what is being measured, and this script is deliberately silent on
   it. What it does show is that the choice, once made, removes the
   divergence and fixes ell_max rather than leaving it free.

   Conventions follow derive_kappa3_spin2_helicity.wl:
     P  = sigma_+ + I sigma_x    (helicity +2)
     Phi = Phi00                 (spin 0, real)
   Screen indices A, B run over the two-dimensional screen plane.
*)

ClearAll["Global`*"];
pass = {}; fail = {};
check[name_, cond_] := If[TrueQ[cond], AppendTo[pass, name], AppendTo[fail, name]];

Print["============================================================"];
Print["beam_width_regulator.wl"];
Print["============================================================"];

(* ---------------------------------------------------------------- *)
(* T1. Screen-plane average of a scalar driver                       *)
(* ---------------------------------------------------------------- *)
(* A normalised, axisymmetric window W on the screen plane has
     INT W d^2 xi = 1,  INT W xi^A d^2 xi = 0,
     INT W xi^A xi^B d^2 xi = (sb^2/2) delta^{AB},
   the last line being the definition of sb. Taylor expanding the driver
   about the central ray and averaging term by term gives the correction
   below. Everything here is exact in the expansion, not a model. *)

Clear[sb, k];
(* second moment of a normalised 2D Gaussian of width sb: <xi^A xi^B> = (sb^2/2) delta *)
gauss2D[xi1_, xi2_, sb_] := Exp[-(xi1^2 + xi2^2)/sb^2]/(Pi sb^2);
norm  = Integrate[gauss2D[x1, x2, sb], {x1, -Infinity, Infinity},
                  {x2, -Infinity, Infinity}, Assumptions -> sb > 0];
mom2  = Integrate[x1^2 gauss2D[x1, x2, sb], {x1, -Infinity, Infinity},
                  {x2, -Infinity, Infinity}, Assumptions -> sb > 0];
momXY = Integrate[x1 x2 gauss2D[x1, x2, sb], {x1, -Infinity, Infinity},
                  {x2, -Infinity, Infinity}, Assumptions -> sb > 0];
check["T1a window is normalised", Simplify[norm == 1]];
check["T1b second moment is sb^2/2", Simplify[mom2 == sb^2/2]];
check["T1c window is isotropic (no xy moment)", Simplify[momXY == 0]];
Print["T1  normalisation ", norm, ";  <xi_1^2> = ", Simplify[mom2],
      ";  <xi_1 xi_2> = ", momXY];

(* Averaging the Taylor series: Phibar = Phi + (1/2) <xi^A xi^B> d_A d_B Phi + ...
   With <xi^A xi^B> = (sb^2/2) delta^{AB} this is Phi + (sb^2/4) Lap Phi. *)
coeffLaplacian = (1/2) (sb^2/2);
check["T1d Laplacian coefficient is sb^2/4", Simplify[coeffLaplacian == sb^2/4]];
Print["T1  Phibar = Phi + (", Simplify[coeffLaplacian], ") Lap_perp Phi + O(sb^4)"];

(* ---------------------------------------------------------------- *)
(* T2. The same statement in Fourier space                           *)
(* ---------------------------------------------------------------- *)
(* Averaging is a convolution, so in Fourier it is multiplication by the
   window transform. For the Gaussian above the transform is exp(-k^2 sb^2/4),
   whose small-k expansion must reproduce T1 with Lap_perp -> -k^2. *)
Wtilde[k_, sb_] := Integrate[gauss2D[x1, x2, sb] Exp[I k x1],
                             {x1, -Infinity, Infinity}, {x2, -Infinity, Infinity},
                             Assumptions -> sb > 0 && k \[Element] Reals];
wt = Simplify[Wtilde[k, sb], Assumptions -> sb > 0];
check["T2a window transform is exp(-k^2 sb^2/4)", Simplify[wt == Exp[-k^2 sb^2/4]]];
smallK = Normal[Series[wt, {k, 0, 2}]];
check["T2b small-k expansion matches the Laplacian correction",
      Simplify[smallK == 1 - (sb^2/4) k^2]];
Print["T2  Wtilde(k) = ", wt, ";  expansion ", smallK,
      "  (Lap_perp -> -k^2 reproduces T1)"];

(* ---------------------------------------------------------------- *)
(* T3. An axisymmetric window does not mix spin weights              *)
(* ---------------------------------------------------------------- *)
(* A spin-s field on the screen picks up e^{I s psi} under a rotation of the
   screen basis by psi. Averaging with a window that depends only on |xi|
   commutes with that rotation, so the averaged field carries the same spin
   weight. The check below is the concrete statement: the angular integral
   of e^{I s theta} against an axisymmetric window vanishes unless s = 0,
   so no other spin weight is generated. *)
mixing[s_] := Integrate[Exp[I s th], {th, 0, 2 Pi}]/(2 Pi);
check["T3a spin 0 survives averaging", Simplify[mixing[0] == 1]];
check["T3b spin 2 is not generated from spin 0", Simplify[mixing[2] == 0]];
check["T3c spin -2 is not generated from spin 0", Simplify[mixing[-2] == 0]];
check["T3d spin 4 is not generated", Simplify[mixing[4] == 0]];
Print["T3  angular projection of an axisymmetric window: s=0 -> ", mixing[0],
      ", s=2 -> ", mixing[2], ", s=-2 -> ", mixing[-2]];
Print["T3  so Phi00 stays spin 0 and Psi0 stays spin 2 under beam averaging."];

(* ---------------------------------------------------------------- *)
(* T4. The windowed transverse integral converges for any slope      *)
(* ---------------------------------------------------------------- *)
(* Unwindowed, the hard pair gives INT dk k P(k) with P ~ k^n, which
   diverges at large k unless n < -2. With the window the integrand carries
   Wtilde^2 = exp(-k^2 sb^2/2), which decays faster than any power, so the
   integral converges for EVERY n. That is the whole content of the
   proposal: the divergence is not removed by a better spectrum, it is
   removed by admitting that the beam has finite width. *)
Clear[n];
unwindowed[n_] := Integrate[k^(1 + n), {k, 1, Infinity}];
windowed[n_, sb_] := Integrate[k^(1 + n) Exp[-k^2 sb^2/2], {k, 1, Infinity},
                               Assumptions -> sb > 0];
divergesTree   = unwindowed[-1.15];   (* halofit slope at k = 1 h/Mpc *)
convergesTree  = unwindowed[-2.29];   (* linear slope at k = 1 h/Mpc *)
check["T4a unwindowed integral diverges at the halofit slope",
      Head[divergesTree] =!= Real || divergesTree === Infinity ||
      !FreeQ[divergesTree, Infinity] || !NumericQ[divergesTree]];
check["T4b unwindowed integral converges at the linear slope",
      NumericQ[convergesTree] && convergesTree > 0];
wHalofit = windowed[-1.15, 1];
check["T4c windowed integral converges at the halofit slope",
      NumericQ[N[wHalofit]] && N[wHalofit] > 0];
Print["T4  unwindowed, n = -2.29 (linear at k=1):  ", convergesTree];
Print["T4  unwindowed, n = -1.15 (halofit at k=1): ", divergesTree];
Print["T4  windowed  , n = -1.15, sb = 1:          ", N[wHalofit]];

(* ---------------------------------------------------------------- *)
(* T5. What beam width corresponds to what cutoff                    *)
(* ---------------------------------------------------------------- *)
(* The window suppresses transverse modes above k ~ 2/sb. With the Limber
   map k = ell/chi_h this is an effective multipole cutoff
   ell_eff = 2 chi_h / sb. Evaluated at the shells the table uses. *)
chiH = {292.6, 2291., 4445., 5582.2};        (* Mpc/h, inner to source *)
shellName = {"innermost", "mid", "outer", "source"};
ellEff[sbMpc_, chi_] := 2 chi/sbMpc;
Print["T5  effective cutoff ell = 2 chi_h / sb, for a beam width sb in Mpc/h:"];
Do[
  Print["      sb = ", sb0, " Mpc/h -> ell_eff = ",
        Table[Round[ellEff[sb0, c]], {c, chiH}], "  (", shellName, ")"],
  {sb0, {0.1, 0.5, 1.0, 3.0}}];
check["T5a a beam of order 1 Mpc/h gives ell_eff of a few thousand at the source shell",
      2000 < ellEff[1.0, 5582.2] < 20000];

(* ---------------------------------------------------------------- *)
(* T6. What averaging does NOT do                                    *)
(* ---------------------------------------------------------------- *)
(* The Sachs system is nonlinear in the optical scalars: the expansion and
   shear obey equations containing rho^2 and sigma sigma-bar. Averaging a
   nonlinear equation over the bundle does not commute with the
   nonlinearity, since <rho^2> is not <rho>^2. So "average the drivers" is
   NOT the same operation as "solve for a finite beam", and the difference
   is second order in the optical scalars. The check below is the algebraic
   statement of that non-commutation, recorded so the limitation is not
   forgotten. *)
Clear[r, dr];
avgOfSquare = Expand[(r + dr)^2];
squareOfAvg = r^2;
mismatch = Simplify[avgOfSquare - squareOfAvg /. dr^2 -> var];
check["T6a averaging does not commute with the quadratic term",
      Simplify[(mismatch /. {dr -> 0, var -> 0}) == 0] &&
      Simplify[mismatch =!= 0]];
Print["T6  <(rho + drho)^2> - <rho>^2 = ", mismatch /. dr -> 0,
      "   (nonzero: the bundle variance survives)"];
Print["T6  so beam averaging defines the DRIVER, it does not solve the"];
Print["T6  finite-beam problem; the two differ at second order in the"];
Print["T6  optical scalars."];

(* ---------------------------------------------------------------- *)
Print["============================================================"];
Print["PASS (", Length[pass], "): ", pass];
If[Length[fail] > 0, Print["FAIL (", Length[fail], "): ", fail],
   Print["FAIL (0): none"]];
Print["============================================================"];
