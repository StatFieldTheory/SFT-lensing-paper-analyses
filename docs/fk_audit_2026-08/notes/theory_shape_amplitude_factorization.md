# Why the collapsed FK vertex factorizes into shape times amplitude

This note derives, from the squeezed limit of the tree-level bispectrum,
the empirical facts established in `finding_three_effect_decomposition.md`:
the multipole cutoff predominantly rescales the FK amplitude while the
angular shape is carried by a soft factor the cutoff barely touches, and
the corrected FK curve tracks Order-0 at a nearly constant ratio. Each
step is a lemma with its proof idea, validity domain and verification
status.

Verification status. Three independent layers:

* symbolic: every expansion coefficient is machine-verified by
  `code/verify_squeezed_expansion.wl` (7 checks, all PASS, including
  remainder-scaling assertions and the collapsed-family kernel reduction);
* independent algebra: a review pass re-derived Lemmas 1 and 2, the
  averages and the projection-slice identity in sympy with exact
  rationals, confirming every coefficient and locating the exact origin
  of the n-dependence (the dipole times the O(eps) difference of the two
  hard spectra);
* direct quadrature: `code/verify_factorization_numeric.py` integrates
  the raw vertex against the real linear P(k) and confirms the
  factorization at the percent level inside its validity window, along
  with the sign-change and coefficient predictions (section 7).

Statements about pipeline measurements cite the finding notes. Section 7
separates what was derived, what was predicted and then verified, and
what is explained after the fact.

Redshift is held fixed throughout this note (one shell in sections 0-4,
one source in section 5). The companion note
`theory_redshift_dependence.md` makes both explicit: the per-shell
growth-times-geometry factorization (exact at tree level, certified on
the real table), the shell-redshift dependence of THIS note's
factorization accuracy (percent-level for z >= 2.5 shells, order-one at
z = 0.1), and the measured source-redshift trends of the fold.

Notation: 2D multipole vectors l with magnitude l; chi the shell comoving
distance; gamma the external separation; P(k) the linear matter power at
the shell redshift (suffix NL when nonlinear); n(k) = dln P/dln k.
Sections 0-4 work at one shell; the fold enters in section 5.

## 0. The object

The Limber (HIGH) branch of the collapsed vertex, as implemented and
source-verified in `input_ready_b_to_zeta_spec.md` (section 2), is for
the TTT channel on the collapsed family (1, cos gamma, cos gamma):

    zeta(gamma; chi) = C(chi) INT_0^inf du u INT_0^inf dv v INT_0^2pi dphi
                       W(u, v, w) B(w/chi, u/chi, v/chi) 2pi J_0(v gamma)   (0.1)

with w = sqrt(u^2 + v^2 + 2 u v cos phi), the hard window
W = 1 iff ell_cut < max(u, v, w) <= ell_max, and C(chi) collecting the
per-leg response scalars, (2pi)^-4, chi^-4 and the measure factors, all
independent of gamma and ell_max.

The 2pi J_0(v gamma) factor follows from the spec's general kernel
alpha_TTT = 2pi J_0(|Q|), Q = u r2 + v e^{-i phi} r3, together with the
flat-triangle map (KAPPA3:1732-1744): for the collapsed triple
(1, c, c) the map gives theta_12 = 0 hence r2 = 0, and |r3| =
theta_31 = gamma, so |Q| = v gamma identically in (u, phi). This
reduction is machine-checked in the .wl (Collapsed-family |Q|, PASS).

Equation (0.1) is the flat-sky statement

    zeta ~ < delta^2(x_a) delta(x_b) >,  |x_b - x_a| = chi gamma           (0.2)

for the Limber-projected field: with delta(x) = INT d^2l/(2pi)^2
delta~(l) e^{i l.x} and the momentum delta function l1 + l2 + l3 = 0, the
phase e^{i(l1+l2).x_a} cancels and only e^{i l3.(x_b-x_a)} survives; the
exact integral over the absolute orientation contributes the explicit
factor 2pi J_0(l3 chi gamma / chi) = 2pi J_0(v gamma) on the single leg
that connects the separated point (independently re-derived in the review
pass). Two structural facts follow:

* F1. gamma enters (0.1) only through J_0(v gamma), attached to the
  single leg v. Convergence bookkeeping, stated carefully because it
  matters below: J_0(x) has envelope x^{-1/2}, so INT dv v J_0(v gamma)
  f(v) converges (as an improper integral) iff the spectrum factor f
  falls faster than v^{-1/2}; and because the surviving oscillation of
  the partial integral INT^L decays only like an oscillatory v^{-1}
  envelope, the CAPPED integral approaches its uncapped value only when
  gamma L >> 10 or so. Both facets are used in section 3.
* F2. ell_max enters only through the window on max(u, v, w).

F1 and F2 alone do not give the separation: if the integral were
dominated by u ~ v ~ w (all hard), raising ell_max would add weight at
large v too, and J_0(v gamma) would turn that into a gamma-dependent
change. The separation needs a property of B: sections 1 and 2.

## 1. Lemma 1: the hard-pair coupling is doubly soft

Let q be soft against k, eps = q/k, mu = cos(angle between k and q),
vectors coplanar as appropriate for the Limber integral (the result
depends only on magnitudes and mu, so it is the same in 2D and 3D). With
F2(k1,k2) = 5/7 + (mu12/2)(k1/k2 + k2/k1) + (2/7) mu12^2,

    F2(k, q - k) = eps^2 [ 3/14 - (5/7) mu^2 ] + O(eps^3).               (1.1)

Proof: direct expansion; the O(1) and O(eps) pieces cancel identically.
Verified symbolically (PASS), independently in sympy (O(eps^3)
coefficient mu(13 - 20 mu^2)/14, nonzero, so the remainder order is
sharp), and the remainder-scaling assertion res(eps)/res(eps/2) ~ 2^3
PASSes. Physical origin: mass and momentum conservation force the
second-order density at total wavenumber q -> 0 to vanish as q^2.

Consequence. In B_tree = 2F2(k1,k2)P1P2 + 2F2(k2,k3)P2P3 +
2F2(k3,k1)P3P1 with k3 = k_s soft and (k1, k2) the hard pair, the
hard-hard term is O(eps^2) P_h^2, against soft-hard terms of order
P_s P_h. Their ratio is O(eps^2) P_h/P_s = O(eps^{2-n}), at most
O(eps^2) for the relevant n < 0. This kills the one pairing that would
have coupled the cutoff to the gamma dependence (both P-carrying legs
hard under the window while the soft spectator carries J_0).

## 2. Lemma 2: the soft-hard terms sum to a response times P_s P_h

Same configuration; c = cos(angle between the soft leg and hard leg k1),
n = n(k_h), P locally power-law at the hard scale. Each soft-hard F2
diverges as 1/eps (dipole), but the two permutations carry exactly
opposite poles (+c and -c: sympy), and the O(eps) difference between the
two hard spectra P(|k1|) vs P(|k1 + k_s|) converts the dipole into the
finite n-dependent piece:

    2 F2(k1,k_s) P_1 P_s + 2 F2(k2,k_s) P_2 P_s
        = 2 P(k_s) P(k_h) [ 13/14 + (4/7) c^2 - (n/2) c^2 ] + O(eps).    (2.1)

Verified symbolically (PASS; remainder scaling ~ 2^1 PASS) and in sympy;
the review pass isolated the mechanism by freezing P_2 -> P(k_h), which
removes exactly the -(n/2) c^2 term. Anchors:

* 3D average <c^2> = 1/3: coefficient 47/21 - n/3, the textbook
  squeezed limit of tree SPT and the tree-level separate-universe
  response (PASS).
* Coplanar (Limber) average over the single angle phi, <c^2>_phi = 1/2:

      Rbar(k_h) = 17/7 - n(k_h)/2 .                                      (2.2)

  (PASS.) This 2D coefficient is what the vertex integral uses; the
  n-independent part alone differs from the 3D value by 8.5%, so
  substituting the 3D coefficient would be a ten-percent-level error.

So in the squeezed regime

    B_tree = R(c; k_h) P(k_s) P(k_h) [1 + O(eps) + O(eps^{2-n})],        (2.3)

with R an O(1) kernel whose phi-average is (2.2).

## 3. Proposition: factorization of the collapsed vertex

Insert (2.3) into (0.1) with k_s = v/chi (the J_0 leg), k_h = u/chi ~
w/chi. In the region v << u the phi integral acts only on the
c-dependence of R (J_0 is phi-independent; the window's phi-dependence is
O(v/u)), and INT dphi R = 2pi Rbar. The double integral separates:

    zeta(gamma; chi) = C~(chi) S_cap(gamma; chi; ell_max) H(chi; ell_max)
                       x [ 1 + E ],      C~ = (2pi)^2 C,                  (3.1)

    S_cap(gamma; chi; L) = INT_0^L dv v J_0(v gamma) P(v/chi)            (3.2)
    H(chi; L)            = INT_{ell_cut}^L du u Rbar(u/chi) P(u/chi).    (3.3)

The (2pi)^2 is one factor from the explicit 2pi J_0 and one from
INT dphi Rbar; the first draft of this note omitted the second and the
direct quadrature caught it. The soft factor carries the window cap: by
F1 the capped integral reaches the uncapped one only when
gamma L >> ~10, so writing INT_0^inf is only legitimate in that limit.
In that limit the soft factor has a name: by the projection-slice
identity

    (1/2pi) INT k dk J_0(k r) P(k) = INT_{-inf}^{inf} dDelta
                                     xi_L( sqrt(r^2 + Delta^2) )          (3.4)

(one line: insert P = FT of xi and integrate the line-of-sight
wavenumber against INT dDelta e^{i k_par Delta} = 2pi delta(k_par);
verified analytically and numerically to 1e-11 in the review pass),

    S_inf(gamma; chi) = 2 pi chi^2 Sigma(chi gamma),
    Sigma(r_perp) = INT dDelta xi_L( sqrt(r_perp^2 + Delta^2) ),          (3.5)

the transverse projected linear correlation, the w_p(r_p) of galaxy
surveys. H is the response-weighted transverse variance of the linear
field between the cuts. Shape lives in S, cutoff dependence in H and in
the cap of S; the separation is clean exactly to the extent the cap of S
is saturated.

Validity domain and error budget E:

* (i) Squeezedness needs the hard integral dominated by u >> 1/gamma.
  The per-log hard integrand u^2 Rbar P rises while n > -2, i.e. up to
  k ~ 0.5 h/Mpc (`finding_ell_max_resolved.md`), so the hard support
  sits at u ~ min(ell_max, ell_to) with ell_to ~ 0.5 chi_h (about 2200
  at the mid shell). Squeezedness requires gamma >> 1/min(ell_max,
  ell_to); tightening the window by raising the cutoff stops at ell_to.
* (ii) The O(eps) remainder of (2.1) plus the O(v/u) window and P(w)
  expansions give E = O(<v>/<u>) = O(1/(gamma ell_*)), ell_* the scale
  dominating H.
* (iii) The hard-hard term: E += O(eps^{2-n}) by Lemma 1.
* (iv) Soft-cap saturation: replacing S_cap by S_inf is only valid for
  gamma ell_max >> ~10; below that the truncation oscillation of the
  J_0 tail (v^{-1} envelope) makes S itself cutoff-dependent at the
  tens-of-percent level. This is a REAL, quantified effect: at
  ell_max = 1000 the direct quadrature finds 26-44% shape shifts
  between cutoffs at gamma = 5-40 arcmin, reproduced to 4-5% by the
  S_cap ratio alone (section 7).
* (v) The swapped region (u below ell_cut, v hard in the window) is an
  additive contribution not of the form S x H, roughly
  [INT_0^{ell_cut} u du P] x [INT_window v dv Rbar P J_0(v gamma)], in
  which gamma and ell_max sit in the same integral. It is suppressed by
  the small transverse variance below ell_cut relative to H and decays
  with gamma through its own Bessel, but it is a candidate source of
  residual cutoff-ratio drift inside the window and is not bounded by
  (i)-(iv).
* (vi) The exact-3j LOW branch (all legs <= ell_cut) is outside this
  derivation. It is ell_max-independent by construction, so it adds no
  cutoff dependence of its own; its admixture matters where it is
  comparable (section 6).

## 4. Corollary A: what the cutoff can and cannot do

Within the validity domain, d/d ell_max acts on H and on the cap of S.
Where the cap is saturated (gamma ell_max >> 10), the shape

    zeta(gamma) / zeta(gamma_0)   is independent of ell_max,             (4.1)

and the amplitude carries the cumulative hard integral. Precisely: the
octave band dH/dln L = L^2 Rbar(L/chi) P(L/chi) has local log-slope
2 + n (up to the slowly varying Rbar); d ln H/d ln L = band/H approaches
2 + n only while H is top-band dominated (n > -2, L well above ell_cut)
and rolls to zero once n < -2. The band statement is exactly the
octave-band scaling measured in `finding_ell_max_resolved.md`, including
the n = -2 turnover criterion.

Response reading: (3.1) is the position-dependent power spectrum
statement <delta_S^2(x_a) delta(x_b)> = response x sigma_S^2(cuts) x
xi_L(separation), the separate-universe picture of the squeezed
bispectrum. The vertex correlates local small-scale power with the
large-scale field; the small-scale variance carries the cutoff, the
large-scale correlation the separation.

Nonlinear bispectrum: the consistency structure survives beyond tree
level with P(k_h) -> the nonlinear response dP_NL/ddelta_L, the soft leg
still linear. Hence the explanation of the measured behaviour (section
7; the x1.66 measurement predates this note): a nonlinear B moves H, and
the shape only within the same E budget, i.e. at the band-drift level,
not exactly.

## 5. Corollary B: the fold, and why FK tracks Order-0

The FK fold multiplies the vertex by response factors that carry no
angular dependence (the FK diagram has four R lines and no C lines; R
depends on times only; spec section 0). By the same per-leg Limber
cancellation that the spec verifies for the three-point kernel, applied
to the two-point kernel, the Phi00 two-point spectrum per shell is
proportional to P(l/chi), so both curves fold the same soft transform:

    xi_O0(gamma) = INT dchi W_2(chi) S(gamma; chi),
    xi_FK(gamma) = INT dchi W_3(chi) H(chi) S(gamma; chi) [1 + E],       (5.1)

with W_2, W_3 the gamma-independent fold weights of the two diagrams.
Hence

    xi_FK / xi_O0 (gamma) = INT W_3 H S / INT W_2 S ,                    (5.2)

a weighted integral of H whose gamma dependence enters only through the
chi-profile of S(gamma; .): the per-shell gamma dependence cancels. The
ratio is exactly constant iff W_3 H is proportional to W_2. Its drift is
controlled by the MISMATCH between the two weight profiles, not merely
by the spread of H: even a constant H drifts if W_3 and W_2 weight the
shells differently (review counterexample), and near a zero of S the
ratio is ill-conditioned. What the structure does guarantee is that the
vertex contributes no angular feature of its own: at tree level FK is a
nearly multiplicative dressing of the Gaussian curve, of relative size
~ Rbar x (transverse hard variance), and a crossover with Order-0 cannot
be generated by the vertex within this regime.

Caveat on the measured comparison: the production Order-0 is the exact
non-Limber corr_op propagator, not the Limber form used in (5.1), so the
measured FK/O0 carries an additional Limber-projection error in the
denominator, largest at wide angles.

## 6. Where the clean picture degrades

Four regimes, all visible in the measurements:

* gamma <~ few/ell_max: not squeezed; shape responds to the cutoff.
  Measured (`finding_three_effect_decomposition.md`, conv8192/cut1000):
  4.06 at 0.5', 2.32 at 3.3', 2.04 at 21.9', 1.72 at 56.3', 1.90 at
  144.7', i.e. a strong small-gamma excess and a non-monotonic tail.
* gamma ell_max <~ 10: the soft cap is unsaturated and S itself is
  cutoff-dependent (error item (iv)); at ell_max = 1000 this covers much
  of the arcminute-to-degree range and quantitatively accounts for the
  drift between the small-gamma excess and the asymptotic hard-factor
  ratio H(4000)/H(1000) = 2.13.
* Near the zero of S: Sigma changes sign at a few hundred Mpc; ratios at
  a zero crossing are ill-conditioned, and any reweighting (different
  H(chi), different cutoff, the truncation oscillation of S_cap) moves
  the folded zero. Large local swings there carry no information.
* Wide angles where the fully-soft LOW component is comparable: the
  archived band table (`coverage_ell_z_mapping.md`) puts the high-ell
  share of zeta_TTT at the far shell at 40% (42'), 26% (117'), 6%
  (326'); sharper still, `finding_ell_cut_seam.md` measures LOW = 2.9x
  the TOTAL at 42' on the far shell, with HIGH cancelling it at nearly
  equal magnitude. LOW is cutoff-independent, HIGH grows with ell_max,
  so their mix, and hence the total shape, shifts with the cutoff
  wherever both matter.

## 7. Confrontation with the measurements, honestly labeled

Derived here and independently measured earlier (agreement, but the
explanation came after the data):

* the cutoff acts dominantly on the amplitude, with growth given by the
  cumulative hard integral and the n = -2 turnover;
* the nonlinear swap moves the amplitude (x1.66) and leaves the shape
  within the band tolerance;
* near-constant FK/O0 within the 1-2 percent band with a factor ~2
  drift.

Predicted by the first draft as falsifiable, then verified during the
review round by direct quadrature of the raw integrand against the real
P(k) (`code/verify_factorization_numeric.py`; independent review run at
chi_h = 4445 Mpc/h, grid-convergence checked):

* factorization: raw / [(2pi)^2 S_cap H] = 1.00 +/- 0.06 over
  gamma = 0.5-20 arcmin at ell_max = 1000, degrading exactly where the
  error budget says (near the zero of S; at the smallest gamma for the
  larger cutoff, where raw(4000)/raw(1000) = 5.2-5.5 at 0.5-1 arcmin
  against the hard-factor ratio 2.13, the ~30% excess being the
  non-squeezed contribution);
* the sign change of the vertex lands on the zero of the soft factor:
  raw crosses zero at 79.8-80.1 arcmin, the uncapped Sigma at 83.4
  arcmin (4%, the size of the stated O(1/(gamma ell_*)) error;
  sub-percent coincidence at ell_max = 4000). The capped S_cap's own
  zero is displaced to ~61 arcmin by its truncation oscillation,
  illustrating the ill-conditioning of the zero neighbourhood rather
  than contradicting the identification;
* the coefficient (2.2): the fitted constant recovers 17/7 - n/2 to
  ~1% when evaluated as an H-weighted average over the hard support
  (3.08-3.14 against 3.090); evaluating it at the single peak scale of
  u^2 P overestimates by 11% because n(k) runs, which quantifies the
  locally-power-law caveat of section 8.

Explained in this language, previously established by measurement:

* the deployed-grid artifact is the replacement S(gamma) -> S(0): the
  interpolation returns the zero-separation value of the soft factor,
  deleting the only gamma-dependent factor in (3.1)
  (`finding_grid_artifact_measured.md`);
* the MC cross-check is blind to all of this because both sides share
  the table; the fold multiplies the vertex by gamma-independent
  factors, preserving (3.1) (`finding_mc_crosscheck.md`).

Scope of the confrontation: the measured kk curves mix TTT with the
spin-2 channels (K011/K022 enter the kk reconstruction), so the fold
comparisons above test the structural factorization, not the
channel-specific coefficients (section 8).

## 8. What this note does not claim

* Nothing about the absolute amplitude: H is underconverged at the
  deployed cutoff, and its converged value is dominated by k of several
  h/Mpc where the bispectrum model is least reliable
  (`finding_ell_max_resolved.md`).
* The spin-2 channels are not derived. Their collapsed-family kernels
  replace J_0(v gamma) by J_|S|(v gamma) with spin phases e^{i s3 phi}
  coupling to the phi-moments of R(c); the factorization logic is
  structurally unchanged, but the channel analogues of (2.2) are
  unverified.
* The equal-shell collapse, measure bookkeeping and Limber cancellation
  are taken from the verified spec, not re-derived; the two-point
  Limber form used in section 5 is the same cancellation applied to the
  two-point kernel, and the production Order-0 it is compared against
  is non-Limber (caveat in section 5).
* (2.1) uses a locally power-law P at the hard scale; the running of
  n(k) across the hard support is an O(10%) effect on the pointwise
  coefficient (measured: 11%), absorbed when Rbar is used under the
  integral as in (3.3).
