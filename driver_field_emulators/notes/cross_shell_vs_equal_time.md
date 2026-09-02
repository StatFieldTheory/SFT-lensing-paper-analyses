# Cross-shell vs equal-time: does a single-z bispectrum emulator suffice?

**Question.** Public emulators output equal-time statistics, e.g.
B_delta(k1, k2, k3; z) at one redshift. The SFT lensing formalism is written in
terms of cross-shell cumulants, e.g. the three-point cumulant
K^(3)_abc(n1, n2, n3; lambda1, lambda2, lambda3) with three independent affine
labels. Is a single-z emulator therefore structurally insufficient?

**Answer: no, for three reasons, ordered by how the pipeline actually works.**

## 1. The production 3-point vertex is ALREADY equal-shell collapsed

The paper (sections/appendix.tex, `append: driving-field spectra`) replaces the
cross-shell cumulant by its diagonal before anything is computed:

    K^(3)_abc(gamma; lambda1, lambda2, lambda3)
      ~= zeta_abc(gamma, lambda1) * delta(lambda1 - lambda2) * delta(lambda1 - lambda3)

(Eq. `appendix equal shell approx`). This is the three-point counterpart of the
Limber approximation: the lensing response kernels are radially broad compared
with the coherence length of the matter correlators (~tens of Mpc), so
unequal-time corrections sit well below the percent level at our scales
(the paper's wording; the standard extended-Limber scaling is
O((Delta chi / chi_s)^2)) [LemosEtAl2017, KitchingNonLimber2017,
FangEtAl2020].

Concretely, the deployed callable is literally named
`callables/kappa3_vertex/equal_time_limber/`: it consumes a per-shell,
single-redshift matter bispectrum, currently tree-level SPT with
B_delta(z) = D(z)^4 B_delta(0) (Eq. `appendix bdelta growth`). So an equal-time
emulator B_emu(k1, k2, k3; z), evaluated shell by shell at z(lambda), slots into
Eq. `appendix limber reduced bispectrum` by replacing

    D(z)^4 * B_tree(l1/chi, l2/chi, l3/chi)   -->   B_emu(l1/chi, l2/chi, l3/chi; z)

with ZERO loss relative to the current approximation level. Everything
downstream (Poisson factor A(a)/k^2 per leg, screen-Hessian multipliers
L^2/chi^2 and sqrt(L^2(L^2-2))/chi^2, Gaunt sum, (1+z)^-4 lambda-measure
Jacobian on the two collapse deltas, h^4 normalization) is untouched.

## 2. Where cross-shell genuinely survives, it is built from equal-time inputs

The two-point propagator C = K^(2)(gamma; lambda', lambda'') IS kept fully
two-time in production (`callables/C_propagator/corr_op/`). But its input is a
single z=0 power spectrum (`PCAMB_pyccl_stf_fid_z0.txt`); the two times enter
through per-leg linear transfer operators, i.e.

    P(k; z1, z2) = D(z1) D(z2) P(k; 0)      (exact at linear order).

So "cross-shell" in this pipeline has never meant "cross-shell measured
statistic"; it means "equal-time input + per-leg time transfer". The same
pattern extends to the 3-point vertex if we ever relax the equal-shell
collapse (next section).

Caveat for nonlinear 2-point inputs: the exact factorization breaks beyond
linear theory. The standard, tested fix is the geometric-mean ansatz

    P_NL(k; z1, z2) ~= sqrt( P_NL(k; z1) * P_NL(k; z2) ),

whose error on lensing observables is far sub-percent (Kitching & Heavens 2017,
Phys. Rev. D 95, 063522, arXiv:1612.00770; de la Bella, Tessore & Bridle 2021,
arXiv:2011.06185). An
equal-time P(k) emulator (EuclidEmulator2, bacco, ...) therefore feeds the
cross-shell C without structural obstruction.

## 3. Lifting an equal-time bispectrum emulator to cross-shell, if ever needed

The FK 2-point topology pins the three vertex legs to the TWO external
directions separated by gamma (two legs coincident on one ray, the third on
the other: the (1, cos gamma, cos gamma) family; a 2026-08-25 code probe
confirmed production queries exactly this, refuting the earlier 2026-05-14
"one spatial point" colocation reading, which was a rounding artifact of a
single 0.5' probe). The three affine times are independent, so the
beyond-equal-shell object is the unequal-time bispectrum near the time
diagonal, B_delta(k1, k2, k3; z1, z2, z3) with z1 ~= z2 ~= z3. Hierarchy of
lifts from equal-time products:

(a) Tree level: closed form, no emulator needed (exact within the
    separable-growth / EdS-kernel approximation the paper already adopts in
    Eq. `appendix bdelta growth`; in LCDM the second-order growth deviates
    from D^2 and F2 acquires weak time dependence, both percent-level). With
    delta(k, z) = D(z) delta1(k) + D(z)^2 delta2(k) + ...,

    B_tree(k1,k2,k3; z1,z2,z3)
      = 2 F2(k1,k2) D(z1) D(z2) D(z3)^2 P(k1) P(k2) + 2 cyclic,

    which reduces to D^4 B(0) on the diagonal. The cross-shell tree-level
    vertex is therefore fully determined by D(z) alone.

(b) Nonlinear, near-diagonal: per-permutation weighted lift, NOT the naive
    cube-root geometric mean. Because the tree bispectrum is a SUM of three
    cyclic terms with unequal growth exponents (1,1,2), the cube-root ansatz
    [B(z1)B(z2)B(z3)]^(1/3) has a residual FIRST order in the time-offset
    differences whenever the permutations are unequal (i.e. precisely for
    squeezed shapes); it improves to second order only after the symmetric
    radial fold. The clean ansatz applies the per-leg growth of item (a) to
    each cyclic term separately (exact at tree level by construction) and
    rescales each term by its equal-time nonlinear-to-tree ratio;
    alternatively Taylor-expand in the time offsets with response-theory
    coefficients.

(c) True cross-shell nonlinear emulator: does not exist publicly (matter
    statistics are measured per snapshot in simulations; unequal-time
    correlators would require cross-snapshot correlation of the same modes,
    which no public emulator provides). If ever required, route (b) calibrated
    against a small set of dedicated cross-snapshot measurements is the
    realistic path.

## Power conservation: property of the projection, not of the input

Worry: the pipeline's equal-shell collapse "conserves power"; does an emulator
input preserve that?

The delta-function representation conserves power iff the delta coefficient
equals the radial integral of the true cross-shell cumulant,

    zeta(gamma, chi) = INT dDelta2 dDelta3
                       zeta_3D(r_perp; Delta2, Delta3; z(chi)),

NOT the coincident-shell (on-diagonal) value, which is a variance-like,
UV-divergent quantity and would not conserve anything. In Fourier space the two
radial integrals project the bispectrum onto k_parallel = 0; by isotropy that
slice is exactly the equal-time 3D bispectrum on closed transverse triangles,
B(l1/chi, l2/chi, l3/chi; z).

The pipeline's zeta is built precisely this way: Eq. `appendix limber reduced
bispectrum` has B^XYZ = (D^4/chi^4) K B_delta(l_i/chi), where the 1/chi^4 (the
3-point analogue of the 2-point Limber 1/chi^2), the (1+z)^-4 lambda-measure
Jacobian on the two collapse deltas, and the h^4 normalization ARE the
integrated-weight bookkeeping. Power conservation therefore lives entirely in
the projection formula and its measure factors; the input B is only asked for
its value on the k_par = 0 slice. Swapping B_tree -> B_emu leaves the
conservation structure untouched, for ANY B.

Historical confirmation: the (1+z)^-4 radial-measure bug and the h^6 -> h^4
bug (both fixed 2026-06-09/10) were exactly violations of this bookkeeping,
caught by h-invariance / boundary tests. Those same tests remain valid
conservation checks after an emulator swap.

Emulator-specific caveats (fidelity, not conservation):
1. Squeezed-triangle coverage: the FK collapsed geometry weights the
   near-degenerate hard pair l1 ~= l2, with the soft leg l3 carrying
   P_{l3}(cos gamma) (Gaunt sum over all parity-even triples; the paper's
   index convention). If the emulator clips or
   extrapolates the squeezed corner of its training domain, the low-ell /
   wide-angle part of zeta is distorted.
2. k-range clipping: k_max ~ ell_max/chi_min may exceed the training domain;
   blend to tree-level/BiHalofit outside the domain rather than hard-clipping,
   otherwise the k_par=0 slice power is artificially altered.
3. Time variation across the coherence length: on the lightcone, Delta chi also
   means Delta t; the fixed-time radial integral has error
   O(xi_r * dlnD/dchi) ~ 1%. Identical for tree and emulator inputs; part of
   the standard Limber budget.

Validation recipe after a swap: (a) closure test, i.e. an independent Limber
projection of the SAME B_emu compared against the pipeline zeta (the fastnc
input-validation pattern); (b) boundary tests at the emulator domain edges
(both endpoints + extreme corners, per the boundary-validation methodology).

## Where the approximation actually bites

The equal-shell collapse is a Limber-class approximation, so its error
concentrates at LOW multipoles / wide angles. That is also where the FK
enhancement lives (crossover ell ~ 10, gamma >~ 2 deg), so the honest statement
is: at ell <~ 10 the collapse itself, not the equal-time restriction of the
emulator, is the leading systematic. (The figures ~0.6% at ell ~ 10 and
10-15% at ell ~ 2 come from the internal low-ell assessment of 2026-06-15,
not from the paper; the paper's own low-ell caveat, conclusion.tex, is
qualitative and concerns a different systematic, the near-horizon gauge
terms. Both sit below the tree-level-bispectrum model error.) An emulator
swap improves the k-space (nonlinearity) accuracy of the
vertex at all angular scales (wide angles included, through the coincident
hard-leg pair of the collapsed configuration; see coverage_ell_z_mapping.md)
without changing this radial-collapse error budget.

## Bottom line

- Equal-time emulator == exactly what `equal_time_limber` consumes today.
- Cross-shell 2-point: equal-time P(k) + per-leg growth (exact linear) or
  geometric mean (nonlinear, sub-percent).
- Cross-shell 3-point beyond equal-shell: tree-level exact via D(z); nonlinear
  via cube-root geometric mean; no public product does better, and none is
  needed at the paper's current error budget.
