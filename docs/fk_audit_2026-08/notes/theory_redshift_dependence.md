# Redshift dependence of the FK leakage: per shell and through the fold

> **Correction, 2026-09-02.** The absolute FK amplitudes in this note are
> `ell_max = 1000` values obtained with the sorting callable and the production
> linear-in-lambda fold. The manuscript's FK is evaluated at `ell_max = 15360`
> with the permutation-aware callable, where the same quantity is `+1.9480e-5`
> at `gamma = 0.5'` (2.31% of Order-0), a factor 4.5 larger. Ratios between
> variants at a fixed cutoff are unaffected. See
> [`../FK_BASELINE_NUMBERS.md`](../FK_BASELINE_NUMBERS.md) for the full key and
> `sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/` for the
> current product.

Companion to `theory_shape_amplitude_factorization.md`, which treats one
shell at fixed redshift and fixed source. This note makes the redshift
dependence explicit at both levels: the per-shell z dependence of the
collapsed vertex, and the source-redshift (z_s) dependence of the folded
xi_FK and of the ratio FK/O0. Every symbolic statement is machine-derived
(`code/verify_redshift_scalings.wl`, 15 checks, all PASS); every claim
about the pipeline is a measurement on the real dense tables
(`code/redshift_pershell.py`, `code/redshift_fold_trends.py`,
`code/probe_lambda_interp.py`), with the certification and control runs
described where they are used.

Notation as in the companion note; additionally D_+(z) is the linear
growth normalised to D_+(0) = 1 (canoes convention, the same function the
table builder applies), D_J = a chi is the background Jacobi map that
builds the fold windows, lambda the affine distance from the observer in
physical Mpc, and lambda_f = lambda(z_s) the source. The distinction
between D_+ and D_J matters: THE STATISTICS GROW WITH D_+, THE RESPONSES
PROPAGATE WITH D_J. All growth sits inside the per-shell table; the fold
weights are pure background geometry.

## 1. Per shell: the exact tree statement

### 1.1 Growth times geometry, with every explicit power listed

The HIGH branch of the table builder evaluates the tree bispectrum at
z = 0 and multiplies growth^4 (KAPPA3:3619-3626), so per shell, exactly
as implemented,

    zeta(gamma; lambda) = pref(z, chi) * D_+(z)^4 * Q_0(gamma; chi),   (1.1)

    pref(z, chi) = [A(a) (1+z)^4]^3 * (1+z)^-4 * (2pi)^-4 * chi^-4 * h^4,
    A(a) = -(3/2) Omega_m H_0^2 (1+z),                                 (1.2)

with Q_0 the z = 0 quadrature at geometry chi (all of section 0-3 of the
companion note), chi in Mpc/h, H_0 = 100/299792.458 h/Mpc. The origin of
each factor, source-read, not derived here: three Phi00 legs contribute
A(a)(1+z)^4 each (Poisson response times the Sachs photon-energy weight,
KAPPA3:1696-1714; the (1+z) factors are photon-energy bookkeeping in the
paper's E-primitive convention, not a Born substitution); (1+z)^-4 is the
radial-measure Jacobian (d lambda/d chi)^2 of the two equal-shell
collapse deltas; chi^-4 is the Limber geometry; h^4 converts the
(Mpc/h)-native spectral density to physical units.

Equation (1.1) factorizes REDSHIFT AWAY FROM GEOMETRY: z enters only
through pref and D_+^4, gamma and the cutoff only through Q_0. At tree
level this is not an approximation. It holds because linear growth is
scale-independent, B_tree(k; z) = D_+(z)^4 B_tree(k; 0), which the builder
implements literally.

Power-law limit (P = A k^n; wl checks V1-V3): combining (1.1) with
Q_0 ~ (2pi)^2 S H,

    zeta ~ (1+z)^11 D_+^4 chi^(-4-2n) gamma^(-(2+n)),                  (1.3)

the 11 being 3x5 from the legs minus 4 from the Jacobian. The vertex
grows steeply toward the source: between the innermost shell (z = 0.10)
and the outermost (z = 5.7) the measured zeta_TTT at 0.5 arcmin grows by
a factor ~4e4, the (1+z)^11 and chi^(-4-2n) factors pulling against
D_+^4 and winning.

### 1.2 Certification against the real table

`redshift_pershell.py` [A2] compares the table against
pref * D_+^4 * (brute-force z = 0 quadrature), i.e. (1.1) with NO
squeezed factorization anywhere. Result: 1.000 +/- 0.005 at
gamma <= 5 arcmin across ALL 16 shells, z = 0.10 to 5.7; within 1% at
10 arcmin (0.991 on the outermost shell); drifting to 0.94 at 19 arcmin
on the far shells (quadrature resolution plus the LOW admixture). A
single mis-read power of (1+z) would tilt this ratio by a factor
6.7/1.1 = 6.1 across the table; the observed tilt is under 0.5% at small
gamma. Every explicit redshift factor in (1.2), the growth power, and
the unit conversions are thereby certified end to end. (An independent
review quadrature, written from (1.1)-(1.2) without the project's code,
reproduced the table to 0.02-0.5% at the z = 1.0 and z = 3.1 shells.)

### 1.3 The factorized form degrades toward low-redshift shells

The same script's [A] column replaces the brute-force Q_0 by the
factorized (2pi)^2 S_cap H. The ratio table/prediction is 1.00 +/- 0.07
for shells at z >= 2.5 but rises smoothly to 1.85 at the z = 0.10 shell.
With [A2] flat, this drift is entirely the squeezed-factorization error
E, and its origin is geometric: at fixed ell in (60, 1000] the shell maps
the window to k in (60/chi, 1000/chi], which at chi = 293 Mpc/h means
k = 0.20-3.4 h/Mpc, most of which lies beyond the point k ~ 0.5 h/Mpc
where the per-log hard integrand u^2 Rbar P turns over (n = -2). Past
that point the hard support slides down to the ell_cut edge, the scale
hierarchy between the J_0 leg and the hard pair compresses, and the
neglected pieces of the error budget (the hard-hard pairing, the swapped
u-v region) turn on together. So THE
ACCURACY OF THE SHAPE-TIMES-AMPLITUDE PICTURE IS ITSELF A FUNCTION OF
SHELL REDSHIFT: excellent for z >~ 2, order-one at z <~ 0.25. Statements
built on the factorization inherit this domain.

### 1.4 Nonlinear input breaks growth-geometry factorization, low z first

For a nonlinear B the builder bypasses growth^4 and evaluates B at
(k, z_shell) directly, because no D_+^4 scaling exists. Measured
(script [B], BiHalofit vs tree, same grid, same cutoff, floor-matched):

    zeta_NL / zeta_tree at 0.5 arcmin:
      z = 0.10: 21.6     z = 0.45: 2.18    z = 1.0: 1.20
      z = 1.9:   1.04    z >= 3:   1.01

The nonlinear enhancement of the collapsed vertex is a LOW-REDSHIFT-SHELL
phenomenon: at fixed ell window the near shells sample k of order
1 h/Mpc where the one-halo term dominates, the far shells sample
k <~ 0.2 h/Mpc where tree is adequate. Equivalently, the effective growth
exponent d ln zeta_NL / d ln D_+ exceeds the tree value by ~18 at the
z = 0.1 shell and by < 0.1 for z >= 2.5. This single per-shell fact
drives the z_s dependence of the nonlinear correction in section 3.

## 2. The fold: geometry-weighted shells

### 2.1 The exact weights

The equal-time FK fold, transcribed from the validated independent
assembly (`fk_analytic.py`, eq. FK; reproduces the production sft-wick
fold to 0.1%):

    xi_FK(gamma; z_s) = INT_0^{lambda_f} dlam W(lam) u(lam)^4
                        V(lam) G_00(lam; gamma),                       (2.1)
    W(lam)  = INT_lam^{lambda_f} [D_J(lam)/D_J(t)]^2 dt,
    u(lam)  = D_J(lam)/D_J(lambda_f),
    V(lam)  = INT_lam^{lambda_f} W(t) u(t)^-4 dt,
    G_00    = -2 (zeta_TTT + zeta_Bmod),

and xi_O0 = INT dlam W(lam)^2 C_2(lam; gamma) for the two-point side.
(The implemented outer domain is the table support [lambda_min,
lambda_f], lambda_min the shell floor of section 4 item 3; the vertex is
zero below it by construction.)
All redshift dependence of the WEIGHTS is the background Jacobi map
D_J = a chi; all redshift dependence of the STATISTICS is inside G_00 via
section 1. (For kk the two spin-2 channels collapse into zeta_Bmod, and
the table shows zeta_Bmod = zeta_TTT to 0.6% at gamma <= 5 arcmin, the
helicity-sum identity, so the kk fold is TTT physics up to normalisation.)

`redshift_fold_trends.py` implements (2.1) directly on the dense tables.
Controls: (i) quadrature converged (identical to six digits for 96, 192
and 384 nodes); (ii) with the production's linear-in-lambda table
interpolation it reproduces the archived production fold (T1, next
subsection); (iii) the closure identity, that the weight-profile average
of the per-shell NL ratio reproduces the folded NL ratio exactly,
confirms the extracted weights are THE weights of the fold. (The closure
is exact by construction, a consistency identity of the weight
extraction, not an independent measurement.)

### 2.2 A measured systematic: the lambda interpolation of the shell grid

Where the fold puts its weight decides how accurately the 16-shell table
resolves it. Measured weight percentiles of the 0.5-arcmin FK fold
(z_s = 5): 25/50/75% of the fold lies below shell z = 0.60/0.88/1.25.
The fold is MID-PATH dominated, and the mid-path is exactly where the
GL-derived shell grid is sparsest (gaps of ~130-450 Mpc over which the
vertex grows by factors of ~2-3 per gap, a factor ~6 across the two
largest gaps combined).

Across such gaps the production callable interpolates LINEARLY in
lambda; a chord above a near-exponential curve overestimates the gap
interior. `probe_lambda_interp.py` settles which interpolation is right
by computing the TRUTH at five mid-gap lambdas with the certified chain
of section 1.2 (valid to 0.5% there):

    linear / truth   = 1.04-1.17  (growing with lambda, i.e. with the
                       local steepness of the vertex, and decreasing
                       with gamma)
    logPCHIP / truth = 0.989-1.003

Folded consequence (T1, tree table, production t_final = 2313.029,
compared at gammas present in BOTH the production sweep grid and the
table rows; the first version of this table paired gammas 11% apart from
the two interleaved grids and overstated the mid-angle inflation, caught
in the review round):

    gamma:               0.5'    4.2'    21.9'    56.3'   144.7'
    production / mine:   1.115   1.114   1.069    1.041    0.896

So the production xi_FK at z_s = 5 carries a lambda-interpolation
inflation of +11% at arcminutes, decaying monotonically to +7% at
~20 arcmin and +4% at ~1 degree, passing through zero near 2 degrees and
reaching -10% at 145 arcmin (where the vertex's lambda profile carries
sign structure); the grid-converged tree baseline at the production
conventions is +3.84e-6, not +4.28e-6. This is a property of the shell
GRID plus interpolation rule, not of the vertex; it was invisible to the
MC cross-check because both sides read the table through the same
interpolation. All trend numbers below use the logPCHIP fold on both
sides of every ratio (interpolation rule validated against truth at
mid-gap shells with z <= 1.7 and gamma <= 20 arcmin; above z ~ 1.4 the
shell gaps shrink from 127 Mpc to a few Mpc, so the uncertified high-z
gaps are also the smallest).

Two source-end conventions, recorded because they matter at this
precision: the production t_final = 2313.029 corresponds to z_s = 4.70 on
the closed-form background (lambda(5) = 2318.0), and the fold is
sensitive to the source convention at ~1.1% per Mpc of lambda_f (the
change reweights u^4 and the windows globally, not the last shells; the
quadrature is converged and the last 100 Mpc carry only 1.8% of the
fold). The trend rows below use the multiz lambda_f(z) convention
consistently on both sides of every ratio.

## 3. Measured source-redshift trends

Fold the SAME dense tables at lambda_f(z_s) for the 20-plane grid's
z_s = 0.5 to 5 (the tables' shells extend past lambda_f(5), so no new
builds are needed); O0 is the exact production 20-plane computation from
the multi-redshift figure data (`multiz_components.npz`). Selected rows;
full tables in `products/redshift/fold_trends.npz`.

### 3.1 Tree level: FK tracks O0 in z_s almost as tightly as in gamma

Amplitudes at 0.5 arcmin normalised to z_s = 5:

    z_s      O0        FK(tree)    FK/O0 relative
    0.50    0.0145     0.0121      0.839
    0.97    0.0763     0.0772      1.012
    1.68    0.2361     0.2471      1.047
    3.11    0.6211     0.6165      0.993
    5.00    1.0000     1.0000      1.000

O0 grows by a factor 69 and FK by 82 from z_s = 0.5 to 5; their RATIO
moves by -16% to +5%, non-monotonically (maximum near z_s ~ 1.7). In
absolute terms the tree FK/O0 at 0.5 arcmin is 0.40-0.50% at every
source redshift. The companion note's conclusion, that FK is a nearly
multiplicative dressing of the Gaussian signal, extends from the angular
to the source-redshift direction: TO A FIRST APPROXIMATION THE TREE FK
FRACTION IS A CONSTANT OF THE SURVEY, not a growing function of source
distance.

Why: in the fold, FK differs from O0 by one extra factor of the
transverse hard variance D_+(z)^2 H(chi)|_{z=0} (the hard factor of the
companion note, evaluated at z = 0; growth split off explicitly) under a
similar but not identical weight profile. The growth factor D_+^2 pulls
the ratio down as z_s grows; the geometry pulls it up. The power-law
anchor makes the GEOMETRIC half of this ledger exact bookkeeping: in the
low-z (growth set to 1), fixed-window, power-law limit the fold ratio
scales as chi_s^(-(n+1)) with NO gamma dependence (wl V4; closed forms
and Beta-function exponents verified against quadrature on both sides of
n = -1, V5a-c). The exponent decomposes as (-n) from the chi-rise of the
hard factor, (-2) from the three-point versus two-point Limber measure
(chi^-4 against chi^-2), and (+1) from the extra line-of-sight response
integral V, an accounting confirmed term by term in the review round;
growth is ABSENT from the anchor by construction, so its down-pull is
the additional real-universe effect the anchor does not capture. Within
the anchor the trend's sign is set by whether the effective slope n at
the sampled scales is above or below -1. Validity domains, stated once:
the soft-factor closed form needs -2 < n < -1/2; the two fold closed
forms hold for -2 < n < 0 (FK needs n < 0, O0 needs n < 1). The measured
local slopes d ln(FK/O0)/d ln chi_s at 0.5 arcmin run between -0.3 and
+0.5, i.e. an effective n scatter around -1, consistent with the mixed
scales the real fold samples. This anchor is a consistency illustration,
not a verification: the real fold mixes z and scale in both factors, and
growth is exactly what it omits.

### 3.2 The gamma direction, per source redshift

The corrected FK/O0 is flat-to-mildly-varying in gamma at EVERY z_s
(0.30-0.78% across 0.5-145 arcmin, all ratios taken at gammas shared by
both grids), while the deployed-table ratio rises with gamma to 26%
(z_s = 0.5) and 192% (z_s = 5) at 145 arcmin. The
apparent "FK becomes more important at wide angles, and more so at high
z_s" of the deployed figure is the frozen-corner artifact meeting a
decaying O0; nothing of it survives the dense grid at any source
redshift.

### 3.3 Nonlinear input: the one genuinely growing-with-1/z_s effect

Folded NL/tree ratio at 0.5 arcmin (floor-matched tables):

    z_s:    0.50   0.97   1.68   2.39   3.11   4.05   5.00
    NL/tree 7.18   3.83   2.55   2.12   1.91   1.76   1.68

combined with 3.1, the NONLINEAR FK fraction falls with source redshift:
FK_NL/O0 at 0.5 arcmin is ~2.7% at z_s = 0.5 against ~0.8% at z_s = 5.
(The 2.7% is the floor-397 BiHalofit fold divided by the full-path O0,
so it inherits the -6.6% lambda_min suppression at z_s = 0.5; the
floor-consistent value is ~2.9%. Both are at ell <= 1000, so
cutoff-conditioned in absolute terms; the DIRECTION is robust because it
reproduces the per-shell fact 1.4 through the closure identity, and
deepening the cutoff only adds more low-redshift small-scale power.) This is the physically meaningful
redshift statement of the corrected pipeline, and it points the OPPOSITE
way from the deployed figure's impression: the three-point leakage
matters most for LOW-redshift sources, because their kernels sit in the
late-time, nonlinear universe.

### 3.4 Where the fold weight sits

Weight percentiles of the 0.5-arcmin FK fold: for z_s = 5, half the fold
lies below shell z = 0.88; for z_s = 0.5, below z = 0.23. Even the
z_s = 5 observable is a late-universe integral: the steep vertex growth
toward the source (1.3) is tamed by the u^4 response suppression, so the
fold never becomes source-shell dominated. This is why per-shell
statements at z <~ 1.5 (where factorization degrades, 1.3, and
nonlinearity lives, 1.4) control the observable even for the highest
source redshifts.

## 4. The deployed multi-redshift trend: four measured distortions

The draft's multi-z statement ("the FK leakage strengthens steadily with
source redshift") was previously flagged untested. It can now be scored
against measurements. Four mechanisms distort the deployed trend, each
now measured rather than estimated:

1. THE COSINE-GRID ARTIFACT, dominant at gamma >~ 20 arcmin at every
   z_s (section 3.2): the deployed wide-angle FK/O0 rise with z_s
   (26% -> 192% at 145 arcmin) is entirely fake.
2. THE FIXED-ELL CUTOFF TILTS THE TREND. The per-shell undercount
   U = H(k <= 10 h/Mpc)/H(ell <= 1000) runs from 1.19 (z = 0.1) to 4.06
   (z = 5.7) because k_cut = 1000/chi falls with distance; fold-weighted,
   <U> = 1.39 at z_s = 0.5 against 2.14 at z_s = 5 (hard factor alone, a
   lower bound on the amplitude correction). A fixed-ell cutoff therefore
   suppresses high-z_s slices ~1.5x more than low-z_s slices: the
   deployed tree trend is TILTED DOWN toward high z_s, and any cutoff
   chosen in ell (rather than k) will redshift-bias a multi-z figure.
3. THE LAMBDA_MIN FLOOR, measured by folding the floor-50 against the
   floor-397 table: -6.6% at z_s = 0.5 falling to -0.1% at z_s = 5 at
   0.5 arcmin (up to -15% at 145 arcmin, z_s = 0.5). Suppresses low-z_s
   slices, the same direction as the draft's claimed trend.
4. THE SHELL-GRID INTERPOLATION (section 2.2): +11% at arcminutes
   decaying to +4% at ~1 degree and -10% at 145 arcmin at z_s = 5;
   smaller at low z_s where the fold sits partly on the denser
   low-lambda extension (logPCHIP/linear = 0.94 at z_s = 0.5, 0.5
   arcmin). Also tilts the trend, upward at high z_s.

At small gamma the deployed AMPLITUDE trend happens to land close to the
corrected tree trend (the four distortions partially cancel there:
compare the deployed 0.0411 with the corrected 0.0121 at z_s = 0.5,
normalised at z_s = 5, so the deployed multi-z figure OVERSTATES the
lowest-z_s FK by ~3.4x while roughly tracking above z_s ~ 1). The
qualitative reading that survives: absolute FK grows steeply with z_s
(as does O0); the FRACTION does not, at tree level, and falls with z_s
once the input is nonlinear.

## 5. What this note establishes, and what it does not

Established, with the verification layer in brackets:

* the exact per-shell redshift factorization (1.1)-(1.2) of the tree
  table, certified end to end on the real table [A2, 0.5%];
* the power-law exponents (1.3) and the low-z fold exponents, including
  FK/O0 ~ chi_s^(-(n+1)) and its gamma-independence, machine-verified on
  both sides of n = -1 [wl, 15 PASS];
* the shell-z dependence of the factorization accuracy [A] and of the
  nonlinear excess [B];
* the fold weights, their z_s-dependent profile, and the closure
  identity [T1/T4 controls];
* the interpolation systematic of the production fold, settled by a
  truth probe at mid-gap shells [probe; +11% at arcminutes, decaying
  through zero near 2 degrees];
* the measured z_s trends of xi_FK, NL/tree and FK/O0 [T2/T3].

Not claimed:

* No absolute amplitude: every FK number here is at ell <= 1000; the
  cutoff correction is itself z_s-dependent (section 4, item 2), which
  is one more reason the production rebuild must fix the cutoff in k,
  not in ell.
* The O0 denominator is the exact production 20-plane computation; its
  lambda_f convention differs from the production FK baseline's at the
  percent level (section 2.2), and FK/O0 ratios inherit that.
* Spin-2 channels beyond the kk-relevant zeta_Bmod (= zeta_TTT to 0.6%
  at small gamma) are untouched, as in the companion note.
* The truth probe certifies gamma <= 20 arcmin and z(lambda) <= 1.7
  (HIGH-dominated cells); the wide-angle interpolation numbers at 56 and
  145 arcmin (+4% and -10%) are measured at matched gammas but not
  truth-decomposed, since the LOW admixture is outside the brute-force
  chain; and the logPCHIP rule is extrapolated to the z > 1.7 shell
  gaps, which are however the smallest gaps of the grid (127 Mpc down
  to a few Mpc).
* The multi-z GIF's z_s < 1 deployed values are PCHIP extrapolations in
  the talk data; the 3.4x overstatement at z_s = 0.5 is relative to
  that product, not to a from-scratch deployed fold.
