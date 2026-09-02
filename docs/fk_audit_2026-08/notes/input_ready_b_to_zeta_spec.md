# Input-ready spec: B_model(k1,k2,k3;z) -> configuration-space zeta table -> sft-wick coupling_fn

> **Correction, 2026-09-02.** The absolute FK amplitudes in this note are
> `ell_max = 1000` values obtained with the sorting callable and the production
> linear-in-lambda fold. The manuscript's FK is evaluated at `ell_max = 15360`
> with the permutation-aware callable, where the same quantity is `+1.9480e-5`
> at `gamma = 0.5'` (2.31% of Order-0), a factor 4.5 larger. Ratios between
> variants at a fixed cutoff are unaffected. See
> [`../FK_BASELINE_NUMBERS.md`](../FK_BASELINE_NUMBERS.md) for the full key and
> `sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/` for the
> current product.

What sft-wick consumes is NOT the Fourier-space bispectrum: it is the
configuration-space three-point cumulant of the driving fields, served as a
six-channel table zeta_ch(cosine triple; lambda shell) through
`coupling_fn(n_list, t_list) -> (3,3,3)`. This note specifies, at the
formula/file:line level, the complete map from an equal-time matter
bispectrum B_model to that table, so the rebuild can be coded directly.
Everything below was extracted from source on 2026-08-25 (canoes at
`/Users/zzhang/projects/angular_statistics/canoes`, src layout; sft-wick at
`/Users/zzhang/projects/SFT/sft-wick`) plus an instrumented production-config
probe. KAPPA3 = `src/canoes/sachs/kappa3.py` (5191 lines); CALLABLE =
`equal_time_limber_kappa3_callable.py`; BUILD =
`build_equal_time_limber_table.py` (both in
`SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber/`).

## 0. Division of labor (empirically probed, production FK config)

sft-wick does, per (gamma, component pair, FK term):
- a 4-D tensor-product Gauss-Legendre integral, n_gauss=24 -> 24^4 = 331,776
  nodes per call (2 integrated external times + F-vertex time + K-vertex
  collapsed time); 2 FK terms x 6 pairs x 40 gammas = 480 batched calls.
- equal_time aliasing: the K vertex's three leg times are one variable
  (verified bit-identical across legs); ONE `int dlambda` for the vertex.
- multiplies in itself: the response R(t,lambda') = Theta (D(lambda')/D(t))^2
  per retained R line (the FK diagram terms have 4 R lines and ZERO C
  lines: c_product = 1, verified by re-running the symbolic expansion with
  production flags; the corr_op C feeds other diagrams, not the FK
  integrand) and the MSR prefactor +i/6. Phase accounting (verified): the
  K vertex's three response lines carry (-i)^3, the F line its own -i;
  combined with +i/6 and the response phase (-i)^4 the net factor is REAL
  +1/6, so a real positive table value enters with no residual phase; the
  physical sign is carried by the (all-negative) F_bare tensor and the
  table values themselves. sft-wick enforces a reality check on the
  contracted coupling.

Therefore the TABLE VALUES must be: the BARE equal-shell collapsed 3-leg
cumulant as a LAMBDA density (physical Mpc), with the (dlambda/dchi)^2 =
(1+z)^-4 collapse Jacobian and the h^4 normalization already folded in, and
NOTHING else (no R, no i/6, no extra measure factors). The callable's only
runtime math is interpolation + channel-to-tensor reconstruction.

## 1. The target object: six channels, general-triple angular kernel

Channels stored (each real float64, shape (n_triples, n_shells)):

| Channel | Cumulant | spins (s1,s2,s3) | collapsed-family reduced fn of c |
|---|---|---|---|
| zeta_TTT | <Phi00^3> | (0,0,0) | P_{l3}(c) |
| zeta_TTP | <Phi00^2 Psi0> | (0,0,2) | d^{l3}_{0,-2}(arccos c) |
| zeta_TPP | <Phi00 Psi0 Psi0> | (0,2,2) | d^{l3}_{2,-2} |
| zeta_PPP | <Psi0^3> | (2,2,2) | d^{l3}_{4,-2} |
| zeta_Bmod | <Phi00 Psi0 conj(Psi0)> | (0,2,-2) | d^{l3}_{2,+2} |
| zeta_Dmod | <Psi0 Psi0 conj(Psi0)> | (2,2,-2) | d^{l3}_{4,+2} |

The angular kernel (LOW: SpinAwareZetaEvaluator kernel with canonical frame
n1-at-pole, `_spin_aware_three_pt.py:566-712`; HIGH: flat-triangle +
spin-phase Bessel alpha kernel, KAPPA3:1732-1785) consumes ONLY the cosine
triple, ell list, and spins - it is provably B-independent, so the nonlinear
swap cannot change the angular layer. Conventions that must be kept as-is:
the W3j_0 parity anchor on all channels (matches
derive_kappa3_cross_cumulants.wl), the d^l_{m,-s} index order, phi3 in
[0, pi] with Re[] projection. Bmod/Dmod differ from TPP/PPP only in the
third spin sign (helicity SUM vs difference); downstream reconstruction:
K000=TTT, K001+perms=TTP, K011=(TPP+Bmod)/2, K022=(Bmod-TPP)/2,
K111=(PPP+3 Dmod)/4, K122=(Dmod-PPP)/4, odd-#-of-index-2 entries exactly 0.

## 2. HIGH branch (Limber, max leg in (60, 1000]): where B_model plugs in

Per cosine triple m and shell (chi in Mpc/h, z), per channel X
(`compute_kappa3_sigma3_high`, KAPPA3:3445-3723; mod variant :3924-4065):

    zeta_X_h(m, shell) = jac_lambda * (2pi)^-4 * R_X(z) * chi^-4
        * SUM_{i,j,p} w_i w_j w_p u_i v_j W_high * B_delta(k1,k2,k3;z)
        * alpha_X(u_i, v_j, phi_p; m)
    then * h^4  (physical units)

- Quadrature: n_ell=96 Gauss-Legendre nodes LINEARLY mapped on
  [1e-3, 1000] (NOT log) for u and v; n_phi=64 GL on [0, 2pi]; grid
  (96,96,64) = 589,824 nodes; third leg w = sqrt(u^2+v^2+2uv cos phi);
  hard window W_high = 1 iff 60 < max(u,v,w) <= 1000 (KAPPA3:3548-3565).
- Per-leg kernel: the screen-Hessian ell^2/chi^2 CANCELS the Poisson 1/k^2
  exactly under Limber (k = ell/chi), leaving per leg
  Phi00: +A(a)(1+z)^4, Psi0: -A(a)(1+z)^4, with
  A(a) = -1.5 Omega_m H0^2 (1+z), H0 = 100/299792.458 h/Mpc
  (KAPPA3:1696-1714, 3537-3539). So
  R_X = (A(1+z)^4)^3 * (-1)^{#Psi0 legs}: one scalar per (channel, shell);
  NO explicit ell factors survive in the kernel.
- Angular factor alpha_X: map the triple to a flat 2D triangle
  (KAPPA3:1732-1744), then the orientation-averaged spin-phase Bessel
  kernel: S = s1+s2+s3; S=0 -> 2pi J_0(|Q|); S!=0 ->
  Re[2pi i^S e^{iS arg Q} J_S(|Q|) e^{i(s1 arg(-(u+v e^{i phi})) + s3 phi)}]
  with Q = u r2 + v e^{-i phi} r3 (KAPPA3:1747-1785). General open
  triangles, no collapsed special-casing.
- B_delta evaluation (KAPPA3:3607-3618): ONCE PER SHELL, shared by ALL 1671
  triples and all six channels: k1 = w/chi, k2 = u/chi, k3 = v/chi
  (h/Mpc; safe-floored at ell_min=1e-3), tree value
  b = 2(F2 P P + cyc) at z=0, then growth**4 (D = M[0]/(1+z), D(0)=1) at
  KAPPA3:3619-3622.

### Insertion point (minimal diff, BOTH sites together)

Site 1: KAPPA3:3612-3626; site 2 (mod): KAPPA3:4023-4036. Replace

    p1,p2,p3 = pk(k1..k3); f.. = F2(..); b_delta_today = 2(f12 p1 p2 + ...)
    radial_base = measure_weights * growth**4 * b_delta_today / chi^4 * pref

by

    b = b_delta_fn(k1, k2, k3, z_shell)      # (96,96,64), (Mpc/h)^6, k h/Mpc
    radial_base = measure_weights * b / chi^4 * pref   # growth**4 BYPASSED

b_delta_fn contract: accepts three float64 arrays (96,96,64) + scalar z,
returns same shape; must be symmetric in its k arguments (no sorting done);
may return garbage where W_high=0 (multiplied by zero; optional mask).
Domain actually exercised where W_high=1: k <= 1000/chi_h, worst case
~3.8 h/Mpc at the z=0.1 shell (chi_h ~ 266 Mpc/h); k_min ~ 4e-7 h/Mpc
(safe floor). Evaluation budget: 589,824 x 16 shells ~ 9.4e6 B calls per
full table (NOT x1671; triple dependence is entirely in the precomputed
alpha kernels). b_delta_fn=None must reproduce today's path bit-for-bit.

## 3. LOW branch (exact 3j, all legs <= 60): separability constraint

The LOW reduced bispectrum is computed by per-leg 1D FFTlog k-integrals
that EXPLOIT tree-level separability: the 18-monomial F2 expansion
(`_phi_accumulator.py:53-106`, layout `_phi_fftlog_data.py:68-99`), P-leg
FFTlog power-law decomposition + hypergeometric double-Bessel kernels
(`_phi_fftlog_legs.py:162-309, 590-597`), closure-leg Gamma formula
(:312-445), per-term product contraction (:650-662), and the z-separable
D^3->D^4 spectator promotion (KAPPA3:4854-4865). A general non-separable
B_model cannot enter here pointwise.

Options: (1) KEEP TREE in LOW (recommended; all legs <= 60 means
k <= 0.026 h/Mpc for z >= 0.3 shells, only the z=0.1 shell reaches
0.21 h/Mpc with an 18% one-leg boost); (2) brute-force triple k-integral:
~1e4-1e5 x current LOW cost, spot-check oracle only; (3) separable
re-expansion of B_NL in SC01/GM12 effective-kernel form
(B_NL = 2 F2_eff P_NL P_NL + cyc with per-leg dressings a,b,c), reusing the
monomial machinery with three dressed FFTlog decompositions per shell,
~3x current LOW cost - the fallback if the seam test fails at z <= 0.25.

## 4. The ell_cut = 60 seam

Hard complementary windows (LOW: all legs <= 60; HIGH: max leg in
(60,1000]), plain channel-wise sum in `kappa3_combine_low_high`
(KAPPA3:4068-4132) after guards on radial_measure/units/grids. CAVEAT: the
Bmod/Dmod channels bypass the guarded combiner (bare `+` in BUILD); the
rebuild must extend the guards to them. Seam boundary test (bypass the
combiner): band-difference LOW(60)-LOW(50) vs HIGH(ell_cut=50,
ell_high_max=60), all six channels, per triple and shell, rel_err gates
1e-3/1e-1, including the (1,c,c) family at both gamma extremes, open
triangles, and both extreme shells; run tree-B first as control, then
nonlinear (double duty: dispatch AND separable-expansion fidelity). Nearest
existing machinery: `ell_band_decomp.py:110-131` (self-consistency median
7e-4).

## 5. Units and measure ledger (with deployed-table provenance caveats)

- P(k) input: `PCAMBz0.txt`, k in h/Mpc, P in (Mpc/h)^3, z=0
  (`_local_cosmo_pk.py:30-54`). NOT the 2-point side's
  PCAMB_pyccl_stf_fid_z0.txt (which integrates to sigma8 = 0.810 exactly).
- Internal chi in Mpc/h; native zeta density (h/Mpc)^4; h^4 -> physical
  Mpc^-4 (KAPPA3:3650-3654 HIGH, 4881-4891 LOW, 3916-3919 / 4060-4063
  mods); lambda-density Jacobian (dlambda/dchi)^2 = (1+z)^-4 applied at
  build via `_kappa3_radial_density_conversion` (KAPPA3:1788-1844).
- PROVENANCE CAVEAT (verified numerically): the deployed NPZ's VALUES are
  post-Jacobian (per-shell deployed/PREJAC ratio = (1+z)^-4 to machine
  precision) but its meta strings and README still carry pre-fix text,
  including a stale `fk_kk_baseline = '+3.086e-5'`; the BUILD docstring's
  '+1.219e-4' is ALSO stale. The true tree FK kk(0.5') baseline of the
  deployed convention is +4.29e-6. A rebuild from current source
  reproduces the deployed values directly and writes correct meta; do NOT
  re-apply (1+z)^-4 on top, and re-derive the baseline, never copy stamps.

## 6. The actual query set, and the grid the rebuild must provide

- Cosine demand: EXACTLY the (1, c, c) family. All 480 production calls
  contain one unique sorted triple each; the 40 production c values run
  from 1 - 1.06e-8 (gamma=0.5') to 0.116 (gamma=5000') (full list in
  config_L2.yaml:61-100).
- CRITICAL GRID FINDING: covgrid16's cosine grid is 15 uniform values,
  spacing 1/7 ~ 0.143. Lookup is cKDTree k=4 inverse-distance weighting.
  Against the production queries (independently recomputed 2026-08-25):
  for gamma = 0.5'..232' (exactly 27 of 40 points) the IDW weight on the
  single (1,1,1) corner cell is >= 0.96 (>= 0.999 below 28'); the
  exact-hit branch triggers for 0 of 40 queries. Since the FK integrand
  has no C lines and R depends on times only, the computed kk FK curve
  below ~232' is grid-flat BY CONSTRUCTION up to the <= 4% IDW weight
  drift: the table cannot resolve any true variation of zeta(1,c,c) over
  1-c in [1e-8, 4.6e-3] there. Whether the TRUE curve is flat is exactly
  what the densified tree control rebuild decides. For gamma >= 294' the
  lookup carries 3.4-34.3% off-family (non-(1,x,x)) contamination,
  exceeding 10% from gamma ~600' and peaking near 3100'.
  REQUIREMENT: the rebuilt table must include (1,c,c) rows
  at (at least) the 40 production c values - ideally log-spaced in (1-c)
  over [1e-8, 0.9] - so every query's neighborhood is on-family. A
  TREE-LEVEL rebuild on the densified grid is REQUIRED as control, to
  separate grid effect from nonlinear effect in the FK curves. (Note:
  stored triples are rounded to 8 decimals, so the exact-hit branch
  (dist < 1e-10) will not trigger; on-family nearest neighbors with
  weight ~1 are sufficient.)
- Lambda demand: queries densely fill ~[5.6, 2307.5] Mpc (verified: this
  is exactly the GL-24 node span on [0, t_final = 2313.0289]); the table's
  shells start at 396.63 Mpc and everything below is ZEROED by the
  out-of-range guard - 6 of the 24 GL nodes (~25%) fall below the floor,
  so the near-observer vertex contribution is silently dropped, a
  gamma-independent systematic on the same curves. Either extend the shell
  grid downward or document the zeroing as intended (quantify first: zeta
  is smallest at low z but the decision should be explicit).

## 7. Replacement-table contract checklist (enforced fail-loud)

1. Arrays: `cosine_triples` (n,3); `lambda_shells_Mpc` ascending, physical
   Mpc, lambda_project; six channels (n, N_lam) real float64 - a missing
   channel is silently ZERO-FILLED (CALLABLE:133-135), so never omit one.
2. `cosmo_meta` JSON (uint8 bytes; json.loads(d['cosmo_meta'].tobytes()))
   with `already_R_contracted: false`, `equal_time: true`,
   `radial_density_measure: 'lambda'`, `has_modulus_channels: true`
   (each guarded, CALLABLE:90-122).
3. Values: bare, R-free, lambda-density, (1+z)^-4 and h^4 folded in, real.
4. Module contract: `coupling_fn(n_list, t_list) -> (3,3,3)` and
   `coupling_fn_batch(n_arr(3,ns,3), t_arr(3,ns)) -> (ns,3,3,3)` with
   `.vectorized = True`; wired BY PATH in the sft-wick YAML
   (coupling_module/coupling_attr/coupling_vectorized/equal_time).
5. Grid coverage: cosine rows making every (1,c,c) query on-family (sec 6);
   lambda span vs [5.6, 2307.5] Mpc decided explicitly.

## 8. Newly discovered issues the rebuild must address (beyond the swap)

1. Corner-dominated interpolation at 27/40 gammas: kk FK curves below
   ~232' are vertex-side grid-flat by construction (sec 6) - affects
   CURRENT paper FK curves; densified tree control rebuild required.
2. Lambda floor zeroing below 396.63 Mpc drops 6/24 GL nodes (sec 6).
3. Bmod/Dmod combine bypasses the convention guards (sec 4).
4. Stale meta/baseline strings on the deployed table and BUILD docstring
   (sec 5); tree-regression gate is +4.29e-6, not +1.219e-4.
5. ell_high_max=1000 convergence unproven for tree AND boosted under the
   P^1 toy scaling (see coverage_ell_z_mapping.md); sweep {1000..4000}.
