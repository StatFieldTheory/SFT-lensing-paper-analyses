# MC Sachs 2-point: numerical quantities, definitions, and analysis-3 comparison

Brief reference for `analyses/mc_sachs_2pt/`. Every quantity is given with its
exact mathematical definition and the code symbol that implements it. All lengths
are physical Mpc; the affine parameter `lambda` is 0 at the observer and increases
toward the source `lambda_f`. Cross-checked against the code (see the final
"Consistency" section).

Component convention everywhere: `s = (s_0, s_1, s_2) = (Phi_00, Re Psi_0, Im Psi_0)
= (kappa-rate, gamma_+-rate, gamma_x-rate)`. A two-ray bundle uses directions
`n_1, n_2` with `n_1 . n_2 = cos gamma`; the 6-vector is `f = (f(n_1), f(n_2))`,
indices `0..2` = ray 1, `3..5` = ray 2.


## 0. Dynamics and the observable

Stochastic Sachs fluctuation ODE (per ray; `F` is the local quadratic vertex):

    d s_a / d lambda = -2 theta_sa(lambda) s_a + F_abc s_b s_c + f_a(n, lambda).   (0.1)

Convergence/shear observable (the `integrate_over='all'` weak-lensing observable):

    kappa_a(n) = int_0^{lambda_f} s_a(n, lambda) dlambda.                          (0.2)

F-vertex (paper Table 1 = `scripts/inputs/F_tensor.npy`; nonzero entries):

    F_000 = -1, F_011 = -1, F_022 = -1, F_101 = -2, F_202 = -2,                    (0.3)

so `F_abc s_b s_c` = `( -(s_0^2+s_1^2+s_2^2), -2 s_0 s_1, -2 s_0 s_2 )`.
Code: `sachs_mc_core._F`, `sachs_mc_core._f_vertex`.

The 2-point observable is the cross-ray correlation
`xi_ab(gamma) = <kappa_a(n_1) kappa_b(n_2)>`; all estimators return the (3,3)
cross-ray block.


## 1. Background quantities  (`background.py`, over `scripts/D_callable.py`)

| Quantity | Definition | Code |
|---|---|---|
| Jacobi amplitude | `D(lambda) = a(lambda) chi(lambda)` (flat FLRW closed form) | `Background.D` |
| Focusing rate | `theta_sa(lambda) = d ln D / dlambda = D'/D` | `Background.theta_sa` |
| Linear response (Green's fn) | `R(t, lambda') = Theta(t - lambda') [D(lambda')/D(t)]^2` | `Background.response_R` |
| LOS window | `W(lambda) = int_lambda^{lambda_f} [D(lambda)/D(t)]^2 dt = int_lambda^{lambda_f} R(t,lambda) dt` | `driver_stats.order0_window` |

`R` is the exact Green's function of the linear (F=0) part of (0.1): `dR/dt =
-2 theta_sa(t) R`, `R(lambda',lambda') = 1`. The linear solution of (0.1) is
`s_a^(1)(n,lambda) = int_0^lambda R(lambda,lambda') f_a(n,lambda') dlambda'`, hence
from (0.2) `kappa_a^(1)(n) = int_0^{lambda_f} W(lambda') f_a(n,lambda') dlambda'`
-- `W` is the window with which the observable sees a driving-field impulse.
Cosmology hardcoded in `D_callable`: flat LCDM, `Omega_m = 0.31609, h = 0.6711`,
`z_s = 5` -> `lambda_f = 2313.029` (analysis-3 baseline) or `2317.96` (D_callable).


## 2. Bare driving-field cumulants  (`driver_stats.py`)

The driving field is treated as white in `lambda` with angular cumulant densities.

### 2a. 2-cumulant `Sigma2`  (`Sigma2Builder`, `sigma2_matrix`)

    <f_a(n_1, lambda) f_b(n_2, lambda')>_c = Sigma2_ab(n_1.n_2; lambda) delta(lambda - lambda').  (2.1)

`Sigma2` is extracted by DE-WINDOWING the corr_op propagator
`C_ab(cos; lambda_1, lambda_2) = <Phi_a(n_1,lambda_1) Phi_b(n_2,lambda_2)>` (the
R-windowed 2-cumulant). Using `R = [D/D]^2` (scalar, iso_R):

    Sigma2_ab(cos; lambda) = D(lambda)^{-4} d/dlambda [ D(lambda)^4 C_ab(cos; lambda, lambda) ].  (2.2)

Identity: `C_ab(cos; t,t) = int_0^t [D(la)/D(t)]^4 Sigma2_ab(cos; la) dla`, so (2.2)
is its exact inverse. Code: `Sigma2Builder.matrix` splines `D^4 C_ab(t,t)` over
`lambda` (corr_op queried at `C_fn(n_1,t,n_2,t)`) and differentiates. Optional
anchor `ANCHOR_C0 = 1/0.881 = 1.135` (`apply_c0`) rescales `Sigma2` so the
white-noise Order-0 matches analysis-3 in absolute normalization (Sec. 5).

Two-ray (6,6) assembly `Sigma2_6x6` (`sigma2_6x6`, `_assemble_6x6`): diagonal 3x3
blocks = within-ray `Sigma2(cos=1)`, off-diagonal blocks = cross-ray
`Sigma2(cos=gamma)`.

### 2b. Order-0 analytic prediction `O0`  (`order0_mc`)

The linear (F=0) limit of `<kappa kappa>`; with (2.1) the double LOS integral
collapses to a single one:

    O0_ab(gamma) = int_{lambda_min}^{lambda_f} Sigma2_ab(cos gamma; lambda) W(lambda)^2 dlambda.  (2.3)

Code: `order0_mc` (Gauss-Legendre, `lambda_min = lam_lo = 406`, the corr_op support
floor). This is the MC's EXACT linear target.

### 2c. 3-cumulant `zeta`  (`zeta_tensor`, `zeta6`)

Equal-time third cumulant (delta-correlated in lambda):

    <f_a(n_1,lambda) f_b(n_2,lambda') f_c(n_3,lambda'')>_c
        = zeta_abc(n_1,n_2,n_3; lambda) delta(lambda-lambda') delta(lambda-lambda'').  (2.4)

`zeta_tensor(n_1,n_2,n_3,lambda)` wraps the project `equal_time_limber` kappa3
vertex: a symmetric (3,3,3) tensor indexed by the cosine-triple `(n_i . n_j)` and
`lambda` (`ALREADY_R_CONTRACTED=False`, `EQUAL_TIME=True`, lambda-density). The
two-ray (6,6,6) `zeta6` places each leg on ray 1 or 2 per its 6-index.

### 2d. Skewness deformation `Q`  (`solve_Q`)

For a node Gaussian covariance `M2` and a target node 3-cumulant `Z`, solve

    Z_abc = Q_amn M2_mb M2_nc + Q_bmn M2_ma M2_nc + Q_cmn M2_ma M2_nb,             (2.5)

for `Q` symmetric in its last two indices (M2-eigenbasis diagonalisation; reduces
to the demo2 scalar `Q = z/(3 m^2)`). Used to inject skewness into the MC field
(Sec. 3b). Code: `solve_Q`; forward map `cum3_from_Q`.


## 3. Monte-Carlo engine  (`sachs_mc_core.py`)

### 3a. Exact-propagator discretization

Uniform grid `lambda_k = lambda_min + k dlambda`, `k = 0..N-1`,
`dlambda = (lambda_f - lambda_min)/(N-1)`, `lambda_min = 406`. Discrete update of
(0.1) (response applied INTO node k, then drift + noise injected with unit
response -- exact for the linear part, O(dlambda) F-vertex):

    resp[k] = [D(lambda_{k-1})/D(lambda_k)]^2,  resp[0] = 1,                       (3.1)
    s_k = resp[k] s_{k-1} + (F_abc s_b s_c)|_{k-1} dlambda + (noise)_k,            (3.2)
    kappa = sum_k w_k s_k,  w_k = dlambda (1/2 dlambda at k=0, N-1)  [trapezoid].  (3.3)

Code: `_step` (the `resp_k * s + fv*dlam + f*dlam` line), `_precompute`,
`_field_precompute`. Telescoping (3.2) with F=0 reproduces
`s_k = sum_{j<=k} [D(lambda_j)/D(lambda_k)]^2 dW_j`, i.e. the exact response (3.1)
of Sec. 1, so the Gaussian MC reproduces (2.3) by construction.

### 3b. Noise models

WHITE (`_precompute`, used by `simulate`, `simulate_ff_crn`): the per-step
increment `dW_k` is Gaussian with `Cov(dW_k) = Sigma2_6x6(cos gamma; lambda_k)
dlambda`; `(noise)_k = L_k xi_k`, `L_k L_k^T = Sigma2_6x6 dlambda`.

COLORED AR(1) (`_field_precompute`, used by `simulate_crn`, `simulate_fk_pathA`,
`simulate_fk_vr`): a field `z_k` with stationary node covariance and correlation
length `sigma_lambda`,

    V_k = Sigma2_6x6(cos gamma; lambda_k)/(2 sigma_lambda),  rho = exp(-dlambda/sigma_lambda),  (3.4)
    z_k = rho z_{k-1} + sqrt(1-rho^2) L_k xi_k,  L_k L_k^T = V_k,                  (3.5)

injected into (3.2) as `(noise)_k = f_k dlambda`. The normalization `V_k = Sigma2
/(2 sigma)` enforces `int <z(lambda) z(lambda')> dDelta = 2 sigma_lambda V = Sigma2`
(one-delta collapse), so the colored field reproduces the white density (2.1) as
`sigma_lambda -> 0`.

SKEW INJECTION (`skew_scale != 0`): deform the colored field by a local quadratic,

    f_s,k = z_k + (1/2) Q_k (z_k (x) z_k - V_k),  Q_k = skew_scale * solve_Q(V_k, Zeta6_k/(2 sigma_lambda)^2).  (3.6)

The node target `mu3_k = Zeta6_k/(2 sigma_lambda)^2` (two-delta collapse) gives
`int int <f_s f_s f_s> dDelta_1 dDelta_2 = Zeta6_k`, reproducing the equal-time
density (2.4). `Q_k ~ Zeta6/(3 Sigma2^2) ~ -1.4e7` is SIGMA-INDEPENDENT; the field
skewness `~ Zeta6_node/V^1.5 ~ 4-6` grows as `sigma_lambda^{-0.5}` (Levy).

### 3c. `MCConfig`

`n_real, batch_size, n_lambda, lam_min=406, lam_source=2313.029, seed,
use_f_vertex, apply_anchor (=apply_c0 on Sigma2), skew_scale, sigma_lambda, blowup`.


## 4. Channels and estimators

The MC yields only FULL (all-orders-in-F) correlations, so each is compared to the
sft-wick partial sum truncated at order 2.

### 4a. Gaussian full `<kappa kappa>`  (`simulate` -> `MCResult`)

White noise, F on, no skew. `xi = xi6[0:3,3:6]` with
`xi6 = (1/N_good) sum_real kappa (x) kappa`. Contains O0 + FF + FFF + ... .
`o0_analytic = order0_mc` (2.3) is the exact linear target.

### 4b. FF channel via F-toggle CRN  (`simulate_ff_crn` -> `FFResult`)

Two arms share the SAME white noise; arm A has F on, arm B has F off. Per
realization difference of products (the O0 floor cancels):

    FF_moment_ab = <kappa_a kappa_b>_{F-on} - <kappa_a kappa_b>_{F-off}   (all orders of F)  (4.1)
    FF_disc_ab   = <kappa_a>_{F-on} <kappa_b>_{F-on} - <kappa_a>_{F-off} <kappa_b>_{F-off}     (4.2)
    FF_conn_ab   = FF_moment_ab - FF_disc_ab.                                                  (4.3)

`<kappa>_{F-off} = 0` (linear response of a zero-mean field), so `FF_disc = <kappa>_2
(x) <kappa>_2`, the square of the mean 2nd-order convergence
`<kappa>_2 = int W F_abc <s_b s_c>`. Fields: `ff_moment, ff_conn, ff_disc,
ff_moment_err, mean_kappa_on`.

### 4c. Analytic equal-time FK  (`fk_analytic.py` -> `fk_analytic(cos_gamma, lambda_f)`)

The Order-2 FK diagram (one F-vertex + one zeta) has 4 R-propagators, 0 C: obs_a
--R--> F-vertex ==R,R==> two zeta legs; third zeta leg --R--> obs_b. The two F-side
zeta legs are CO-LOCATED on ray n_a (so the cosine-triple is `(1, cos gamma,
cos gamma)`); the equal-time delta of (2.4) collapses the demo2 4-D time integral
to 2-D, and `integrate_over='all'` turns each external R into the window `W`:

    xi_FK_ab(gamma) = int_0^{lambda_f} dtau int_0^{tau} dlambda_v
        W(tau) W(lambda_v) [D(lambda_v)/D(tau)]^4 G_ab(lambda_v; gamma),          (4.4)
    G_ab = (M + M^T)_ab,  M_ab = sum_{m,n} F_amn Z_bmn,                            (4.5)
    Z = zeta_tensor(n_1, n_1, n_2; lambda_v)  [cosine-triple (1, cos gamma, cos gamma)].

(For kk: `G_00 = -2 (zeta_TTT + zeta_Bmod)`.) The MSR/combinatorial prefactor is 1
(calibrated on the validated demo2 `_fk_spatial_integral`). The external index `a`
rides the F-vertex (`F_amn`), only `b` rides one zeta leg (`Z_bmn`); `M + M^T` is
the x<->y symmetrisation. Code reduces (4.4) to two 1-D quadratures via
`u = D/D(lambda_f)`: `xi_FK_ab = int dlambda_v W(lambda_v) u^4 G_ab(lambda_v)
H(lambda_v)`, `H(lambda_v) = int_{lambda_v}^{lambda_f} W(tau) u(tau)^{-4} dtau`.

### 4d. FK Monte-Carlo, variance-reduced  (`simulate_fk_vr` -> `FKVRResult`)

FK is linear in `Q`, and obs_b connects to its zeta leg F-FREE, so the no-F
windowed response `kappa_d^(0)_B = int W delta_f_B` suffices (no `s_d` ODE). Pull
`Q` OUTSIDE the average:

    T_AB = (1/2) dlambda sum_k W(lambda_k) sum_ij Q_k[B,i,j] C_{A,ij}(k),          (4.6)
    C_{A,ij}(k) = Cov( kappa_g^F_A , z_i(lambda_k) z_j(lambda_k) )  [SAMPLE covariance],
    FK_ab = (T + T^T)[a, 3+b]   (cross-ray block).                                 (4.7)

Three ingredients (all required; see `FK_NOTES.md`):
1. `Q_k` (deterministic) multiplies the AVERAGED `C`, never the per-sample `w_ij =
   z_i z_j - V` -- so the large `Q ~ -1.4e7` does not amplify per-sample variance.
2. `C` is the MATCHED SAMPLE covariance `<kg z_i z_j>_batch - <kg>_batch <z_i z_j>
   _batch` (NOT minus the theoretical `V_k`): the disconnected piece is ~1e7x the
   connected signal, so the subtraction must use the same-sample `<z_i z_j>`.
3. `kappa_g^F = kappa_g(F-on) - kappa_g(F-off)` (shared-noise CRN): only the
   F-induced part L2+ couples to `z z`; dropping the linear L1 cuts the variance
   ~1000x (and fixes the sign). Both arms evolved in the forward loop.

`simulate_fk_pathA` (`PathAResult`) is the equivalent NAIVE estimator
`FK = <kappa_g (x) kappa_d + kappa_d (x) kappa_g>` with `kappa_d` the linearised-F
response to `delta_f` (3.6); blow-up-free but statistics-hungry -- kept for
reference, do not trust its number.


## 5. Comparison with analysis-3

Analysis-3 npz: `sftwick_outputs/2PCF/{C_corr_op_O0, C_corr_op_K_limber_FF,
C_corr_op_K_limber_FK}/xi_*.npz`. The kk channel = rows with `a==0 & b==0`;
`gamma = degrees(arccos(x_hat . y_hat)) * 60` arcmin; `lambda_f = 2313.029`. All
inputs (Sigma2, zeta, F, R, observable) are byte-for-byte the SAME objects the MC
uses, so the cross-check is honest.

| Channel | MC quantity | analysis-3 | Result |
|---|---|---|---|
| Order 0 | `order0_mc` (2.3) | `C_corr_op_O0` | ratio = 0.881 (gamma-indep; `ANCHOR_C0` absorbs it; the 12% is the equal-time collapse dropping the cross-shell tail) |
| FF (order 2) | `simulate_ff_crn` (4.1-4.3) | `C_corr_op_K_limber_FF` | MC FF_moment tracks analysis-3 FF at all gamma; analysis-3's flat large-gamma floor (5.785e-7) = the DISCONNECTED `<kappa>_2^2` (MC FF_disc = 6.1e-7, ratio 1.06); MC FF_conn DECAYS and overlays analysis-3 (FF - floor). FF ~ Sigma2^2 (O0 ~ Sigma2^1), so one `ANCHOR_C0` cannot overlay both. |
| FK (order 2) | `fk_analytic` (4.4) | `C_corr_op_K_limber_FK` | ratio = 0.992, FLAT across the FULL (3,3) matrix (kk peak 5.449e-5 vs 5.489e-5; parity-zero pairs (0,2),(1,2) exactly 0; FK(1,1)=FK(2,2)). The 0.8% is sft-wick's coarser internal discretization (stable under refinement). |
| FK (order 2) | `simulate_fk_vr` (4.6) | via `fk_analytic` | MC FK / fk_analytic ~ 1.04-1.23 (flat over gamma 1'..114', SE ~few%), SIGN correct. Independent stochastic-ODE confirmation of the FK channel to ~10-25%. Grid convergence is NON-monotonic (fine-grid heavy-tailed upward bias); a clean <10% sigma->0 limit is future work. |

Figures (`figures/`): `xi_kappa_channels.pdf` (full O0/FF/FK decomposition, MC vs
analysis-3), `fk_matrix_validation.pdf` (FK (3,3) matrix + selection rules),
`fk_mc_vs_analytic.pdf` (FK MC gamma-sweep), `ff_channel_mc_vs_analysis3.pdf` (FF
decomposition), `anchor_FK_check.pdf`, `anchor_O0_check.pdf`. (The earlier
apples/oranges `mc_vs_analysis3_xi_kappa` is archived under `figures/_archive_*`;
it was made by `compare_to_sft.py`, now superseded by `fig_xi_channels.py`.) Detail
notes: `FF_NOTES.md` (Task 1 verdict), `FK_NOTES.md` (FK anchor + MC story),
`DESIGN.md` (architecture).


## Consistency (doc <-> code, verified 2026-06-06)

- (0.3) F-vertex: `sachs_mc_core._F` and `_f_vertex` give exactly
  `(-(s0^2+s1^2+s2^2), -2 s0 s1, -2 s0 s2)`; matches `inputs/F_tensor.npy`.
- (1) `R, W`: `Background.response_R` = `[D(la')/D(t)]^2 Theta`; `order0_window`
  integrates `[D(la)/D(t)]^2` over `t in [la, lambda_f]`. `theta_sa` = `d lnD/dlam`.
- (2.2) `Sigma2`: `Sigma2Builder.matrix` returns `d/dlam[D^4 C(t,t)]/D^4` (spline
  derivative), corr_op queried equal-time `C_fn(N1,t,n2,t)`. `apply_c0` toggles 1.135.
- (2.3) `O0`: `order0_mc` = `int Sigma2_00 W^2 dlambda`, `W = order0_window`,
  lower limit `lam_lo = 406`.
- (2.5) `solve_Q`: implements (2.5); `cum3_from_Q` round-trips it.
- (3.1-3.3) MC: `_step` = `resp*s + F*dlam + f*dlam`; trapezoid `kap += w*s`,
  `w = 0.5 dlam` at endpoints; `resp[k] = (D[k-1]/D[k])^2`.
- (3.4-3.6) colored/skew: `_field_precompute` sets `V = Sigma2/(2 sigma)`,
  `rho = exp(-dlam/sigma)`, `Q = skew_scale * solve_Q(V, Zeta6/(2 sigma)^2)`;
  `simulate_crn` injects `f_s = z + 0.5 Q (z(x)z - V)`.
- (4.1-4.3) FF: `simulate_ff_crn` accumulates the per-realization product
  difference (moment), the arm means (disconnected), `conn = moment - disc`.
- (4.4-4.5) FK analytic: `fk_analytic` builds `M = einsum('amn,bmn->ab', F, Z)`,
  `G = M + M.T`, `Z = zeta_tensor(N1,N1,n2,lam_v)`, and the swapped 1-D `H`-integral
  with `u = D/D(lambda_f)`.
- (4.6-4.7) FK MC: `simulate_fk_vr` -> `T += 0.5*dlam*W[k]*einsum('Bij,Aij->AB',
  Q[k], Cc_k)`, `Cc_k = C_k - meankap (x) zz_k` (sample), `kg = (kap_g - kap_g0)`
  (F-induced CRN), `FK6 = T + T.T`, returns `FK6[0:3,3:6]`.

NUMERICALLY VERIFIED (2026-06-06, doc formula vs code, independent quadrature):
1. `response_R(t,la') == [D(la')/D(t)]^2`  (exact).
2. `order0_window(la) == int_la^lf [D(la)/D(t)]^2 dt`  (rel 2e-10 vs scipy.quad).
3. `order0_mc(gamma) == int Sigma2_00 W^2 dla`  (rel 2e-4; GL vs trapezoid).
4. `fk_analytic == ` brute-force 2-D quadrature of (4.4)  (rel 2e-2; coarse
   reference grid -- confirms the integrand structure, not a formula mismatch).
5. colored field `2 sigma V_00 == Sigma2_00`, i.e. `int <z z> dDelta == Sigma2`  (exact).
