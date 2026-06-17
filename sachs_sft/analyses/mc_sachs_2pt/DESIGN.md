# DESIGN — Monte-Carlo Sachs evolution for xi_kappa(gamma)

Status: in progress (2026-06-05). Independent MC ground-truth cross-check of
the analytic sft-wick L2 pipeline (analysis 3 O0/FF/FK decomposition).

## 1. What we compute

For each angular separation `gamma`, pick two sky directions `n1, n2`
(`n1.n2 = cos gamma`), evolve the Sachs fluctuation along the two rays under a
shared, correlated random driving field, and Monte-Carlo average to get the
two-point correlation `<s_a(n1, lam_f) s_b(n2, lam_f)>`.  Component map (matches
the sft-wick output npz `a,b in {0,1,2}`):

    s_0 = kappa (Phi_00 / Ricci focus),  s_1 = gamma_+ (Re Psi_0),  s_2 = gamma_x (Im Psi_0)

Headline observable: `xi_kappa(gamma) = <s_0(n1,lam_f) s_0(n2,lam_f)>` (the sign
in `kappa = -int X` cancels in the auto-correlation).

## 2. Dynamics (optical-scalar fluctuation = the user's "linear Jacobi D")

The Sachs fluctuation 3-vector `s = (s_0, s_1, s_2)` obeys

    ds_a/dlambda = -2 theta_sa(lambda) s_a + F_abc s_b s_c + f_a(n, lambda),   s_a(0)=0,

with `theta_sa(lambda) = D'/D = d ln D/dlambda`, `D = a*chi` (closed form,
reuse `D_callable` via `background.Background`), and the paper's F-vertex
(Table 1, `sections/sachs_dynamics.tex`):

    F_111 = -1, F_122 = -1, F_133 = -1, F_212 = -2, F_313 = -2  (others 0).

**Why this is exactly the user's "Linear Jacobi 2x2 D".** The 3-vector `s` holds
the 3 independent components of the symmetric 2x2 distortion matrix; the
optical-scalar form and the `D''=T D` Jacobi form are the same object. We use
the `s`-parameterization because:
- Its linear (F=0) Green's function is **identically** the project's
  `response_R(t,lam') = Theta(t-lam')[D(lam')/D(t)]^2` (since
  `d/dt[D(lam')/D(t)]^2 = -2 theta_sa(t)[.]`), so the MC reproduces the analytic
  pipeline's response by construction.
- Its components map 1:1 onto the sft-wick output npz (`a,b in {0,1,2}`).
The F-vertex is tiny in the weak-lensing regime, so blow-ups are absent in
practice; we keep demo2's `|s|>1e6` guard as a safety net.

Readout per realization: `s_a(n_i, lam_f)` directly (no post-processing of a
distortion matrix needed).

## 3. Driving-field statistics (reuse project cumulants; do NOT use corr_op)

The corr_op C_propagator is the **R-windowed** CorrOp
`= int int [D(lam')/D(lam1)]^2 kappa2(lam',lam'') [D(lam'')/D(lam2)]^2` — i.e.
it is the Order-0 *observable*, the analysis-3 comparison target.  The MC must
**not** use it as the driving noise (that squares the response).  The MC samples
the *bare* driver cumulants and lets the ODE Green's function supply the window
exactly once.  Consistency identity:

    <s_a(lam_f) s_b(lam_f)>_linear
      = int int [D(lam')/D(lam_f)]^2 kappa2_ab(lam',lam'') [D(lam'')/D(lam_f)]^2
      = CorrOp_ab(lam_f, lam_f)   (= analysis-3 Order-0).

**Equal-time (white-noise) model is exact here.** Production is the L_cut=0
sigma2-only path: the 2-cumulant is equal-time (delta-correlated in lambda).
So both the 2- and 3-cumulant drivers are white in lambda, with per-step
increment cumulants `(Sigma2 * dlambda, zeta * dlambda)`.  The small finite
lambda-correlation length is the regularization knob `sigma_lambda`; we verify
`sigma_lambda`-insensitivity (white-noise limit).

### 2-cumulant Sigma2_ab(cos; lambda)  [equal-time density, bare]

Source: the **limber C table**
`callables/C_propagator/limber/limber_C_table_zs5_21pt_stf_fid_omega03161_h06711.npz`
(`cl_tensor_convention="bare"`, `delta_measure="lambda"`, delta-collapsed in
lambda).  Read the equal-lambda diagonal `cl_pp[:,j,j]`, `cl_px[:,j,j]`,
`cl_xx[:,j,j]` (Phi-Phi, Phi-Psi, Psi-Psi angular power at shell `lam_j`) and
build the (3,3) angular covariance via the great-circle spin-2 Legendre/Wigner-d
sum (reuse `spin_rotation` from `scripts/`):

    Sigma2_00 = sum_l (2l+1)/(4pi) P_l(cos) cl_pp[l,j]
    Sigma2 spin-2 block (1,1),(1,2),(2,2) via d^l_{2,2}, d^l_{2,-2} and psi-rotation
    Sigma2 cross (0,1),(0,2) via d^l_{0,2}

evaluated at `cos = 1` (within-ray) and `cos = cos gamma` (cross-ray).

**Numerical anchor (gates everything):** the analytic linear order
`O0_MC(gamma) = int_0^{lam_f} [D(lam')/D(lam_f)]^4 Sigma2_00(cos gamma; lam') dlam'`
must reproduce analysis-3 `C_corr_op_O0` kk within a few %.  If a constant
measure factor appears, it is identified and recorded here (expected O(1)).

### 3-cumulant zeta_abc(cos-triple; lambda)  [equal-time density, bare]

Source: `callables/kappa3_vertex/equal_time_limber/` `coupling_fn(n_list,t_list)
-> (3,3,3)` (`ALREADY_R_CONTRACTED=False`, `EQUAL_TIME=True`,
`RADIAL_MEASURE="lambda"`, (1+z)^4/leg baked in).  For the 2-ray bundle the
needed direction triples are `(n1n1n1),(n1n1n2),(n1n2n2),(n2n2n2)` with
`cos11=1, cos12=cos gamma`.

### Skewness injection (tensor generalization of demo2)

demo2 uses `eta_tilde = eta + alpha(eta^2 - lam)`; its 3rd cumulant is
`kappa3 = 2 alpha lam^2 delta_ab delta_bc [...]` (`demo2/k3_coupling.py`).  We
invert: given the target `zeta` and the per-step Gaussian covariance
`M2 = Sigma2 * dlambda`, solve the closed-form linear map for the deformation
tensor `Q` (symmetric in its last two indices) such that

    f_tilde = f + (1/2) Q (f ox f - <f ox f>)
    => Cum3[f_tilde]_abc = Q_amn M2_mb M2_nc + Q_bmn M2_ma M2_nc + Q_cmn M2_ma M2_nb
       (set equal to zeta_abc * dlambda; solve for Q).

Done per (gamma, lambda) on the 6-dim two-ray space (3 comps x 2 rays).  Note
the white-noise skewness is finite only for finite step `dlambda = sigma_lambda`
(the Levy limit has divergent per-step skewness); this is the regularization.

## 4. Engines (both)

Shared NumPy `driver_stats` precompute (Sigma2 6x6 cov + Q over the lambda grid)
feeds:
- `sachs_mc_core.py` — NumPy fixed-grid Heun, realizations as array columns,
  alpha-toggle common-random-numbers, blow-up guard (demo2-parity).
- `sachs_mc_jax.py` — JAX + diffrax, vmap over realizations, CRN via shared key
  (needs `pip install diffrax` in the sft-wick env; jax 0.10 present).

## 5. Comparison (analysis 3) and the two MC runs

The MC yields only the full (all-orders) correlation. Run twice:

| MC driver | full output | analysis-3 target (order <= 2) |
|---|---|---|
| Gaussian (Sigma2 only) | xi^G(gamma) | O0 + FF  (`C_corr_op_O0` + `C_corr_op_K_limber_FF`, a=b=0) |
| skewed (Sigma2 + zeta) | xi^NG(gamma) | O0 + FF + FK (+ `C_corr_op_K_limber_FK`) |

CRN difference `xi^NG - xi^G` (shared base Gaussian draws, only skewness
toggled) isolates the 3-cumulant effect at low variance; compare to analysis-3
FK (order 2).  Overlay on `analyses/analysis3` curves; gamma grid = 40 log pts
0.5..5000 arcmin, `gamma = degrees(arccos(x.y))*60`, z_s=5 (`lam_f`).

## 6. Validation ladder (each gates the next)

1. Background: `D=a*chi`, `theta_sa=D'/D`, spline response == `D_callable.response_R`. (DONE)
2. Linear-response / O0 anchor: `int [D/D]^4 Sigma2_00 dlam` == analysis-3 O0 (few %).
3. Gaussian MC full == O0 + FF.
4. Skewed MC full == O0 + FF + FK; CRN diff == FK (order 2).
5. sigma_lambda-insensitivity (white-noise limit): halve sigma_lambda, xi unchanged.
6. Boundary rule (`~/.claude/rules/common/boundary-validation.md`): extreme gamma
   (0.5', 5000') and lambda corners; promote step 2 + one pinned xi to
   `tests/test_mc_boundaries.py` (isfinite + rel_err < tol).

## 7. Files

`background.py` (done), `driver_stats.py`, `sachs_mc_core.py`,
`sachs_mc_jax.py`, `response_check.py`, `run_mc.py`, `compare_to_sft.py`,
`tests/test_mc_boundaries.py`.

## 8. Key reuse paths

- `scripts/sachs_sft/scripts/D_callable.py` (D, theta_sa, response_R).
- `scripts/sachs_sft/scripts/spin_rotation.py` (great-circle spin-2 basis, Wigner-d).
- `callables/C_propagator/limber/limber_C_table_*.npz` (bare equal-time Sigma2).
- `callables/kappa3_vertex/equal_time_limber/` (bare equal-time zeta).
- `sftwick_outputs/2PCF/{C_corr_op_O0,C_corr_op_K_limber_FF,C_corr_op_K_limber_FK}/xi_*.npz` (analysis-3 targets).
- `sft-wick/examples/demo2/{run_simulation,k3_coupling}.py` (MC template + zeta<->deform relation).
