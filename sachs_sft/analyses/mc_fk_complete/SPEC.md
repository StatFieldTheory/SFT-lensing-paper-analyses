# Specification of the stochastic system and the FK observable

*This file is the mathematical statement only. It deliberately contains no
estimator design and no implementation, so that a clean-room re-derivation can
be checked against `fk_expect_exact.py` / `fk_complete_core.py` without having
read them.*

## 1. The stochastic system

State `s[k] in R^6` on a uniform grid `lam[k] = lam_min + k dlam`,
`k = 0 .. N-1`, `lam_min = 406.0`, `lam[N-1] = LAM_SOURCE_BASELINE = 2313.029`.
The 6 components are two rays x three screen components
`(Phi00, Re Psi0, Im Psi0)`; ray 0 is indices 0..2, ray 1 is 3..5.

Background: `D(lam)` from `background.Background(lam_min=120, lam_max=lam_f+5)`.
Response ratios `resp[0] = 1`, `resp[k] = (D[k-1]/D[k])^2`.

Evolution (explicit Euler, F evaluated at the OLD state):

```
s[k] = resp[k] s[k-1] + F(s[k-1]) dlam + f[k] dlam ,     s[-1] = 0
```

`F` acts within each ray on its 3-vector, with the paper's vertex
(`sections/sachs_dynamics.tex`, Table 1, 0-indexed)
`F_111 = F_122 = F_133 = -1`, `F_212 = F_313 = -2`, others zero; i.e.
`F(s)_0 = -(s0^2+s1^2+s2^2)`, `F(s)_1 = -2 s0 s1`, `F(s)_2 = -2 s0 s2`.

Observable: `kappa_a = sum_k w_k s_a[k]`, `w_k = dlam` with `w_0 = w_{N-1} =
dlam/2`.

## 2. The driving field

`f[k] = z[k] + 0.5 Q[k] : (z[k] z[k]^T - C[k])` where

* `z` is a zero-mean Gaussian AR(1) chain:
  `z[0] = L[0] xi_0`, `z[k] = rho z[k-1] + sqrt(1-rho^2) L[k] xi_k`,
  `rho = exp(-dlam/sigma_lambda)`, `L[k] L[k]^T = V[k] = Sigma2(lam_k)/(2 sigma_lambda)`,
  `xi_k` iid standard normal 6-vectors.
* `Sigma2(cos gamma; lam)` is the two-ray (6,6) second-cumulant DENSITY,
  assembled from `driver_stats.Sigma2Builder` (within-ray block on the diagonal,
  cross-ray block off-diagonal).
* `Q[k]` is symmetric in its last two indices and is fixed by requiring the field
  to carry the tabulated three-point cumulant density
  `zeta6(cos gamma; lam)` (`driver_stats.zeta6`, a (6,6,6) tensor). The relation
  solved by `driver_stats.solve_Q(M, Z)` is
  `Z_abc = Q_amn M_mb M_nc + Q_bmn M_ma M_nc + Q_cmn M_ma M_nb`.
* `C[k]` is a subtracted constant (any choice shifts the drive by a constant).

The `sigma_lambda > 0` colored-noise regulator is unavoidable: in the white-noise
limit the deformation per step does not scale down with `dlam` and the
accumulated skew drive diverges.

## 3. The observable to be tested

`FK` is the contribution to the cross-ray block `<kappa_a(n1) kappa_b(n2)>`
at first order in the F vertex AND first order in `zeta`. For the
kappa-kappa entry `(a,b) = (0,0)` the analytic target is the paper's converged
permutation-closed fold,
`SFT-lensing-paper-analyses/sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz`
(order 2, a = b = 0):

| gamma | 0.5' | 1.0' | 2.0' | 5.0' | 17.3' |
|---|---|---|---|---|---|
| FK | +1.948e-05 | +1.598e-05 | +1.223e-05 | +7.038e-06 | +1.966e-06 |

An independent continuum quadrature of the same diagram is in
`../rebuild/fk_kernel_crosscheck.py`:

```
FK(gamma) = 2 int dv  W(v) H(v) [ - sum_c zeta6_{c c 3}(gamma; v) ],
W(la) = int_la^{lam_f} [D(la)/D(t)]^2 dt,
H(v)  = D(v)^4 int_v^{lam_f} W(u)/D(u)^4 du.
```

## 4. Environment

```
export PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
```
Bind the corrected vertex callable exactly as `_bootstrap.wire()` does:
`perm_aware_kappa3_callable` with
`TABLE_PATH = callables/kappa3_vertex/equal_time_limber_cut15360_permaware/table_permclosed.npz`, `_CACHE = None`,
then `driver_stats._k3 = pa`.
