"""Placement-complete, variance-controlled Monte-Carlo of the FK channel.

What is different from ``sachs_mc_core.simulate_fk_vr``
------------------------------------------------------
The three-point cumulant is realised by a local quadratic deformation of the
Gaussian base field, ``f = z + 0.5 Q (z z - <z z>)``.  At O(Q) the induced
three-point function is NON-LOCAL along the ray,

    <f_b[m] f_c[n] f_B[p]> = Q[m]_bij R[m,n]_ic R[m,p]_jB + (Q at n) + (Q at p),

with ``R[k,l] = <z[k] z[l]^T>`` -- three terms, one per leg that carries the
``Q``, each collapsing to ``Q:Sigma2 Sigma2`` only as ``sigma_lambda -> 0``.  In
the FK diagram those three legs go to the far observable and to the two inputs
of the F vertex.  ``simulate_fk_vr`` injects the deformation only on
the observable leg and therefore measures one of the three; the missing two are
recovered here by letting the skew drive propagate through the LINEARISED F
vertex, which is what ``simulate_fk_pathA`` does.

What is different from ``simulate_fk_pathA``
-------------------------------------------
Path A forms ``kappa_g * kappa_d`` per realisation.  Both factors are dominated
by their F-free parts, whose product ``kappa1 * kappa_d1`` has zero mean but
carries the full ``Q^2`` spread -- that is the variance pathology.  Splitting the
estimator into the two placement classes,

    T1_AB = Cov( kappa_g^(2)_A ,  kappa_d^(1)_B )   -- Q on the observable leg
    T2_AB = Cov( kappa_g^(1)_A ,  kappa_d^(2)_B )   -- Q on an F-input leg
    FK    = (T1 + T2) + (T1 + T2)^T

never forms that product: each term pairs an F-free factor with an F-induced
one.  The huge ``Q`` scales signal and noise together, so the relative error is
``Q``-independent.  Two flavours are returned:

* ``pert``  -- strictly order-controlled arms (``kappa^(1)``, ``kappa^(2)``,
  ``kappa_d^(1)``, ``kappa_d^(2)`` evolved separately).  Exactly the F^1 zeta^1
  diagram, nothing else.
* ``full``  -- the nonlinear Sachs arm with the increment field linearised about
  it (Path A's dynamics), split by the F-free control variates
  ``kappa_g - kappa_g0`` and ``kappa_d - kappa_d0``.  Carries F^2 and higher, so
  the difference from ``pert`` measures the higher-order contamination.

Centering
---------
The drive uses the BATCH SAMPLE mean of ``z_u z_v`` at each node, not the
nominal ``V``.  For ``T1`` any constant is removed by the covariance; for ``T2``
it is not -- there the two ``Q`` legs can contract with each other, and that
tadpole survives unless the subtracted constant is the field's actual variance.
The AR(1) chain is not stationary (``V`` varies along the ray), so the nominal
``V`` is the wrong constant; the sample mean is the right one and additionally
removes its own sampling noise.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import fk_expect_exact as ex

_CORE = ex._CORE


@dataclass
class FKCompleteResult:
    gamma_arcmin: float
    cos_gamma: float
    sigma_lambda: float
    n_lambda: int
    n_real: int
    seed: int
    t1: NDArray[np.float64]          # (6,6) perturbative, before symmetrisation
    t2: NDArray[np.float64]          # (6,6)
    t1_full: NDArray[np.float64]     # (6,6) nonlinear-arm flavour
    t2_full: NDArray[np.float64]
    fk: NDArray[np.float64]          # (3,3) cross-ray, perturbative
    fk_full: NDArray[np.float64]     # (3,3) cross-ray, nonlinear arm
    fk_vr_like: NDArray[np.float64]  # (3,3) T1 only -- the published structure
    fk_batch_se: NDArray[np.float64]  # (3,3) batch-scatter SE, perturbative
    n_good: int
    n_batches: int

    @property
    def kk(self) -> float:
        return float(self.fk[0, 0])


def _f_vertex6(s: NDArray[np.float64]) -> NDArray[np.float64]:
    """F_abc s_b s_c for the stacked two-ray state, s shaped (6, m)."""
    out = np.empty_like(s)
    for r in (0, 3):
        s0, s1, s2 = s[r], s[r + 1], s[r + 2]
        out[r] = -(s0 * s0 + s1 * s1 + s2 * s2)
        out[r + 1] = -2.0 * s0 * s1
        out[r + 2] = -2.0 * s0 * s2
    return out


def _f_vertex6_lin(sg: NDArray[np.float64],
                   sd: NDArray[np.float64]) -> NDArray[np.float64]:
    """DF(sg).sd = F_abc (sg_b sd_c + sd_b sg_c), stacked two-ray, (6, m)."""
    out = np.empty_like(sd)
    for r in (0, 3):
        g0, g1, g2 = sg[r], sg[r + 1], sg[r + 2]
        d0, d1, d2 = sd[r], sd[r + 1], sd[r + 2]
        out[r] = -2.0 * (g0 * d0 + g1 * d1 + g2 * d2)
        out[r + 1] = -2.0 * (g0 * d1 + d0 * g1)
        out[r + 2] = -2.0 * (g0 * d2 + d0 * g2)
    return out


def simulate_fk_complete(gamma_arcmin: float, sigma_lambda: float = 8.0,
                         n_lambda: int = 1000, n_real: int = 24_000,
                         batch_size: int = 2_000, seed: int = 12345,
                         lam_min: float = 406.0, grid=None, nodes=None,
                         calibrate: str = "smeared",
                         blowup: float = 1e6) -> FKCompleteResult:
    cos_gamma = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    grid = grid if grid is not None else ex.build_grid(n_lambda=n_lambda,
                                                       lam_min=lam_min)
    V, A, Q, rho = nodes if nodes is not None else ex.node_stats(
        grid, cos_gamma, sigma_lambda, calibrate=calibrate)
    N = grid.lam.size
    dlam, resp, Wd = grid.dlam, grid.resp, grid.Wd
    Lf = np.array([_CORE._cholesky_psd(V[k]) for k in range(N)])
    Qflat = Q.reshape(N, 6, 36)                       # (N, a, uv) for the drive
    sqrt1mr2 = float(np.sqrt(max(1.0 - rho * rho, 0.0)))
    w_node = np.full(N, dlam)
    w_node[0] = w_node[-1] = 0.5 * dlam

    rng = np.random.default_rng(seed)
    n_batches = max(1, n_real // batch_size)
    t1_b, t2_b, t1f_b, t2f_b = [], [], [], []
    n_good = 0

    for _ in range(n_batches):
        m = batch_size
        # --- pass 1: the AR(1) base field, stored ---------------------------
        zstore = np.empty((N, 6, m))
        z = np.zeros((6, m))
        for k in range(N):
            innov = Lf[k] @ rng.standard_normal((6, m))
            z = innov if k == 0 else (rho * z + sqrt1mr2 * innov)
            zstore[k] = z
        # per-node sample second moment (the centering constant, see module doc)
        c_node = np.einsum("kim,kjm->kij", zstore, zstore, optimize=True) / m

        # --- pass 2: the arms ----------------------------------------------
        s1 = np.zeros((6, m)); s2 = np.zeros((6, m))
        sd1 = np.zeros((6, m)); sd2 = np.zeros((6, m))
        sg = np.zeros((6, m)); sd = np.zeros((6, m))
        k1 = np.zeros((6, m)); k2 = np.zeros((6, m))
        kd1 = np.zeros((6, m)); kd2 = np.zeros((6, m))
        kg = np.zeros((6, m)); kd = np.zeros((6, m))
        blown = np.zeros(m, dtype=bool)
        for k in range(N):
            zk = zstore[k]
            ww = np.einsum("im,jm->ijm", zk, zk, optimize=True) - c_node[k][:, :, None]
            df = 0.5 * (Qflat[k] @ ww.reshape(36, m))          # (6, m)
            # strictly order-controlled arms
            fv1 = _f_vertex6(s1)
            fvx = _f_vertex6_lin(s1, sd1)
            s2n = resp[k] * s2 + fv1 * dlam
            sd2n = resp[k] * sd2 + fvx * dlam
            s1n = resp[k] * s1 + zk * dlam
            sd1n = resp[k] * sd1 + df * dlam
            # nonlinear arm + increment field linearised about it (Path A)
            fvg = _f_vertex6(sg)
            fvd = _f_vertex6_lin(sg, sd)
            sgn = resp[k] * sg + fvg * dlam + zk * dlam
            sdn = resp[k] * sd + fvd * dlam + df * dlam
            s1, s2, sd1, sd2, sg, sd = s1n, s2n, sd1n, sd2n, sgn, sdn
            wk = w_node[k]
            k1 += wk * s1; k2 += wk * s2
            kd1 += wk * sd1; kd2 += wk * sd2
            kg += wk * sg; kd += wk * sd
            blown |= np.any(np.abs(sg) > blowup, axis=0)
        good = ~blown
        mb = int(good.sum())

        def cov(x, y):
            xg = x[:, good]; yg = y[:, good]
            return (xg @ yg.T) / mb - np.outer(xg.mean(axis=1), yg.mean(axis=1))

        t1_b.append(cov(k2, kd1))
        t2_b.append(cov(k1, kd2))
        t1f_b.append(cov(kg - k1, kd1))
        t2f_b.append(cov(k1, kd - kd1))
        n_good += mb

    t1 = np.mean(t1_b, axis=0); t2 = np.mean(t2_b, axis=0)
    t1f = np.mean(t1f_b, axis=0); t2f = np.mean(t2f_b, axis=0)
    tot = np.array(t1_b) + np.array(t2_b)
    fk6_b = tot + np.transpose(tot, (0, 2, 1))
    fk6 = (t1 + t2) + (t1 + t2).T
    fk6_full = (t1f + t2f) + (t1f + t2f).T
    vr6 = t1 + t1.T
    se6 = fk6_b.std(axis=0, ddof=1) / np.sqrt(len(t1_b)) if len(t1_b) > 1 \
        else np.zeros((6, 6))
    return FKCompleteResult(
        gamma_arcmin, cos_gamma, float(sigma_lambda), N, n_real, seed,
        t1, t2, t1f, t2f, fk6[0:3, 3:6], fk6_full[0:3, 3:6], vr6[0:3, 3:6],
        se6[0:3, 3:6], n_good, len(t1_b))
