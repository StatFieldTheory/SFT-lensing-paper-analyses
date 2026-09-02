"""Monte-Carlo of the stochastic Sachs system with a defect in the DYNAMICS.

A copy of ``fk_complete_core.simulate_fk_complete`` restricted to the strictly
order-controlled arms, with three switchable errors in the ODE itself:

  ``F_new_state``   F evaluated at the NEW state (Euler off-by-one)
  ``resp_pow1``     response ratio (D_{k-1}/D_k)^2 -> ^1
  ``F_no_leg_sym``  the linearised vertex drops one of its two terms

The point is that these live in the simulated dynamics, not in any reference
kernel, so what they test is the object the paper claims to have validated.
"""
import _env
import numpy as np
import fk_expect_exact as ex
import inj_exact as ij

_CORE = ex._CORE


def _fv(s):
    out = np.empty_like(s)
    for r in (0, 3):
        s0, s1, s2 = s[r], s[r + 1], s[r + 2]
        out[r] = -(s0 * s0 + s1 * s1 + s2 * s2)
        out[r + 1] = -2.0 * s0 * s1
        out[r + 2] = -2.0 * s0 * s2
    return out


def _fvlin(sg, sd, half=False):
    out = np.empty_like(sd)
    for r in (0, 3):
        g0, g1, g2 = sg[r], sg[r + 1], sg[r + 2]
        d0, d1, d2 = sd[r], sd[r + 1], sd[r + 2]
        if half:      # F_abc sg_b sd_c only (drop the sd_b sg_c placement)
            out[r] = -1.0 * (g0 * d0 + g1 * d1 + g2 * d2)
            out[r + 1] = -2.0 * (g0 * d1)
            out[r + 2] = -2.0 * (g0 * d2)
        else:
            out[r] = -2.0 * (g0 * d0 + g1 * d1 + g2 * d2)
            out[r + 1] = -2.0 * (g0 * d1 + d0 * g1)
            out[r + 2] = -2.0 * (g0 * d2 + d0 * g2)
    return out


def simulate(cos_gamma, grid, nodes, sigma_lambda=8.0, n_real=16000,
             batch_size=2000, seed=12345, defect="none"):
    V, A, Q, rho = nodes
    N = grid.lam.size
    dlam, resp, Wd = grid.dlam, grid.resp, grid.Wd
    Lf = np.array([_CORE._cholesky_psd(V[k]) for k in range(N)])
    Qflat = Q.reshape(N, 6, 36)
    s1mr2 = float(np.sqrt(max(1.0 - rho * rho, 0.0)))
    w = np.full(N, dlam); w[0] = w[-1] = 0.5 * dlam
    rng = np.random.default_rng(seed)
    nb = max(1, n_real // batch_size)
    t1b, t2b = [], []
    for _ in range(nb):
        m = batch_size
        zst = np.empty((N, 6, m)); z = np.zeros((6, m))
        for k in range(N):
            innov = Lf[k] @ rng.standard_normal((6, m))
            z = innov if k == 0 else (rho * z + s1mr2 * innov)
            zst[k] = z
        cn = np.einsum("kim,kjm->kij", zst, zst, optimize=True) / m
        s1 = np.zeros((6, m)); s2 = np.zeros((6, m))
        sd1 = np.zeros((6, m)); sd2 = np.zeros((6, m))
        k1 = np.zeros((6, m)); k2 = np.zeros((6, m))
        kd1 = np.zeros((6, m)); kd2 = np.zeros((6, m))
        half = (defect == "F_no_leg_sym")
        for k in range(N):
            zk = zst[k]
            ww = np.einsum("im,jm->ijm", zk, zk, optimize=True) - cn[k][:, :, None]
            df = 0.5 * (Qflat[k] @ ww.reshape(36, m))
            s1n = resp[k] * s1 + zk * dlam
            sd1n = resp[k] * sd1 + df * dlam
            if defect == "F_new_state":
                fv1 = _fv(s1n); fvx = _fvlin(s1n, sd1n, half)
            else:
                fv1 = _fv(s1); fvx = _fvlin(s1, sd1, half)
            s2 = resp[k] * s2 + fv1 * dlam
            sd2 = resp[k] * sd2 + fvx * dlam
            s1, sd1 = s1n, sd1n
            k1 += w[k] * s1; k2 += w[k] * s2
            kd1 += w[k] * sd1; kd2 += w[k] * sd2

        def cov(x, y):
            return (x @ y.T) / m - np.outer(x.mean(1), y.mean(1))
        t1b.append(cov(k2, kd1)); t2b.append(cov(k1, kd2))
    t1 = np.mean(t1b, 0); t2 = np.mean(t2b, 0)
    return float(((t1 + t2) + (t1 + t2).T)[0, 3])
