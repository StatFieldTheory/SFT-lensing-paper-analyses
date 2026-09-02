"""Instrumented clone of fk_complete_core.simulate_fk_complete.

Byte-for-byte the same arithmetic when ``center="sample"`` and
``rng_mode="int"``; adds switches that isolate the candidate mechanisms:

* ``center``  -- "sample" (production: in-sample batch second moment),
                 "exact"  (the constant A[k] the exact Wick expectation assumes),
                 "pilot"  (an INDEPENDENT pilot batch's second moment: unbiased
                           centering with the same statistical noise level).
* ``rng_mode`` -- "int"   (production: default_rng(seed)),
                  "spawn" (SeedSequence(base).spawn(n)[s], guaranteed independent).
* returns the PER-BATCH t1/t2 and the blow-up counts.
"""
from __future__ import annotations
import sys
import numpy as np

sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/"
                   "driver_field_emulators/code/mc_fk_complete")
import fk_expect_exact as ex                      # noqa: E402
from fk_complete_core import _f_vertex6, _f_vertex6_lin   # noqa: E402

_CORE = ex._CORE


def make_rng(seed=None, base=None, index=None, mode="int", nmax=4096):
    if mode == "int":
        return np.random.default_rng(seed)
    if mode == "spawn":
        return np.random.default_rng(np.random.SeedSequence(base).spawn(nmax)[index])
    raise ValueError(mode)


def run(gamma_arcmin, sigma_lambda=8.0, n_lambda=1000, n_real=24_000,
        batch_size=2_000, seed=12345, grid=None, nodes=None,
        calibrate="smeared", blowup=1e6, center="sample", rng=None,
        c_ext=None, pilot_size=None):
    cos_gamma = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    grid = grid if grid is not None else ex.build_grid(n_lambda=n_lambda)
    V, A, Q, rho = nodes if nodes is not None else ex.node_stats(
        grid, cos_gamma, sigma_lambda, calibrate=calibrate)
    N = grid.lam.size
    dlam, resp = grid.dlam, grid.resp
    Lf = np.array([_CORE._cholesky_psd(V[k]) for k in range(N)])
    Qflat = Q.reshape(N, 6, 36)
    sqrt1mr2 = float(np.sqrt(max(1.0 - rho * rho, 0.0)))
    w_node = np.full(N, dlam)
    w_node[0] = w_node[-1] = 0.5 * dlam

    rng = rng if rng is not None else np.random.default_rng(seed)
    n_batches = max(1, n_real // batch_size)
    t1_b, t2_b, nb_good = [], [], []

    for _ in range(n_batches):
        m = batch_size
        zstore = np.empty((N, 6, m))
        z = np.zeros((6, m))
        for k in range(N):
            innov = Lf[k] @ rng.standard_normal((6, m))
            z = innov if k == 0 else (rho * z + sqrt1mr2 * innov)
            zstore[k] = z
        if center == "sample":
            c_node = np.einsum("kim,kjm->kij", zstore, zstore, optimize=True) / m
        elif center == "exact":
            c_node = A
        elif center == "given":
            c_node = c_ext
        elif center == "pilot":
            mp = pilot_size or m
            zp = np.empty((N, 6, mp)); zz = np.zeros((6, mp))
            for k in range(N):
                iv = Lf[k] @ rng.standard_normal((6, mp))
                zz = iv if k == 0 else (rho * zz + sqrt1mr2 * iv)
                zp[k] = zz
            c_node = np.einsum("kim,kjm->kij", zp, zp, optimize=True) / mp
        else:
            raise ValueError(center)

        s1 = np.zeros((6, m)); s2 = np.zeros((6, m))
        sd1 = np.zeros((6, m)); sd2 = np.zeros((6, m))
        sg = np.zeros((6, m)); sd = np.zeros((6, m))
        k1 = np.zeros((6, m)); k2 = np.zeros((6, m))
        kd1 = np.zeros((6, m)); kd2 = np.zeros((6, m))
        blown = np.zeros(m, dtype=bool)
        for k in range(N):
            zk = zstore[k]
            ww = np.einsum("im,jm->ijm", zk, zk, optimize=True) - c_node[k][:, :, None]
            df = 0.5 * (Qflat[k] @ ww.reshape(36, m))
            fv1 = _f_vertex6(s1)
            fvx = _f_vertex6_lin(s1, sd1)
            s2n = resp[k] * s2 + fv1 * dlam
            sd2n = resp[k] * sd2 + fvx * dlam
            s1n = resp[k] * s1 + zk * dlam
            sd1n = resp[k] * sd1 + df * dlam
            fvg = _f_vertex6(sg)
            fvd = _f_vertex6_lin(sg, sd)
            sgn = resp[k] * sg + fvg * dlam + zk * dlam
            sdn = resp[k] * sd + fvd * dlam + df * dlam
            s1, s2, sd1, sd2, sg, sd = s1n, s2n, sd1n, sd2n, sgn, sdn
            wk = w_node[k]
            k1 += wk * s1; k2 += wk * s2
            kd1 += wk * sd1; kd2 += wk * sd2
            blown |= np.any(np.abs(sg) > blowup, axis=0)
        good = ~blown
        mb = int(good.sum())

        def cov(x, y):
            xg = x[:, good]; yg = y[:, good]
            return (xg @ yg.T) / mb - np.outer(xg.mean(axis=1), yg.mean(axis=1))

        t1_b.append(cov(k2, kd1)); t2_b.append(cov(k1, kd2)); nb_good.append(mb)

    t1_b = np.array(t1_b); t2_b = np.array(t2_b)
    tot = t1_b + t2_b
    fk6_b = tot + np.transpose(tot, (0, 2, 1))
    return dict(kk_batches=fk6_b[:, 0, 3], kk=float(fk6_b[:, 0, 3].mean()),
                t1=float(t1_b.mean(axis=0)[0, 3]),
                t2=float(t2_b.mean(axis=0)[0, 3]),
                n_good=np.array(nb_good), n_batches=len(nb_good))
