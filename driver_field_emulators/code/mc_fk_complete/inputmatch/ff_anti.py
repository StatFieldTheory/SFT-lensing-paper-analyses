"""Antithetic-pair variance reduction for the FF CRN estimator.

This is NOT a re-implementation of the physics: it calls ``sachs_mc_core``'s own
``_precompute``, ``_f_vertex`` and the same Euler recursion, and only changes how
the standard normals are drawn -- realisation j and realisation j + m/2 of every
batch share the SAME z at every lambda node, with opposite sign.

Why it helps.  Write the F-on state as s = s0 + s1 + s2 + ... with s_n carrying
n powers of the vertex, so s_n has noise-parity (-1)^(n+1).  The CRN estimand is

    on.on - off.off = 2 s0 s1 + s1 s1 + 2 s0 s2 + O(F^3),

and <2 s0 s1> = 0 (odd Gaussian moment): the leading term is pure noise, and it
is exactly what z -> -z flips.  Averaging the antithetic pair cancels it
identically per realisation, leaving the signal terms untouched.  The estimator
is unbiased, and the pair average is the same random variable the production
estimator averages, so the two must agree within errors (checked in t_anti.py).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np


def ff_moment_anti(mc, cfg, gamma_arcmin: float):
    """Return (ff_moment_kk, se_kk, n_pairs) for the (0,0) two-ray kk element."""
    cfg = replace(cfg, skew_scale=0.0, sigma_lambda=0.0)
    cos_gamma = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    bg, lam, dlam, resp, L_white, _Q, _sig, builder = mc._precompute(cfg, cos_gamma)
    rng = np.random.default_rng(cfg.seed)
    N = cfg.n_lambda
    m = int(cfg.batch_size)
    if m % 2:
        raise ValueError("batch_size must be even for antithetic pairing")
    h = m // 2
    n_batches = max(1, cfg.n_real // m)

    tot = 0.0
    tot2 = 0.0
    npair = 0
    nblow = 0
    for _ in range(n_batches):
        s_on = np.zeros((6, m)); s_off = np.zeros((6, m))
        kap_on = np.zeros((6, m)); kap_off = np.zeros((6, m))
        blown = np.zeros(m, dtype=bool)
        for k in range(N):
            zh = rng.standard_normal((6, h))
            z = np.concatenate([zh, -zh], axis=1)          # antithetic pairing
            noise = L_white[k] @ z
            s_cr = s_on.reshape(2, 3, m).transpose(1, 0, 2)
            fv = mc._f_vertex(s_cr).transpose(1, 0, 2).reshape(6, m)
            s_on = resp[k] * s_on + fv * dlam + noise
            s_off = resp[k] * s_off + noise
            w = 0.5 * dlam if (k == 0 or k == N - 1) else dlam
            kap_on += w * s_on; kap_off += w * s_off
            blown |= np.any(np.abs(s_on) > cfg.blowup, axis=0)
        # kk element: ray-0 kappa (index 0) x ray-1 kappa (index 3)
        d = kap_on[0] * kap_on[3] - kap_off[0] * kap_off[3]
        pair = 0.5 * (d[:h] + d[h:])                        # antithetic average
        bad = blown[:h] | blown[h:]
        pair = pair[~bad]
        nblow += int(bad.sum())
        tot += pair.sum(); tot2 += float(pair @ pair); npair += pair.size

    mean = tot / max(npair, 1)
    var = tot2 / max(npair, 1) - mean * mean
    se = float(np.sqrt(max(var, 0.0) / max(npair, 1)))
    return float(mean), se, npair, nblow
