"""Continuum limit of the clean-room discrete kernel, with controllable range.

Same diagram as fk_ref_discrete.py, evaluated as

    FK = 2 int_lo^hi dv W(v) H(v) [ -sum_c zeta6_{cc3}(gamma; v) ],
    W(la) = int_la^{lam_f} [D(la)/D(t)]^2 dt,
    H(v)  = D(v)^4 int_v^{lam_f} W(u)/D(u)^4 du,

with W from driver_stats.order0_window and H by accurate Gauss-Legendre
quadrature from v itself (not from the first node of a shared sub-grid).
Used (a) to confirm that the discrete kernels of fk_ref_discrete.py converge
to these, and (b) to isolate the effect of the integration range, which is
the one place ../rebuild/fk_kernel_crosscheck.py differs.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bootstrap import wire, FOLD  # noqa: E402
from fk_ref_discrete import load_fold_kk, zeta_contraction  # noqa: E402


def H_exact(ds, bg, v, lam_f, n_gauss=64, n_panel=8):
    """D(v)^4 int_v^{lam_f} W(u)/D(u)^4 du, panelled Gauss-Legendre from v."""
    nodes, wts = leggauss(n_gauss)
    v = np.asarray(v, float)
    out = np.empty_like(v)
    for i, x in enumerate(v):
        if x >= lam_f:
            out[i] = 0.0
            continue
        edges = np.linspace(x, lam_f, n_panel + 1)
        tot = 0.0
        for e0, e1 in zip(edges[:-1], edges[1:]):
            u = 0.5 * (e1 - e0) * nodes + 0.5 * (e1 + e0)
            jw = 0.5 * (e1 - e0) * wts
            Wu = ds.order0_window(bg, u, lam_f, n_gauss=n_gauss)
            Du = np.asarray(bg.D(u), float)
            tot += float(np.sum(jw * Wu / Du**4))
        out[i] = float(np.asarray(bg.D(x), float)) ** 4 * tot
    return out


def fk_cont(ds, bg, pa, gamma, lam_f, lo, hi, n_out=400, n_gauss=64):
    v = np.linspace(lo, hi, n_out)
    W = ds.order0_window(bg, v, lam_f, n_gauss=n_gauss)
    H = H_exact(ds, bg, v, lam_f, n_gauss=n_gauss)
    tru = zeta_contraction(pa, math.cos(math.radians(gamma / 60.0)), v)
    return 2.0 * float(np.trapezoid(W * H * tru, v))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.5, 1.0, 2.0, 5.0, 17.3])
    ap.add_argument("--ranges", type=str, nargs="+",
                    default=["406.0:2313.029", "430.0:2280.0"])
    ap.add_argument("--n-out", type=int, default=400)
    a = ap.parse_args()

    ds, bgm, _ = wire()
    import perm_aware_kappa3_callable as pa
    lam_f = bgm.LAM_SOURCE_BASELINE
    bg = bgm.Background(lam_min=120.0, lam_max=lam_f + 5.0)
    g_fold, v_fold = load_fold_kk(FOLD)

    rngs = [tuple(float(t) for t in r.split(":")) for r in a.ranges]
    print(f"{'gamma':>7} " + "  ".join(f"{r:>16}" for r in a.ranges)
          + f"  {'fold':>13}")
    for g in a.gammas:
        vals = [fk_cont(ds, bg, pa, g, lam_f, lo, hi, n_out=a.n_out)
                for lo, hi in rngs]
        fold = float(np.interp(g, g_fold, v_fold))
        print(f"{g:7.2f} " + "  ".join(f"{v:16.6e}" for v in vals)
              + f"  {fold:13.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
