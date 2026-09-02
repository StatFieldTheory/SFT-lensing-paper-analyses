"""Deterministic cross-check of the folded FK channel against a direct contraction.

The FK contribution to xi_kappa(gamma) is one diagram: the third cumulant of
the driving field, with two legs closed through the F vertex on one ray and
the remaining leg carried to the observable on the other ray.  Written out,
the kernel of that diagram factorises into the Order-0 window W(v) on the
F-free leg and the response integral

    H(v) = D(v)^4 * int_v^{lam_f} W(u) / D(u)^4 du

on the F-carrying ray (R(u,v) = [D(v)/D(u)]^2 is the Sachs response), so

    FK(gamma) = 2 * int dv  W(v) H(v) * [ -sum_c zeta6_{c c B}(gamma; v) ],

the factor 2 counting the two ray orderings and zeta6 the two-ray (6,6,6)
third-cumulant density built from the tabulated vertex.  This script
evaluates that integral directly -- plain quadrature against the tabulated
cumulant, no sft-wick machinery -- and compares it with the production fold.
The two computations share the vertex table by construction; everything else
(the kernel, the integration) is independent, so agreement validates the
workflow's assembly of the diagram.

Usage::

    <PyCCL python, PYTHONPATH=canoes/src> fk_kernel_crosscheck.py \
        --table <vertex table npz> --fold <folded xi npz> [--n-lam 48]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"
_FIX = _REPO / "driver_field_emulators" / "code" / "callable_fixed"


def load_fold_kk(path: Path):
    """gamma [arcmin] and FK value for the kappa-kappa order-2 entries."""
    d = np.load(path, allow_pickle=True)
    m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)), -1, 1)))) * 60)
    g = np.asarray(g)
    v = np.asarray(d["value"], float)[m]
    order = np.argsort(g)
    return g[order], v[order]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--fold", type=Path, required=True)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.5, 0.8, 1.3, 2.0, 3.2, 5.0, 8.0, 12.0, 17.3,
                             28.0, 44.4, 70.0])
    ap.add_argument("--n-lam", type=int, default=48)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    sys.path.insert(0, str(_MC))
    sys.path.insert(0, str(_FIX))
    import driver_stats as ds
    import perm_aware_kappa3_callable as pa
    import background as bgm
    pa.TABLE_PATH = a.table.resolve()
    pa._CACHE = None
    ds._k3 = pa

    lam_f = bgm.LAM_SOURCE_BASELINE
    bg = bgm.Background(lam_min=120.0, lam_max=lam_f + 5.0)

    # The vertex table's shells span ~[430, 2280]; integrate over that support.
    lam = np.linspace(430.0, 2280.0, a.n_lam)
    W = ds.order0_window(bg, lam, lam_f, n_gauss=64)
    D = np.asarray(bg.D(lam), float)

    lam_u = np.linspace(430.0, lam_f, 600)
    Wu = ds.order0_window(bg, lam_u, lam_f, n_gauss=64)
    Du = np.asarray(bg.D(lam_u), float)
    gu = Wu / Du**4
    H = np.array([D[i]**4 * np.trapezoid(gu[lam_u >= v], lam_u[lam_u >= v])
                  for i, v in enumerate(lam)])
    w_ker = W * H

    g_fold, v_fold = load_fold_kk(a.fold)

    print(f"[xcheck] table={a.table.name}, fold={a.fold.name}, "
          f"{a.n_lam} lambda nodes")
    print(f"{'gamma':>7} {'direct':>12} {'fold':>12} {'direct/fold':>11}")
    rows = []
    for g in a.gammas:
        cosg = math.cos(math.radians(g / 60.0))
        tru = np.array([-sum(ds.zeta6(cosg, float(v))[c, c, 3] for c in range(3))
                        for v in lam])
        direct = 2.0 * np.trapezoid(w_ker * tru, lam)
        fold = float(np.interp(g, g_fold, v_fold))
        rows.append((g, direct, fold, direct / fold))
        print(f"{g:7.2f} {direct:12.4e} {fold:12.4e} {direct/fold:11.4f}",
              flush=True)

    arr = np.array(rows)
    ratios = arr[:, 3]
    print(f"\n[xcheck] direct/fold: median {np.median(ratios):.3f}, "
          f"range [{ratios.min():.3f}, {ratios.max():.3f}]")
    if a.out:
        np.savez(a.out, gamma=arr[:, 0], direct=arr[:, 1], fold=arr[:, 2],
                 ratio=ratios, n_lam=a.n_lam)
        print(f"[xcheck] -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
