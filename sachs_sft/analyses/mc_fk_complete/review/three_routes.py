"""Side-by-side of the three routes to the local-limit FK reference.

A  clean-room discrete kernel of fk_ref_discrete.py (SPEC.md sec 1 only),
   N-sweep + Richardson extrapolation in 1/N;
B  its continuum limit, evaluated by accurate quadrature over the full
   grid span [406, 2313.029] with H computed from v itself;
C  ../rebuild/fk_kernel_crosscheck.py, run separately (values pasted from
   its own stdout at --n-lam 800, where it is converged to 5 digits).

Also reports the permutation-symmetrised variant of A/B: the single ordering
(n1,n1,n2) that the reference queries is replaced by the mean over the three
placements of the ray-2 (observable) leg.  The tabulated vertex is not exactly
closed under that exchange, and the spread is the dominant residual.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bootstrap import wire, FOLD  # noqa: E402
import fk_ref_discrete as R  # noqa: E402
from fk_ref_continuum import H_exact  # noqa: E402

GAMMAS = [0.5, 1.0, 2.0, 5.0, 17.3]
NS = [250, 500, 1000, 2000, 4000, 8000]
# route C, from fk_kernel_crosscheck.py --n-lam 800 (identical at --n-lam 200)
XCHECK = {0.5: 1.9272e-05, 1.0: 1.5912e-05, 2.0: 1.2257e-05,
          5.0: 7.1319e-06, 17.3: 2.1003e-06}


def contraction(pa, cos_gamma, lam, order):
    n1 = np.array([0.0, 0.0, 1.0])
    c = float(np.clip(cos_gamma, -1.0, 1.0))
    n2 = np.array([math.sqrt(max(0.0, 1.0 - c * c)), 0.0, c])
    rays = {0: n1, 1: n2}
    n = len(lam)
    n_arr = np.stack([np.repeat(rays[o][None, :], n, 0) for o in order], 0)
    z = pa.coupling_fn_batch(n_arr, np.stack([lam, lam, lam], 0))
    j = order.index(1)
    ii = [k for k in range(3) if k != j]
    tot = np.zeros(n)
    for cc in range(3):
        idx = [0, 0, 0]
        idx[ii[0]] = cc; idx[ii[1]] = cc; idx[j] = 0
        tot += z[:, idx[0], idx[1], idx[2]]
    return -tot


def main() -> int:
    ds, bgm, _ = wire()
    import perm_aware_kappa3_callable as pa
    lam_f = bgm.LAM_SOURCE_BASELINE
    bg = bgm.Background(lam_min=120.0, lam_max=lam_f + 5.0)
    g_fold, v_fold = R.load_fold_kk(FOLD)

    v = np.linspace(R.LAM_MIN, lam_f, 400)
    W = ds.order0_window(bg, v, lam_f, n_gauss=64)
    H = H_exact(ds, bg, v, lam_f)

    print("=" * 108)
    print("A  discrete N-sweep (clean-room kernel)")
    print("=" * 108)
    print(f"{'gamma':>6} " + " ".join(f"{'N=%d' % N:>13}" for N in NS)
          + f" {'A(Rich)':>13}")
    A = {}
    for g in GAMMAS:
        vals = [R.fk_ref(bg, pa, g, N, lam_f) for N in NS]
        A[g] = (vals, 2 * vals[-1] - vals[-2])
        print(f"{g:6.2f} " + " ".join(f"{x:13.6e}" for x in vals)
              + f" {A[g][1]:13.6e}")

    print()
    print("=" * 108)
    print("Three routes vs the paper fold (fold linearly interpolated onto "
          "the requested gamma, as NOTES.md and the crosscheck do)")
    print("=" * 108)
    print(f"{'gamma':>6} {'A discrete':>13} {'B continuum':>13} "
          f"{'C xcheck':>13} {'fold':>13} | {'A/fold':>8} {'B/fold':>8} "
          f"{'C/fold':>8} {'C/A':>8}")
    rows = []
    for g in GAMMAS:
        cg = math.cos(math.radians(g / 60.0))
        b = 2.0 * float(np.trapezoid(W * H * contraction(pa, cg, v, (0, 0, 1)), v))
        a = A[g][1]
        c = XCHECK[g]
        fo = float(np.interp(g, g_fold, v_fold))
        rows.append((g, a, b, c, fo))
        print(f"{g:6.2f} {a:13.6e} {b:13.6e} {c:13.6e} {fo:13.6e} | "
              f"{a/fo:8.4f} {b/fo:8.4f} {c/fo:8.4f} {c/a:8.4f}")

    print()
    print("=" * 108)
    print("At the fold's own gamma nodes (no interpolation), single ordering "
          "vs permutation-symmetrised")
    print("=" * 108)
    sel = [0, 3, 6, 10, 15]
    print(f"{'gamma_node':>10} {'fold':>13} {'B (n1,n1,n2)':>13} "
          f"{'B sym-mean':>13} | {'single/fold':>11} {'sym/fold':>9} "
          f"{'leg spread':>10}")
    for i in sel:
        g = float(g_fold[i]); fo = float(v_fold[i])
        cg = math.cos(math.radians(g / 60.0))
        vs = [2.0 * float(np.trapezoid(W * H * contraction(pa, cg, v, o), v))
              for o in ((0, 0, 1), (0, 1, 0), (1, 0, 0))]
        m = float(np.mean(vs))
        print(f"{g:10.5f} {fo:13.6e} {vs[0]:13.6e} {m:13.6e} | "
              f"{vs[0]/fo:11.5f} {m/fo:9.5f} {(max(vs)-min(vs))/m:10.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
