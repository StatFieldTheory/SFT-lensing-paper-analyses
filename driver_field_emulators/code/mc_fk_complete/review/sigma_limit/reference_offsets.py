"""Where the sigma-INDEPENDENT offsets between the estimator's local reference
and the paper's analytic FK live.  (None of this depends on sigma_lambda; it is
the part of the error budget the sigma->0 extrapolation does NOT address.)"""
from __future__ import annotations
import sys
import numpy as np

sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/"
                   "driver_field_emulators/code/mc_fk_complete")
import fk_expect_exact as ex   # noqa: E402

PAPER = {0.5: 1.948e-05, 1.0: 1.598e-05, 2.0: 1.223e-05,
         5.0: 7.038e-06, 17.3: 1.966e-06}
PAPER_FULL = {1.0: 1.597716e-05}

for gam in (1.0, 5.0):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    print("=" * 88)
    print(f"gamma = {gam}'   paper analytic FK = {PAPER[gam]:.6e}")
    print("=" * 88)
    print(f"  {'N':>6} {'dlam':>8} {'discrete fold, FULL grid':>26} "
          f"{'restricted to table support':>30}")
    prev = None
    for N in (1000, 1999, 3997, 8000):
        grid = ex.build_grid(n_lambda=N)
        zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
        kern = grid.dlam * grid.Wd * grid.Hd
        full = np.einsum("m,Abc,mbcB->AB", kern, ex.F6, zt, optimize=True)
        full = (full + full.T)[0, 3]
        m = (grid.lam >= 430.0) & (grid.lam <= 2280.0)
        res = np.einsum("m,Abc,mbcB->AB", kern * m, ex.F6, zt, optimize=True)
        res = (res + res.T)[0, 3]
        print(f"  {N:>6} {grid.dlam:>8.4f} {full:>26.6e} {res:>30.6e}")
        prev = (full, res, grid.dlam)
    # O(dlam) Richardson from the last two
    print(f"  ratio to paper (finest grid): full {prev[0]/PAPER[gam]:.4f}   "
          f"restricted {prev[1]/PAPER[gam]:.4f}")
    # kernel identity: does the discrete Wd match driver_stats.order0_window?
    grid = ex.build_grid(n_lambda=2000)
    W0 = ex._DS.order0_window(grid.bg, grid.lam, float(ex._BG.LAM_SOURCE_BASELINE),
                              n_gauss=64)
    d = np.abs(grid.Wd - W0) / np.abs(W0).max()
    print(f"  discrete Wd vs driver_stats.order0_window: max rel dev "
          f"{d.max():.2e}")
