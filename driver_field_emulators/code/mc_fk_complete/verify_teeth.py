"""How large an error in the FK assembly would this comparison detect?

Perturb the discrete reference kernel in physically motivated ways and see how
far the ratio to the paper's fold moves. The measured agreement is
1.002 +/- 0.006 at gamma = 1', so any perturbation moving the ratio by more than
about 1.5% would have been caught.
"""
import numpy as np, fk_expect_exact as ex
from run_gate2 import analytic_fk

N = 1000
grid = ex.build_grid(n_lambda=N)
D, Wd, dlam = grid.D, grid.Wd, grid.dlam
print(f"{'perturbation':<52} {'gamma=1':>9} {'gamma=5':>9} {'gamma=17.3':>11}")


def ref_with(Hd, zeta, scale=1.0):
    kern = dlam * Wd * Hd * scale
    M = np.einsum("m,Abc,mbcB->AB", kern, ex.F6, zeta, optimize=True)
    return (M + M.T)[0, 3]


rows = {}
for gam in (1.0, 5.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    zs = (zt + zt.transpose(0, 1, 3, 2) + zt.transpose(0, 2, 1, 3)
          + zt.transpose(0, 2, 3, 1) + zt.transpose(0, 3, 1, 2)
          + zt.transpose(0, 3, 2, 1)) / 6.0
    ana = float(analytic_fk([gam])[0])
    base = ref_with(grid.Hd, zs) / ana
    # 1. no causal restriction on the F node
    g_all = dlam * Wd[1:] / D[:-1] ** 4
    Hd_all = D ** 4 * np.full(N, g_all.sum())
    # 2. one propagator instead of two on the F leg (D^2 not D^4)
    g2 = dlam * Wd[1:] / D[:-1] ** 2
    t2 = np.zeros(N); t2[:-1] = np.cumsum(g2[::-1])[::-1]
    Hd_p2 = D ** 2 * t2
    # 3. trapezoid end weights dropped
    w = np.full(N, dlam)
    tail = np.cumsum((w / D ** 2)[::-1])[::-1]
    Wd_flat = D ** 2 * tail
    gflat = dlam * Wd_flat[1:] / D[:-1] ** 4
    tf = np.zeros(N); tf[:-1] = np.cumsum(gflat[::-1])[::-1]
    kern_flat = dlam * Wd_flat * D ** 4 * tf
    Mf = np.einsum("m,Abc,mbcB->AB", kern_flat, ex.F6, zs, optimize=True)
    rows.setdefault("baseline (symmetrised reference)", {})[gam] = base
    rows.setdefault("  drop causality: F node before the drives", {})[gam] = \
        ref_with(Hd_all, zs) / ana
    rows.setdefault("  D^2 instead of D^4 on the F leg", {})[gam] = \
        ref_with(Hd_p2, zs) / ana
    rows.setdefault("  drop trapezoid end weights", {})[gam] = (Mf + Mf.T)[0, 3] / ana
    rows.setdefault("  omit the A<->B symmetrisation (factor 2)", {})[gam] = base / 2
    rows.setdefault("  scale zeta by 1.05", {})[gam] = base * 1.05
    rows.setdefault("  use raw (unsymmetrised) zeta6 legs", {})[gam] = \
        ref_with(grid.Hd, zt) / ana
for k, v in rows.items():
    print(f"{k:<52} {v[1.0]:>9.4f} {v[5.0]:>9.4f} {v[17.3]:>11.4f}")
print("\nmeasured MC agreement at gamma=1': 1.0022 +/- 0.0056  ->  the comparison")
print("resolves a fractional error in the assembly of about 1.5% (2 sigma).")
