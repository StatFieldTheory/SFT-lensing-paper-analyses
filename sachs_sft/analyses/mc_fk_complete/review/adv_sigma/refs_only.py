"""Cheap: only the LOCAL references (no O(N^2) expectation).

Gives, per (N, sigma):
  ref_true  = discrete local fold of the TABULATED zeta
  ref_inj   = discrete local fold of the cumulant the deformation ACTUALLY injects
  fold      = the paper's analytic FK
so we can see where the sigma_lambda -> 0 limit of the estimator actually sits
as a function of the lattice, without any extrapolation at all.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import fk_expect_exact as ex

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--Ns", type=int, nargs="+", required=True)
    ap.add_argument("--sigmas", type=float, nargs="+", required=True)
    ap.add_argument("--calib", default="smeared")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    cosg = float(np.cos(np.deg2rad(a.gamma / 60.0)))
    rows = []
    for N in a.Ns:
        grid = ex.build_grid(n_lambda=N)
        zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
        ref_true = ex.fk_reference(grid, zt)[0, 3]
        for s in a.sigmas:
            t0 = time.time()
            V, A, Q, rho = ex.node_stats(grid, cosg, s, calibrate=a.calib)
            M = ex.calib_matrix(grid, V, A, rho, s, a.calib)
            zi = ex.zeta_injected(M, Q)
            ref_inj = ex.fk_reference(grid, zi)[0, 3]
            # injection fidelity in tensor norm, kernel weighted
            kern = grid.dlam * grid.Wd * grid.Hd
            zts = sum(np.transpose(zt, (0,) + tuple(1 + np.array(pp)))
                      for pp in ((0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0))) / 6.0
            num = np.einsum("m,mabc->", np.abs(kern), np.abs(zi - zts))
            den = np.einsum("m,mabc->", np.abs(kern), np.abs(zt))
            rows.append(dict(N=N, dlam=grid.dlam, sigma=s, ref_true=ref_true,
                             ref_inj=ref_inj, fid=ref_inj/ref_true, tnorm=num/den))
            print(f"N={N:>5d} dlam={grid.dlam:>7.4f} sig={s:>6.2f} r={s/grid.dlam:>6.2f} "
                  f"ref_true={ref_true:.8e} ref_inj={ref_inj:.8e} "
                  f"inj/true={ref_inj/ref_true:.6f} |dz|/|z|={num/den:.3e} [{time.time()-t0:.0f}s]",
                  flush=True)
    if a.out:
        np.savez(a.out, **{k: np.array([r[k] for r in rows]) for k in rows[0]})
        print("->", a.out)
