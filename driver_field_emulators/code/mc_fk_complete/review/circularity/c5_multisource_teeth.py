"""Does the multi-source test ADD constraints, or repeat one?

c4 showed the estimator tracks the fold across 7 source distances (fold
amplitude 196x).  That is only informative if a kernel defect would produce a
DIFFERENT ratio at different t_final.  Here the same 7 source distances are run
with deliberate kernel defects: if the defect's ratio varies with t_final, the
seven points are seven constraints on the kernel's lambda profile, not one.
"""
import _env, numpy as np
from pathlib import Path
import fk_expect_exact as ex, inj_exact as ij
import perm_aware_kappa3_callable as pa
from _fold import fold_nodes

PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/products")
pa.TABLE_PATH = PROD / "table_permclosed_cut1000.npz"; pa._CACHE = None
ex._DS._k3 = pa
FOLD = PROD / "table_permclosed_cut1000_permfixmultiz_xi.npz"
TFS = np.unique(np.load(FOLD, allow_pickle=True)["t_final"])
N = 1000
DEF = ["none", "acausal_both", "U_pow1", "F_no_leg_sym", "no_end_weights",
       "resp_pow1"]
gtar = 1.0
print(f"{'t_final':>9} {'lam quartiles of FK integrand':>32}   " +
      "  ".join(f"{d:>14}" for d in DEF))
rows = {}
for tf in TFS:
    gg, vv = fold_nodes(FOLD, t_final=tf)
    j = int(np.argmin(np.abs(gg - gtar)))
    fold = float(vv[j]); cosg = float(np.cos(np.deg2rad(gg[j] / 60.0)))
    grid = ex.build_grid(n_lambda=N, lam_source=float(tf))
    gp1 = ij.build_grid(n_lambda=N, lam_source=float(tf), resp_pow=1.0)
    gnw = ij.build_grid(n_lambda=N, lam_source=float(tf), end_weights=False)
    zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    kern = grid.dlam * grid.Wd * grid.Hd
    Mw = np.einsum("m,Abc,mbcB->mAB", kern, ex.F6, zt, optimize=True)
    w = 2.0 * Mw[:, 0, 3]; w = w / w.sum(); cw = np.cumsum(w)
    qs = [grid.lam[np.searchsorted(cw, f)] for f in (0.25, 0.5, 0.75)]
    r = {}
    for sig in (8.0, 4.0):
        V, A, Q, rho = ex.node_stats(grid, cosg, sig)
        for d in DEF:
            if d == "resp_pow1":
                T1, T2 = ij.expectation(gp1, cosg, sig, V, A, Q, rho)
            elif d == "no_end_weights":
                T1, T2 = ij.expectation(gnw, cosg, sig, V, A, Q, rho)
            else:
                T1, T2 = ij.expectation(grid, cosg, sig, V, A, Q, rho, defect=d)
            r[(sig, d)] = ij.fk_kk(T1, T2, d) / fold
    line = f"{tf:>9.1f} {qs[0]:>10.0f}/{qs[1]:.0f}/{qs[2]:.0f} Mpc        "
    for d in DEF:
        v = 2 * r[(4.0, d)] - r[(8.0, d)]
        rows.setdefault(d, []).append(v)
        line += f"{v:>16.4f}"
    print(line, flush=True)
print("\nspread of each defect's ratio across the 7 source distances:")
for d in DEF:
    a = np.array(rows[d])
    print(f"  {d:<16} min {a.min():.4f}  max {a.max():.4f}  "
          f"peak-to-peak {a.max()-a.min():.4f}")
