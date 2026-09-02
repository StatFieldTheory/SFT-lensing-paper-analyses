"""The lambda-kernel probe the gamma sweep cannot give.

The gamma sweep moves the FK integrand's median lambda only from 1868 to 1706
Mpc (c1), so all nine separations weight the shared table almost identically:
they are close to ONE constraint on the kernel.  Changing the SOURCE distance
does move it: the kernel W(v)H(v) is rebuilt on [406, t_f] for each t_f, while
the shared table content is untouched.  If the estimator tracks the fold across
t_f = 1331 -> 2318 Mpc, the agreement is a statement about the kernel's
lambda profile, which is NOT shared.

Fold: products/table_permclosed_cut1000_permfixmultiz_xi.npz  (7 source
distances, same table + same perm-aware callable as the estimator here).
"""
import _env, time, numpy as np
import fk_expect_exact as ex, inj_exact as ij
import perm_aware_kappa3_callable as pa
from _fold import fold_nodes
from pathlib import Path

PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/products")
pa.TABLE_PATH = PROD / "table_permclosed_cut1000.npz"
pa._CACHE = None
ex._DS._k3 = pa

FOLD = PROD / "table_permclosed_cut1000_permfixmultiz_xi.npz"
d = np.load(FOLD, allow_pickle=True)
TFS = np.unique(d["t_final"])
N = 1000
SIGMAS = (8.0, 4.0)

g0, _ = fold_nodes(FOLD, t_final=TFS[-1])
# pick fold gamma nodes closest to 1, 5, 17.3 arcmin
targets = [1.0, 5.0, 17.3]
gsel = [float(g0[np.argmin(np.abs(g0 - t))]) for t in targets]
print("fold gamma nodes used:", [f"{g:.4f}'" for g in gsel])
print(f"{'t_final':>10} {'gamma':>9} {'fold':>13} {'exact s=8':>12} "
      f"{'s=4':>9} {'r(s->0)':>9}")
rows = []
for tf in TFS:
    gg, vv = fold_nodes(FOLD, t_final=tf)
    grid = ex.build_grid(n_lambda=N, lam_source=float(tf))
    for gs in gsel:
        i = int(np.argmin(np.abs(gg - gs)))
        fold = float(vv[i])
        cosg = float(np.cos(np.deg2rad(gg[i] / 60.0)))
        rr = {}
        for sig in SIGMAS:
            V, A, Q, rho = ex.node_stats(grid, cosg, sig)
            T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
            rr[sig] = float(((T1 + T2) + (T1 + T2).T)[0, 3]) / fold
        r0 = 2 * rr[4.0] - rr[8.0]
        rows.append((tf, gg[i], fold, rr[8.0], rr[4.0], r0))
        print(f"{tf:>10.1f} {gg[i]:>8.3f}' {fold:>13.4e} {rr[8.0]:>12.4f} "
              f"{rr[4.0]:>9.4f} {r0:>9.4f}", flush=True)
a = np.array([r[-1] for r in rows])
print(f"\nsigma->0 ratio across 7 source distances x 3 separations: "
      f"median {np.median(a):.4f}  min {a.min():.4f}  max {a.max():.4f}")
print(f"amplitude range of the fold over t_final: "
      f"{max(r[2] for r in rows)/min(r[2] for r in rows):.1f}x")
np.savez("_c4_multisource.npz", rows=np.array(rows))
