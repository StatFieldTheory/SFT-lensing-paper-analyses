"""Does Sigma2 (or zeta) have lambda-structure on the sigma_lambda scale at the
nodes where the nominal Q blows up?"""
import numpy as np, fk_expect_exact as ex
gam = 17.3
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
grid = ex.build_grid(n_lambda=1000)
V, A, Q, rho = ex.node_stats(grid, cosg, 8.0, calibrate="nominal")
S2 = V * 16.0
d = np.array([np.linalg.eigvalsh(S2[k]) for k in range(1000)])
Z = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
zc = np.array([-sum(Z[k][c, c, 3] for c in range(3)) for k in range(1000)])
print("lambda window around the worst node (1307 Mpc); d0 = smallest eigenvalue")
i0 = int(np.argmin(np.abs(grid.lam - 1307.0)))
print(f"{'lam':>8} {'d0':>12} {'d0/d0_smooth':>13} {'zeta_c':>13} {'|Q|max':>11}")
sm = np.convolve(d[:, 0], np.ones(21) / 21, mode="same")
for k in range(i0 - 6, i0 + 7):
    qm = np.abs(Q[k]).max()
    print(f"{grid.lam[k]:>8.1f} {d[k,0]:>12.4e} {d[k,0]/sm[k]:>13.4f} "
          f"{zc[k]:>13.4e} {qm:>11.3e}")
# how many nodes are anomalous, and what do they carry
ratio = d[:, 0] / np.convolve(d[:, 0], np.ones(21) / 21, mode="same")
kern = grid.dlam * grid.Wd * grid.Hd
bad = np.where(ratio < 0.5)[0]
print(f"\nnodes with d0 below half the local 21-node mean: {bad.size} of 1000, "
      f"lam = {np.round(grid.lam[bad][:10],1)}")
print(f"they carry {kern[bad].sum()/kern.sum():.4f} of the FK kernel weight")
