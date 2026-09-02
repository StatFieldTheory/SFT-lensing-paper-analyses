"""Where does B deviate from Sigma2, and does the kernel weight it?"""
import numpy as np, fk_expect_exact as ex
gam, N, sig = 1.0, 1000, 8.0
cosg = float(np.cos(np.deg2rad(gam/60.0)))
grid = ex.build_grid(n_lambda=N)
V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate="smeared")
B = ex.smeared_covariance(grid, A, rho); S2 = V*(2*sig)
r = np.array([np.abs(B[k]).max()/np.abs(S2[k]).max() for k in range(N)])
kern = grid.dlam*grid.Wd*grid.Hd; kern = kern/kern.sum()
print("B/Sigma2 (per-node max-norm ratio):")
for k in (0, 1, 3, 10, 50, 200, 500, 800, 950, 990, 998, 999):
    print(f"  k={k:4d} lam={grid.lam[k]:7.1f}  B/S2={r[k]:.4f}  "
          f"kernel weight={kern[k]:.3e}")
print(f"kernel-weighted mean B/S2 = {(r*kern).sum():.4f}")
mn = np.array([np.linalg.eigvalsh(B[k]).min() for k in range(N)])
bad = np.where(mn <= 0)[0]
print(f"nodes with a non-positive eigenvalue of B: {bad[:12]} (n={bad.size}); "
      f"their total kernel weight = {kern[bad].sum():.3e}")
