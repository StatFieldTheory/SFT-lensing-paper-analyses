"""Is the smeared covariance a valid covariance (PD), and how does it differ?"""
import numpy as np, fk_expect_exact as ex
gam, N = 1.0, 1000
cosg = float(np.cos(np.deg2rad(gam/60.0)))
grid = ex.build_grid(n_lambda=N)
for sig in (2.0, 8.0, 32.0):
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate="smeared")
    B = ex.smeared_covariance(grid, A, rho)
    S2 = V*(2*sig)
    mn = np.array([np.linalg.eigvalsh(B[k]).min() for k in range(N)])
    # subspace alignment: angle between B's and Sigma2's smallest eigenvector
    ang = []
    for k in (100, 500, 900):
        _, Ub = np.linalg.eigh(B[k]); _, Us = np.linalg.eigh(S2[k])
        ang.append(abs(float(Ub[:, 0] @ Us[:, 0])))
    rel = np.abs(B-S2).max()/np.abs(S2).max()
    # relative deviation measured in the Sigma2 metric (what Q amplifies)
    amp = []
    for k in (100, 500, 900):
        d, U = np.linalg.eigh(S2[k])
        Bt = U.T @ B[k] @ U
        amp.append(float(np.abs(Bt[0, :]/np.sqrt(d[0]*d)).max()))
    print(f"sigma={sig:5.1f}  min eig(B) > 0: {bool((mn>0).all())}  "
          f"|B-S2|/|S2| = {rel:.3e}")
    print(f"    |<u_min(B)|u_min(S2)>| = {np.round(ang,4)}   "
          f"max |B_0j|/sqrt(d_0 d_j) = {np.round(amp,3)}")
