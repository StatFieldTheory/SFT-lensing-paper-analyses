"""EXACT (not Monte-Carlo) check of which node covariance the AR(1) chain has.

z[0] = L0 xi0 ;  z[k] = rho z[k-1] + c L[k] xi_k,  c = sqrt(1-rho^2).
Build the explicit linear map z = G xi and form Cov = G G^T exactly, then
compare against rho^|k-l| A[min] (A the exact AR(1) variance recursion) and
against rho^|k-l| V[min] (the nominal per-node V).  A 6-line, assumption-free
adjudication of `smeared` vs `smeared_V`, with no FK anywhere in sight.

Also reports, for the SAME grid, the object the calibration needs:
B[k] = sum_l dlam <z[k] z[l]^T>.
"""
import sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/review/calib_attack")
import numpy as np, altcal
import fk_expect_exact as ex

for (N, sig) in [(200, 8.0), (400, 8.0), (400, 32.0)]:
    gam = 1.0152
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=N)
    V, A, Z, rho = altcal.base_stats(grid, cosg, sig)
    L = np.empty_like(V)
    for k in range(N):
        d, U = np.linalg.eigh(V[k])
        L[k] = U @ np.diag(np.sqrt(np.clip(d, 0.0, None))) @ U.T
    c = np.sqrt(1.0 - rho * rho)
    G = np.zeros((N, 6, N, 6))                    # z[k] = sum_j G[k,:,j,:] xi_j
    G[0, :, 0, :] = L[0]
    for k in range(1, N):
        G[k] = rho * G[k - 1]
        G[k, :, k, :] += c * L[k]
    Gf = G.reshape(N * 6, N * 6)
    Cov = (Gf @ Gf.T).reshape(N, 6, N, 6).transpose(0, 2, 1, 3)   # (N,N,6,6)
    idx = np.arange(N)
    lag = (rho ** np.abs(idx[:, None] - idx[None, :]))[:, :, None, None]
    predA = lag * A[np.minimum(idx[:, None], idx[None, :])]
    predV = lag * V[np.minimum(idx[:, None], idx[None, :])]
    sc = np.abs(Cov).max()
    print(f"N={N} sigma={sig} dlam={grid.dlam:.3f} rho={rho:.5f}")
    print(f"   max |Cov - rho^|k-l| A[min]| / max|Cov| = "
          f"{np.abs(Cov-predA).max()/sc:.3e}    <-- 'smeared' uses A")
    print(f"   max |Cov - rho^|k-l| V[min]| / max|Cov| = "
          f"{np.abs(Cov-predV).max()/sc:.3e}    <-- 'smeared_V' uses V")
    Btrue = grid.dlam * Cov.sum(axis=1)
    Bsm = altcal.calib_matrices(grid, V, A, rho, sig, "smeared")
    BsmV = altcal.calib_matrices(grid, V, A, rho, sig, "smeared_V")
    kern = grid.dlam * grid.Wd * grid.Hd; kern = kern / kern.sum()
    e1 = np.array([np.abs(Bsm[k]-Btrue[k]).max()/np.abs(Btrue[k]).max() for k in range(N)])
    e2 = np.array([np.abs(BsmV[k]-Btrue[k]).max()/np.abs(Btrue[k]).max() for k in range(N)])
    print(f"   B: kernel-weighted mean |smeared   - exact|/|exact| = {float((e1*kern).sum()):.3e}")
    print(f"   B: kernel-weighted mean |smeared_V - exact|/|exact| = {float((e2*kern).sum()):.3e}")
