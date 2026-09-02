"""Is the finite-sigma inflation the AR(1) smearing decohering A[k] from V[k]?"""
import numpy as np, fk_expect_exact as ex
gam, N = 1.0, 1000
cosg = float(np.cos(np.deg2rad(gam/60.0)))
grid = ex.build_grid(n_lambda=N)
for sig in (2.0, 8.0, 32.0):
    V, A, Q, rho = ex.node_stats(grid, cosg, sig)
    k = N//2
    dv = np.linalg.eigvalsh(V[k]); da = np.linalg.eigvalsh(A[k])
    rel = np.abs(A[k]-V[k]).max()/np.abs(V[k]).max()
    # projection of A onto V's small eigendirections
    w, U = np.linalg.eigh(V[k])
    proj = np.array([U[:, i] @ A[k] @ U[:, i] / w[i] for i in range(6)])
    print(f"sigma={sig:5.1f} rho={rho:.4f} |A-V|/|V|max={rel:.3e}")
    print(f"   eig V : {dv}")
    print(f"   eig A : {da}")
    print(f"   <u_i|A|u_i>/d_i (V eigbasis): {np.array2string(proj, precision=3)}")
    # how fast does Sigma2 vary over sigma?
    dV = np.abs(V[k+int(max(sig/grid.dlam,1))]-V[k]).max()/np.abs(V[k]).max()
    print(f"   |V(k+sigma)-V(k)|/|V| = {dV:.3e}")
