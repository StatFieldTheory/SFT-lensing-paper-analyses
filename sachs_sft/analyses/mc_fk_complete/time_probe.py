import time, numpy as np, fk_expect_exact as ex, fk_complete_core as fc
gam, sig, N = 1.0, 8.0, 1000
cosg = float(np.cos(np.deg2rad(gam/60.0)))
t=time.time(); grid = ex.build_grid(n_lambda=N); print(f"grid {time.time()-t:.1f}s dlam={grid.dlam:.3f} sig/dlam={sig/grid.dlam:.2f}")
t=time.time(); nodes = ex.node_stats(grid, cosg, sig); print(f"node_stats {time.time()-t:.1f}s")
V,A,Q,rho = nodes
print("rho =", rho, " |Q|max =", np.abs(Q).max())
zi = ex.zeta_injected(V, Q, sig)
zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
ref_true = ex.fk_reference(grid, zt); ref_inj = ex.fk_reference(grid, zi)
print(f"discrete local-limit FK (true zeta)     kk = {ref_true[0,3]:.6e}")
print(f"discrete local-limit FK (injected zeta) kk = {ref_inj[0,3]:.6e}  "
      f"fidelity {ref_inj[0,3]/ref_true[0,3]:.5f}")
np.savez("_cache_nodes_g1_s8_N1000.npz", V=V, A=A, Q=Q, rho=rho,
         zeta_true=zt, zeta_inj=zi, lam=grid.lam)
t=time.time(); r = fc.simulate_fk_complete(gam, sigma_lambda=sig, n_lambda=N,
        n_real=2000, batch_size=2000, seed=7, grid=grid, nodes=nodes)
print(f"MC 1 batch of 2000: {time.time()-t:.1f}s  kk={r.kk:.4e}")
