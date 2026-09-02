"""Exact expectation of the estimator at production settings, plus references."""
import sys, time
import numpy as np
import fk_expect_exact as ex

gam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
sig = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
N = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
calib = sys.argv[4] if len(sys.argv) > 4 else "smeared"
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
grid = ex.build_grid(n_lambda=N)
t = time.time(); V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
M = ex.calib_matrix(grid, V, A, rho, sig, calib)
zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
zi = ex.zeta_injected(M, Q)
ref_true = ex.fk_reference(grid, zt)[0, 3]
ref_inj = ex.fk_reference(grid, zi)[0, 3]
print(f"[setup {time.time()-t:.0f}s] gamma={gam}' sigma={sig} N={N} "
      f"dlam={grid.dlam:.3f} sig/dlam={sig/grid.dlam:.2f} rho={rho:.4f} calib={calib}")
t = time.time()
T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
tot = (T1 + T2) + (T1 + T2).T
vr = T1 + T1.T
print(f"[expectation {time.time()-t:.0f}s]")
print(f"  local-limit reference, true zeta      : {ref_true:.6e}")
print(f"  local-limit reference, injected zeta  : {ref_inj:.6e} "
      f"({ref_inj/ref_true:.4f} of true)")
print(f"  exact estimator expectation  T1+T2    : {tot[0,3]:.6e} "
      f"({tot[0,3]/ref_inj:.4f} of injected-zeta reference)")
print(f"  exact T1 only (fk_vr structure)       : {vr[0,3]:.6e} "
      f"({vr[0,3]/tot[0,3]:.4f} share)")
np.savez(f"_exact_g{gam}_s{sig}_N{N}.npz", T1=T1, T2=T2,
         ref_true=ref_true, ref_inj=ref_inj, gamma=gam, sigma=sig, N=N)
