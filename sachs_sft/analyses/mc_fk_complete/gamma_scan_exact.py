"""Exact T1/T2 split versus gamma, at fixed sigma_lambda."""
import sys
import numpy as np
import fk_expect_exact as ex

sig = float(sys.argv[1]) if len(sys.argv) > 1 else 8.0
N = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
gams = [float(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 \
    else [0.5, 1.0, 2.0, 5.0, 17.3]
calib = sys.argv[4] if len(sys.argv) > 4 else "smeared"
grid = ex.build_grid(n_lambda=N)
print(f"sigma={sig} N={N} calib={calib}")
print(f"{'gamma':>7} {'cond(V)':>10} {'T1':>12} {'T2':>12} {'T1+T2':>12} "
      f"{'share T1':>9} {'ratio_ref':>10}")
for gam in gams:
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
    M = ex.calib_matrix(grid, V, A, rho, sig, calib)
    ref = ex.fk_reference(grid, ex.zeta_injected(M, Q))[0, 3]
    d = np.linalg.eigvalsh(V[N // 2])
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    t1 = (T1 + T1.T)[0, 3]; tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
    print(f"{gam:>7.2f} {d.max()/d.min():>10.2e} {t1:>12.5e} "
          f"{tot-t1:>12.5e} {tot:>12.5e} {t1/tot:>9.4f} {tot/ref:>10.4f}",
          flush=True)
