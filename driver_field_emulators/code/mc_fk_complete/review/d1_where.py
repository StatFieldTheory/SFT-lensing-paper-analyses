import sys
import numpy as np
sys.path.insert(0, "..")
import fk_expect_exact as ex
from kappa3_norm import two_sided, local_prediction, smeared_B

N = 1000
grid = ex.build_grid(n_lambda=N)
gam = 1.0
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
Zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
np.set_printoptions(precision=4, suppress=False, linewidth=200)
for sig in (16.0, 8.0, 2.0):
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate="nominal")
    S2 = V * (2 * sig)
    B = smeared_B(grid, A, rho)
    print(f"--- sigma={sig} rho={rho:.6f} dlam={grid.dlam:.4f} r={grid.dlam/sig:.4f}")
    for k in (0, 5, 20, 100, 500, 900, 999):
        print(f"  k={k:4d} B00/S200={B[k,0,0]/S2[k,0,0]:8.4f} "
              f"B03/S203={B[k,0,3]/S2[k,0,3]:8.4f} "
              f"A00/V00={A[k,0,0]/V[k,0,0]:8.4f}")
    # node-wise injected-vs-target ratio for channel 000 and 003
    ratio000 = np.array([ex._DS.cum3_from_Q(Q[k], B[k])[0,0,0]/Zt[k,0,0,0] for k in range(N)])
    print("  cum3(Q_nom,B)/zeta  [0,0,0] at k=0,5,20,100,500,900,999:",
          np.round(ratio000[[0,5,20,100,500,900,999]], 4))
    # integrand of the local prediction
    w = grid.dlam * grid.Wd**3
    integ = w * Zt[:, 0, 0, 0]
    cum = np.cumsum(integ)/np.sum(integ)
    print("  local-integrand cumulative fraction at k=5,20,100,500:",
          np.round(cum[[5,20,100,500]], 4), " Wd[0]/Wd[500] =",
          round(float(grid.Wd[0]/grid.Wd[500]),3))
    # exact vs "local with B" prediction
    Y = two_sided(A, rho, grid.dlam, grid.Wd)
    print("  |Y[k] - Wd[k] B[k]| / |Wd B| at k=100,500,900:",
          [round(float(np.abs(Y[k]-grid.Wd[k]*B[k]).max()/np.abs(grid.Wd[k]*B[k]).max()),4)
           for k in (100,500,900)])
