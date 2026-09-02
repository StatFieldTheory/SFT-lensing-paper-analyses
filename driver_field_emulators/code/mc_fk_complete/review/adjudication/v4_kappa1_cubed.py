"""ADJUDICATOR check 4: the F-free calibration discriminator.

<kappa1^3> at O(Q) involves NO F vertex, NO FK kernel and NO comparison with
the paper's fold, so it cannot have been tuned to the FK answer.  It must equal
sum_j dlam Wd_j^3 zeta_sym[j].
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/review")
import fk_expect_exact as ex
import kappa3_norm as kn

N = 1000
for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam/60.0)))
    grid = ex.build_grid(n_lambda=N)
    zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    zs = (zt + zt.transpose(0,1,3,2) + zt.transpose(0,2,1,3) + zt.transpose(0,2,3,1)
          + zt.transpose(0,3,1,2) + zt.transpose(0,3,2,1))/6.0
    loc = kn.local_prediction(grid, zs)
    print(f"\ngamma={gam}'  N={N}   [channels (0,0,0) (0,0,3) (0,3,3) (0,1,1)]")
    for calib in ("nominal", "smeared"):
        for sig in (8.0, 4.0, 2.0):
            V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
            m3 = kn.kappa1_third_moment(grid, A, Q, rho)
            r = [m3[i]/loc[i] for i in ((0,0,0),(0,0,3),(0,3,3),(0,1,1))]
            print(f"  {calib:<8} sigma={sig:>5.1f}  " +
                  "  ".join(f"{x:>9.4f}" for x in r), flush=True)
