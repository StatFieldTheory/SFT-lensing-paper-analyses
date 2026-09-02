"""Two premises of the sigma->0 proof, checked numerically.

(i) The drive's injected three-point cumulant DENSITY is sigma-independent by
    construction under the 'smeared' calibration: cum3_from_Q(solve_Q(B,Z), B)
    = sym(Z) at every sigma.  So the sigma->0 limit cannot move the target.
(ii) B[k] = sum_l dlam <z_k z_l^T>  ->  Sigma2(lam_k) as sigma -> 0, at rate
    O(sigma).  This is what makes the 'smeared' calibration a reparametrisation
    of the same physics rather than a different model (solve_Q is homogeneous).
"""
from __future__ import annotations
import sys
import numpy as np

sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/"
                   "driver_field_emulators/code/mc_fk_complete")
import fk_expect_exact as ex   # noqa: E402

import os
gam, N = 1.0, int(os.environ.get("NL", 400))
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
grid = ex.build_grid(n_lambda=N)
S2 = np.array([ex._CORE._assemble_6x6(grid.builder, cosg, float(l))
               for l in grid.lam])
zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
interior = (grid.lam > 500) & (grid.lam < 2200)
print(f"gamma={gam}'  N={N}  dlam={grid.dlam:.3f}")
print(f"  {'sigma':>8} {'||B-Sigma2||/||Sigma2||':>24} {'/sigma':>10} "
       f"{'||zeta_inj-sym(zeta)||/||zeta||':>32}")
for sig in (32.0, 16.0, 8.0, 4.0, 2.0, 1.0):
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate="smeared")
    B = ex.calib_matrix(grid, V, A, rho, sig, "smeared")
    e = (np.abs(B - S2)[interior].max() / np.abs(S2).max())
    zi = ex.zeta_injected(B, Q)
    zs = (zt + zt.transpose(0, 2, 1, 3) + zt.transpose(0, 3, 2, 1)
          + zt.transpose(0, 1, 3, 2) + zt.transpose(0, 2, 3, 1)
          + zt.transpose(0, 3, 1, 2)) / 6.0
    ez = np.abs(zi - zs).max() / np.abs(zs).max()
    print(f"  {sig:>8.2f} {e:>24.6e} {e/sig:>10.3e} {ez:>32.3e}")
