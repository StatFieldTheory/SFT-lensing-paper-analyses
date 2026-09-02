import sys
import numpy as np
sys.path.insert(0, "..")
import fk_expect_exact as ex
from kappa3_norm import two_sided, local_prediction, smeared_B, kappa1_third_moment

N = 1000
grid = ex.build_grid(n_lambda=N)
gam = 1.0
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
Zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
w = grid.dlam * grid.Wd ** 3
print("Wd:", np.round(grid.Wd[[0,100,300,500,700,900,999]],2))
print("zeta_000 sign along ray:", np.round(Zt[[0,100,300,500,700,900,999],0,0,0],4))
print("sign changes in zeta_000:", int(np.sum(np.diff(np.sign(Zt[:,0,0,0]))!=0)))
L = float(np.sum(w*Zt[:,0,0,0]))
print("local prediction 000 =", L)
print("sum |w zeta| =", float(np.sum(np.abs(w*Zt[:,0,0,0]))))
for sig in (16.0, 8.0, 2.0):
    V,A,Q,rho = ex.node_stats(grid, cosg, sig, calibrate="nominal")
    S2 = V*(2*sig); B = smeared_B(grid, A, rho)
    inj = np.array([ex._DS.cum3_from_Q(Q[k], B[k])[0,0,0] for k in range(N)])
    injS = np.array([ex._DS.cum3_from_Q(Q[k], S2[k])[0,0,0] for k in range(N)])
    m3 = kappa1_third_moment(grid, A, Q, rho)[0,0,0]
    print(f"sig={sig}: exact m3={m3:.5e}  local-with-B={np.sum(w*inj):.5e} "
          f"local-with-S2={np.sum(w*injS):.5e}  L={L:.5e}")
    print("   inj/zeta at k=0,5,20,100,500,900,999:",
          np.round((inj/Zt[:,0,0,0])[[0,5,20,100,500,900,999]],4))
    print("   injS/zeta (should be 1):",
          np.round((injS/Zt[:,0,0,0])[[0,5,20,100,500,900,999]],6))
