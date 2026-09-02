import sys
import numpy as np
sys.path.insert(0, "..")
import fk_expect_exact as ex
from kappa3_norm import smeared_B

N = 1000
grid = ex.build_grid(n_lambda=N)
for gam in (1.0, 30.0):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    Zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    w = grid.dlam * grid.Wd ** 3
    print(f"===== gamma={gam}'")
    for sig in (16.0, 8.0, 2.0):
        V,A,Q,rho = ex.node_stats(grid, cosg, sig, calibrate="nominal")
        S2 = V*(2*sig); B = smeared_B(grid, A, rho)
        inj = np.array([ex._DS.cum3_from_Q(Q[k], B[k])[0,0,0] for k in range(N)])
        rat = inj/Zt[:,0,0,0]
        contrib = w*inj
        tot = contrib.sum()
        j = np.argsort(-np.abs(contrib))[:6]
        print(f" sig={sig}: total={tot:.4e}  ratio to L={tot/np.sum(w*Zt[:,0,0,0]):.3f}")
        print(f"   inj/zeta: min={rat.min():.3f} max={rat.max():.3f} "
              f"argmax={int(np.argmax(rat))} median={np.median(rat):.3f}")
        print(f"   top-6 |contribution| nodes: {j.tolist()}")
        print(f"     lam={np.round(grid.lam[j],1).tolist()}")
        print(f"     inj/zeta there={np.round(rat[j],2).tolist()}")
        print(f"     share of total={np.round(contrib[j]/tot,3).tolist()}")
        # how much of the total comes from the last 30 nodes / first 30 nodes
        print(f"   first30 share={contrib[:30].sum()/tot:.3f}  "
              f"last30 share={contrib[-30:].sum()/tot:.3f}  "
              f"last3 share={contrib[-3:].sum()/tot:.3f}")
