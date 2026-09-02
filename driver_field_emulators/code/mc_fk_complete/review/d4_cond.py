import sys
import numpy as np
sys.path.insert(0, "..")
import fk_expect_exact as ex

N = 1000
grid = ex.build_grid(n_lambda=N)
for gam in (1.0, 17.3, 30.0):
    cosg = float(np.cos(np.deg2rad(gam/60.0)))
    S2 = np.array([ex._CORE._assemble_6x6(grid.builder, cosg, float(l))
                   for l in grid.lam])
    ev = np.array([np.linalg.eigvalsh(S2[k]) for k in range(N)])
    cond = ev[:, -1]/ev[:, 0]
    neg = int(np.sum(ev[:, 0] <= 0))
    print(f"gamma={gam}': cond(Sigma2) median={np.median(cond):.3e} "
          f"max={cond.max():.3e} at k={int(np.argmax(cond))} (lam={grid.lam[np.argmax(cond)]:.1f}); "
          f"#nodes with a non-positive eigenvalue = {neg}")
    j = np.argsort(-cond)[:8]
    print("   worst nodes k:", j.tolist())
    print("   lam        :", np.round(grid.lam[j],1).tolist())
    print("   cond       :", ["%.2e"%c for c in cond[j]])
    print("   d_min/d_max:", ["%.2e"%x for x in (ev[j,0]/ev[j,-1])])
    print("   #nodes cond>1e3:", int(np.sum(cond>1e3)), " >1e4:", int(np.sum(cond>1e4)),
          " >1e6:", int(np.sum(cond>1e6)))
