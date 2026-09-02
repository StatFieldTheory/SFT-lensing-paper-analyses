import sys
import numpy as np
sys.path.insert(0, "..")
import fk_expect_exact as ex

N = 1000
grid = ex.build_grid(n_lambda=N)
for gam in (1.0, 30.0):
    cosg = float(np.cos(np.deg2rad(gam/60.0)))
    S2 = np.array([ex._CORE._assemble_6x6(grid.builder, cosg, float(l))
                   for l in grid.lam])
    ev = np.array([np.linalg.eigvalsh(S2[k]) for k in range(N)])
    bad = np.where(ev[:,0] <= 0)[0]
    print(f"gamma={gam}': non-PD nodes = {bad.tolist()}  lam={np.round(grid.lam[bad],2).tolist()}")
    for k in bad:
        print(f"   k={k} eigs={['%.3e'%x for x in ev[k]]}  d_min/d_max={ev[k,0]/ev[k,-1]:.3e}")
    # how close to zero is the smallest eigenvalue elsewhere
    r = ev[:,0]/ev[:,-1]
    order = np.argsort(r)[:8]
    print("   8 smallest d_min/d_max:", ["%.3e"%x for x in r[order]], "at k=",order.tolist())
    # Sigma2_00 along the ray near the bad nodes
    print("   Sigma2_00 at k=470..475:", ["%.4e"%x for x in S2[470:476,0,0]])
    print("   Sigma2_00 at k=680..685:", ["%.4e"%x for x in S2[680:686,0,0]])
