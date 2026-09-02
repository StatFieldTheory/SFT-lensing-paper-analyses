"""Is the nominal calibration's failure node-localised (isolated near-singular
V) rather than a smooth small-gamma effect?"""
import numpy as np, fk_expect_exact as ex
for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=1000)
    Vn, A, Qn, rho = ex.node_stats(grid, cosg, 8.0, calibrate="nominal")
    _, _, Qs, _ = ex.node_stats(grid, cosg, 8.0, calibrate="smeared")
    B = ex.smeared_covariance(grid, A, rho)
    S2 = Vn * 16.0
    condV = np.array([np.linalg.cond(S2[k]) for k in range(1000)])
    condB = np.array([np.linalg.cond(B[k]) for k in range(1000)])
    qn = np.abs(Qn).reshape(1000, -1).max(axis=1)
    qs = np.abs(Qs).reshape(1000, -1).max(axis=1)
    kern = grid.dlam * grid.Wd * grid.Hd
    print(f"gamma={gam}'  cond(Sigma2): median {np.median(condV):.2f}, "
          f"max {condV.max():.3e} at lam={grid.lam[condV.argmax()]:.1f}")
    print(f"            cond(B):      median {np.median(condB):.2f}, "
          f"max {condB.max():.3e}")
    print(f"            |Q|max nominal: median {np.median(qn):.3e}, "
          f"max {qn.max():.3e} at lam={grid.lam[qn.argmax()]:.1f}")
    print(f"            |Q|max smeared: median {np.median(qs):.3e}, "
          f"max {qs.max():.3e}")
    bad = np.argsort(-qn)[:5]
    print(f"            5 worst nominal nodes: lam={np.round(grid.lam[bad],1)}, "
          f"cond(S2)={np.array2string(condV[bad], precision=1)}, "
          f"kernel weight frac={np.round(kern[bad]/kern.sum(),5)}")
    print(f"            ratio |Q_nom|/|Q_sm| at those nodes: "
          f"{np.array2string(qn[bad]/qs[bad], precision=2)}")
