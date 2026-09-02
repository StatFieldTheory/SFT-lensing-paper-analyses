"""PART 3d -- excise ONLY the handful of near-singular Sigma2 nodes."""
import sys
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                                # noqa: E402
from kappa3_norm import smeared_B                           # noqa: E402

_DS = ex._DS
N, SIG = 1000, 8.0
grid = ex.build_grid(n_lambda=N)
print(f"N={N} sigma={SIG}.  'exact/inj' = exact FK divided by the local FK the")
print("SAME (masked) Q was designed to inject -- 1.0 means faithful injection.")
print(f"{'gamma':>7} {'excised nodes':>26} {'nom exact/inj':>14} "
      f"{'sm exact/inj':>13} {'nom/sm FK':>10} {'sm FK':>12} {'sm FK shift':>12}")
for gam in (1.0, 5.0, 17.3, 30.0):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    S2 = np.array([ex._CORE._assemble_6x6(grid.builder, cosg, float(l))
                   for l in grid.lam])
    Z = np.array([_DS.zeta6(cosg, float(l)) for l in grid.lam])
    V = S2 / (2 * SIG)
    rho = float(np.exp(-grid.dlam / SIG))
    A = np.empty_like(V); A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1 - rho * rho) * V[k]
    B = smeared_B(grid, A, rho)
    base_sm = None
    for label, mask in [("none", []),
                        ("472,682,683", [472, 682, 683]),
                        ("470-476,680-686", list(range(470, 477)) + list(range(680, 687)))]:
        Qn = np.array([_DS.solve_Q(S2[k], Z[k]) for k in range(N)])
        Qs = np.array([_DS.solve_Q(B[k], Z[k]) for k in range(N)])
        for k in mask:
            Qn[k] = 0.0; Qs[k] = 0.0
        rn = ex.fk_reference(grid, ex.zeta_injected(S2, Qn))[0, 3]
        rs = ex.fk_reference(grid, ex.zeta_injected(B, Qs))[0, 3]
        T1, T2 = ex.expectation(grid, cosg, SIG, V=V, A=A, Q=Qn, rho=rho)
        fn = ((T1 + T2) + (T1 + T2).T)[0, 3]
        T1, T2 = ex.expectation(grid, cosg, SIG, V=V, A=A, Q=Qs, rho=rho)
        fs = ((T1 + T2) + (T1 + T2).T)[0, 3]
        if base_sm is None:
            base_sm = fs
        print(f"{gam:>7.2f} {label:>26} {fn/rn:>14.4f} {fs/rs:>13.4f} "
              f"{fn/fs:>10.4f} {fs:>12.5e} {fs/base_sm:>12.4f}", flush=True)
