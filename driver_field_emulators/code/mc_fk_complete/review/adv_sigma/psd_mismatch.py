"""Does the MC develop a sigma-dependent bias because it floors the
non-PSD tabulated covariance while the exact expectation does not?

fk_complete_core samples with ``_cholesky_psd(V[k])`` -- negative eigenvalues
clipped to zero -- whereas fk_expect_exact propagates ``V[k]`` as given through
the AR(1) recursion.  The two therefore describe DIFFERENT fields.  The AR(1)
smearing hides the defect at large sigma_lambda (neighbouring nodes average it
away) and exposes it as the regulator is removed, so the mismatch is a
candidate for a bias that grows exactly where the extrapolation is anchored.

Here the exact evaluator is run twice at each sigma: once with the field the
exact expectation assumes, once with the field the Monte-Carlo actually
simulates (same Q in both).  The difference is the whole systematic.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import fk_expect_exact as ex


def floor_psd(V):
    out = np.empty_like(V)
    nneg = 0
    worst = 0.0
    for k in range(V.shape[0]):
        w, U = np.linalg.eigh(V[k])
        if w.min() < 0.0:
            nneg += 1
            worst = max(worst, -w.min() / max(abs(w).max(), 1e-300))
            out[k] = (U * np.maximum(w, 0.0)) @ U.T
        else:
            out[k] = V[k]
    return out, nneg, worst


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--N", type=int, default=1000)
    ap.add_argument("--sigmas", type=float, nargs="+", default=[8.0, 4.0, 2.0, 1.0])
    a = ap.parse_args()
    cosg = float(np.cos(np.deg2rad(a.gamma / 60.0)))
    grid = ex.build_grid(n_lambda=a.N)
    print(f"gamma={a.gamma}' N={a.N} dlam={grid.dlam:.4f}")
    for s in a.sigmas:
        t0 = time.time()
        V, A, Q, rho = ex.node_stats(grid, cosg, s, calibrate="smeared")
        Vf, nV, wV = floor_psd(V)
        Af, nA, wA = floor_psd(A)
        # the field the MC simulates: AR(1) built from the FLOORED per-node V
        Af2 = np.empty_like(Vf)
        Af2[0] = Vf[0]
        for k in range(1, a.N):
            Af2[k] = rho * rho * Af2[k - 1] + (1.0 - rho * rho) * Vf[k]
        _, nA2, wA2 = floor_psd(Af2)
        T1, T2 = ex.expectation(grid, cosg, s, V=V, A=A, Q=Q, rho=rho, block=128)
        base = ((T1 + T2) + (T1 + T2).T)[0, 3]
        T1f, T2f = ex.expectation(grid, cosg, s, V=V, A=Af2, Q=Q, rho=rho, block=128)
        mcf = ((T1f + T2f) + (T1f + T2f).T)[0, 3]
        print(f"sig={s:>6.2f} r={s/grid.dlam:>6.2f} | V neg-nodes={nV:>3d} (worst {wV:.2e}) "
              f"A neg-nodes={nA:>3d} A_floored neg={nA2:>3d} | exact={base:.8e} "
              f"MC-field={mcf:.8e} ratio={mcf/base:.6f} ({(mcf/base-1)*100:+.3f}%) "
              f"[{time.time()-t0:.0f}s]", flush=True)
