"""PART 3c -- which nodes carry the amplification, and why."""
import sys
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                              # noqa: E402
from kappa3_norm import smeared_B                         # noqa: E402

_DS = ex._DS
N, SIG = 1000, 8.0
grid = ex.build_grid(n_lambda=N)
for gam in (1.0, 30.0):
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
    nrm = np.abs(S2).max(axis=(1, 2))
    med = float(np.median(nrm))
    Q = np.array([_DS.solve_Q(S2[k], Z[k]) for k in range(N)])
    epsn = np.zeros(N); err = np.zeros(N); cond = np.zeros(N)
    for k in range(N):
        d, U = np.linalg.eigh(S2[k])
        cond[k] = d.max() / d.min() if d.min() > 0 else np.inf
        if d.min() <= 0:
            epsn[k] = np.inf
        else:
            Wi = U @ np.diag(d ** -0.5) @ U.T
            epsn[k] = np.abs(np.linalg.eigvalsh(Wi @ (B[k] - S2[k]) @ Wi)).max()
        iM = _DS.cum3_from_Q(Q[k], S2[k]); iB = _DS.cum3_from_Q(Q[k], B[k])
        err[k] = np.abs(iB - iM).max() / max(np.abs(iM).max(), 1e-300)
    # contribution of each node to the FK kernel integral
    kern = grid.dlam * grid.Wd * grid.Hd
    ctrM = np.einsum("m,Abc,mbcB->mAB", kern, ex.F6,
                     np.array([_DS.cum3_from_Q(Q[k], S2[k]) for k in range(N)]),
                     optimize=True)[:, 0, 3]
    ctrB = np.einsum("m,Abc,mbcB->mAB", kern, ex.F6,
                     np.array([_DS.cum3_from_Q(Q[k], B[k]) for k in range(N)]),
                     optimize=True)[:, 0, 3]
    exc = ctrB - ctrM
    order = np.argsort(-np.abs(exc))
    print("=" * 104)
    print(f"gamma = {gam}'  sigma = {SIG}  N = {N}")
    print(f"  ||Sigma2|| percentiles vs median: p1={np.percentile(nrm,1)/med:.4f} "
          f"p5={np.percentile(nrm,5)/med:.4f} p50=1  p95={np.percentile(nrm,95)/med:.3f}")
    print(f"  ||eps||: median={np.median(epsn[np.isfinite(epsn)]):.4f} "
          f"p90={np.percentile(epsn[np.isfinite(epsn)],90):.4f} "
          f"p99={np.percentile(epsn[np.isfinite(epsn)],99):.4f} "
          f"max={epsn[np.isfinite(epsn)].max():.4g}   (#non-PD nodes: {int(np.sum(~np.isfinite(epsn)))})")
    tot = float(np.sum(ctrB)); base = float(np.sum(ctrM))
    print(f"  local FK from Q_nom with Sigma2 = {base:.5e}   with B = {tot:.5e}"
          f"   ratio {tot/base:.4f}")
    print(f"  {'rank':>4} {'k':>5} {'lam':>9} {'|S2|/med':>10} {'cond':>10} "
          f"{'||eps||':>10} {'node err':>10} {'share of excess':>16}")
    for i, k in enumerate(order[:10]):
        print(f"  {i:>4} {k:>5} {grid.lam[k]:>9.1f} {nrm[k]/med:>10.5f} "
              f"{cond[k]:>10.3e} {epsn[k]:>10.3e} {err[k]:>10.3e} "
              f"{exc[k]/(tot-base):>16.4f}")
    cum = np.cumsum(exc[order]) / (tot - base)
    for m in (1, 3, 10, 30, 100):
        print(f"  top-{m:<4d} nodes carry {cum[m-1]:.4f} of the excess")
