"""PART 3b -- the same mechanism, measured on the real driving field.

Two truncations are applied to solve_Q:

* ``rcond`` -- solve_Q's own argument.  It is a RELATIVE, within-node floor:
  it removes eigendirections with d_p < rcond * d_max(node).  This is the
  amplifier that grows as the two rays become parallel (cond(Sigma2) rises from
  3.4 at 30' to 268 at 1').
* a NODE floor -- zero Q at nodes whose whole Sigma2 is small in ABSOLUTE terms.
  This is needed because the extracted equal-time density Sigma2(lambda) is not
  positive along the whole ray: Sigma2_00 crosses zero near lambda = 1309 and
  dips to ~1% of its neighbours near lambda = 1708, so Q ~ zeta/Sigma2^2 has
  isolated near-singular nodes that no within-node relative floor can see.

Both are switched on independently so the two amplifiers can be told apart.
"""
import sys
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                                      # noqa: E402
from kappa3_norm import kappa1_third_moment, local_prediction, smeared_B  # noqa: E402

_DS = ex._DS
N = 1000
SIG = 8.0
GAMMAS = [1.0, 5.0, 17.3, 30.0]
grid = ex.build_grid(n_lambda=N)


def build(cosg, sig):
    V = np.empty((N, 6, 6))
    Z = np.empty((N, 6, 6, 6))
    for k, lk in enumerate(grid.lam):
        V[k] = ex._CORE._assemble_6x6(grid.builder, cosg, float(lk)) / (2.0 * sig)
        Z[k] = _DS.zeta6(cosg, float(lk))
    rho = float(np.exp(-grid.dlam / sig))
    A = np.empty_like(V)
    A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1 - rho * rho) * V[k]
    S2 = V * (2.0 * sig)
    B = smeared_B(grid, A, rho)
    return V, Z, A, rho, S2, B


def solveQ(M, Z, rcond, node_floor):
    nrm = np.array([np.abs(M[k]).max() for k in range(N)])
    med = float(np.median(nrm))
    keep = nrm >= node_floor * med
    Q = np.zeros((N, 6, 6, 6))
    for k in range(N):
        if keep[k]:
            Q[k] = _DS.solve_Q(M[k], Z[k], rcond=rcond)
    return Q, int(N - keep.sum())


print("=" * 112)
print("3b-i  per-node correlation between the whitened mismatch and the "
      "per-node injection error  (sigma = 8, N = 1000)")
print("=" * 112)
print(f"{'gamma':>7} {'cond(S2) med':>13} {'nu med':>9} {'||eps|| med':>12} "
      f"{'node err med':>13} {'Spearman(log||eps||, log err)':>30}")
for gam in GAMMAS:
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    V, Z, A, rho, S2, B = build(cosg, SIG)
    Q, _ = solveQ(S2, Z, 1e-8, 0.0)
    cond, nu, epsn, err = [], [], [], []
    for k in range(N):
        d, U = np.linalg.eigh(S2[k])
        if d.max() <= 0:
            continue
        good = d > 1e-8 * d.max()
        if good.sum() < 6:
            continue
        Wi = U @ np.diag(d ** -0.5) @ U.T
        e = Wi @ (B[k] - S2[k]) @ Wi
        iM = _DS.cum3_from_Q(Q[k], S2[k])
        iB = _DS.cum3_from_Q(Q[k], B[k])
        sc = max(np.abs(iM).max(), 1e-300)
        cond.append(d.max() / d.min())
        nu.append(np.abs(B[k] - S2[k]).max() / np.abs(S2[k]).max())
        epsn.append(np.abs(np.linalg.eigvalsh(e)).max())
        err.append(np.abs(iB - iM).max() / sc)
    cond, nu, epsn, err = map(np.asarray, (cond, nu, epsn, err))
    from scipy.stats import spearmanr
    rs = spearmanr(np.log(epsn + 1e-300), np.log(err + 1e-300)).statistic
    print(f"{gam:>7.2f} {np.median(cond):>13.4e} {np.median(nu):>9.4f} "
          f"{np.median(epsn):>12.4f} {np.median(err):>13.4f} {rs:>30.4f}")

print()
print("=" * 112)
print("3b-ii  observable-level truncation sweep.  'nom err' = exact FK with the "
      "nominal Q, divided by")
print("       the local FK the SAME (truncated) Q was designed to inject.  "
      "1.0 = faithful injection.")
print("=" * 112)
for gam in GAMMAS:
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    V, Z, A, rho, S2, B = build(cosg, SIG)
    cnd = float(np.linalg.cond(S2[N // 2]))
    print(f"\ngamma = {gam}'   cond(Sigma2) mid-ray = {cnd:.4e}")
    print(f"{'rcond':>9} {'nodefloor':>10} {'#nodes cut':>11} "
          f"{'nom err':>10} {'sm err':>9} {'nom/sm FK':>10} "
          f"{'<k1^3> nom':>11} {'<k1^3> sm':>10}")
    for rcond, nf in [(1e-8, 0.0), (1e-3, 0.0), (1e-2, 0.0), (3e-2, 0.0),
                      (1e-1, 0.0), (3e-1, 0.0),
                      (1e-8, 0.05), (1e-8, 0.2), (1e-8, 0.5),
                      (1e-1, 0.2), (3e-1, 0.5)]:
        Qn, cut = solveQ(S2, Z, rcond, nf)
        Qs, _ = solveQ(B, Z, rcond, nf)
        refn = ex.fk_reference(grid, ex.zeta_injected(S2, Qn))[0, 3]
        refs = ex.fk_reference(grid, ex.zeta_injected(B, Qs))[0, 3]
        T1, T2 = ex.expectation(grid, cosg, SIG, V=V, A=A, Q=Qn, rho=rho)
        fn = ((T1 + T2) + (T1 + T2).T)[0, 3]
        T1, T2 = ex.expectation(grid, cosg, SIG, V=V, A=A, Q=Qs, rho=rho)
        fs = ((T1 + T2) + (T1 + T2).T)[0, 3]
        m3n = kappa1_third_moment(grid, A, Qn, rho)
        m3s = kappa1_third_moment(grid, A, Qs, rho)
        ln = local_prediction(grid, ex.zeta_injected(S2, Qn))
        ls = local_prediction(grid, ex.zeta_injected(B, Qs))
        c = (0, 0, 3)
        print(f"{rcond:>9.0e} {nf:>10.2f} {cut:>11d} "
              f"{fn/refn:>10.4f} {fs/refs:>9.4f} {fn/fs:>10.4f} "
              f"{m3n[c]/ln[c]:>11.4f} {m3s[c]/ls[c]:>10.4f}", flush=True)
