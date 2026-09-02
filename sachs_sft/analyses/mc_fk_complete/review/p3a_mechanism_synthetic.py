"""PART 3a -- the amplification mechanism, isolated from the data.

Claim under test.  Write B = Sigma2 + delta and whiten with Sigma2:

    eps = Sigma2^{-1/2} delta Sigma2^{-1/2},      ||eps|| <= (||delta||/||Sigma2||) * cond(Sigma2)

Because  Q = solve_Q(Sigma2, zeta) ~ zeta / d^2  in the Sigma2 eigenbasis, the
cumulant the field actually receives is

    cum3(Q, B) = cum3(Q, Sigma2 (I + eps))  =  zeta * [1 + O(eps)]^2 ,

i.e. the deformation error is controlled by the mismatch measured in the
Sigma2 METRIC, not by the mismatch in norm.  A relative norm mismatch nu is
therefore amplified by up to cond(Sigma2), and it enters a three-leg product,
so an O(1) whitened mismatch already gives an O(few) error.

Test: dial cond(M) = kappa at fixed relative perturbation nu and check that the
injection error tracks nu*kappa; then project the small eigendirection out with
solve_Q's rcond and check that the amplification disappears.
"""
import sys
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                                    # noqa: E402

_DS = ex._DS
rng = np.random.default_rng(11)
n = 6


def sym3(T):
    return sum(T.transpose(p) for p in
               [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]) / 6.0


def one_draw(kappa, nu, rcond, rs):
    U = np.linalg.qr(rs.normal(size=(n, n)))[0]
    d = np.ones(n)
    d[0] = 1.0 / kappa                       # one near-null direction
    M = U @ np.diag(d) @ U.T
    M = 0.5 * (M + M.T)
    Z = sym3(rs.normal(size=(n, n, n)))
    S = rs.normal(size=(n, n))
    S = 0.5 * (S + S.T)
    S /= np.abs(S).max()
    delta = nu * np.abs(M).max() * S
    B = M + delta
    Q = _DS.solve_Q(M, Z, rcond=rcond)
    inj_M = _DS.cum3_from_Q(Q, M)            # what Q was designed to inject
    inj_B = _DS.cum3_from_Q(Q, B)            # what the field actually receives
    err = np.abs(inj_B - inj_M).max() / max(np.abs(inj_M).max(), 1e-300)
    # the whitened mismatch
    w, Uv = np.linalg.eigh(M)
    keep = w > rcond * w.max()
    Wi = Uv[:, keep] @ np.diag(w[keep] ** -0.5) @ Uv[:, keep].T
    eps = Wi @ delta @ Wi
    return err, float(np.abs(np.linalg.eigvalsh(eps)).max())


print("=" * 92)
print("3a-i  injection error vs cond(M) at fixed relative perturbation nu "
      "(rcond = 1e-8, nothing truncated)")
print("=" * 92)
print(f"{'nu':>8} " + " ".join(f"{'k=%g' % k:>11}" for k in
                               (1, 3, 10, 30, 100, 300, 1000)))
for nu in (0.01, 0.03, 0.1):
    med, epsm = [], []
    for kappa in (1, 3, 10, 30, 100, 300, 1000):
        rs = np.random.default_rng(7)
        vals = [one_draw(kappa, nu, 1e-8, rs) for _ in range(40)]
        med.append(np.median([v[0] for v in vals]))
        epsm.append(np.median([v[1] for v in vals]))
    print(f"{nu:>8.2f} " + " ".join(f"{m:>11.4f}" for m in med) + "   <- rel. injection error")
    print(f"{'':>8} " + " ".join(f"{m:>11.4f}" for m in epsm) + "   <- ||eps||  (median)")
    print(f"{'':>8} " + " ".join(f"{nu*k:>11.4f}" for k in
                                 (1, 3, 10, 30, 100, 300, 1000)) + "   <- nu * cond(M)")
    print()

print("=" * 92)
print("3a-ii  same, but the near-null direction is projected out via rcond")
print("=" * 92)
nu = 0.1
print(f"nu = {nu}")
print(f"{'rcond':>10} " + " ".join(f"{'k=%g' % k:>11}" for k in (1, 10, 100, 1000)))
for rcond in (1e-8, 1e-4, 1e-3, 1e-2, 1e-1, 0.5):
    row = []
    for kappa in (1, 10, 100, 1000):
        rs = np.random.default_rng(7)
        vals = [one_draw(kappa, nu, rcond, rs) for _ in range(40)]
        row.append(np.median([v[0] for v in vals]))
    print(f"{rcond:>10.0e} " + " ".join(f"{m:>11.4f}" for m in row))
print()
print("reading: the row rcond=1e-8 keeps the near-null direction and the error")
print("grows ~ nu*kappa; once rcond > 1/kappa that direction is removed from Q")
print("and the error collapses to ~nu, independent of kappa.")
