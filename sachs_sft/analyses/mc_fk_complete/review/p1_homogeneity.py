"""PART 1 -- the two calibrations coincide in the white-noise limit.

Two statements, proved algebraically in REPORT.md and checked here numerically:

(A) ``solve_Q`` is jointly homogeneous:  solve_Q(a M, b Z) = (b/a^2) solve_Q(M, Z)
    for any a > 0 and any b.  Hence, with a = 2 sigma and b = 1,

        solve_Q(Sigma2, Z) = solve_Q(2 sigma V, Z) = solve_Q(V, Z/(2 sigma)^2)

    EXACTLY -- the two ways of writing the "nominal" calibration are one object.

(B) The only difference between the two calibrations is the matrix Q is solved
    against: Sigma2 = 2 sigma V (nominal) vs B (smeared).  Substituting B -> 2 sigma V
    in the smeared branch must reproduce the nominal Q bit-for-bit.  And
    B -> Sigma2 in the joint limit dlam/sigma -> 0, sigma/L_Sigma -> 0, so the
    two calibrations share the same white-noise target.
"""
import sys
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                                   # noqa: E402

_DS = ex._DS
rng = np.random.default_rng(20260828)


def rel(x, y):
    d = np.abs(x - y).max()
    s = max(np.abs(x).max(), np.abs(y).max(), 1e-300)
    return d / s


print("=" * 78)
print("PART 1A -- homogeneity of solve_Q  (random SPD M, random symmetric Z)")
print("=" * 78)
# random well-conditioned SPD M and a fully symmetric Z, n = 6 like the real one
n = 6
G = rng.normal(size=(n, n))
M = G @ G.T + 0.3 * np.eye(n)
Zr = rng.normal(size=(n, n, n))
Z = sum(Zr.transpose(p) for p in
        [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]) / 6.0

Q0 = _DS.solve_Q(M, Z)
print(f"{'a':>10} {'b':>10} {'rel err vs (b/a^2) Q0':>24}")
for a, b in [(2.0, 1.0), (0.5, 1.0), (16.0, 1.0), (1.0, 7.0), (1e3, 1e-4),
             (1e-3, 1e5), (37.0, 37.0 ** 3)]:
    Qab = _DS.solve_Q(a * M, b * Z)
    print(f"{a:>10.4g} {b:>10.4g} {rel(Qab, (b / a ** 2) * Q0):>24.3e}")

print()
print("the identity actually used:  solve_Q(2s V, Z)  ==  solve_Q(V, Z/(2s)^2)")
print(f"{'sigma':>10} {'rel err':>12}")
for sig in (0.5, 2.0, 8.0, 32.0, 1000.0):
    V = M / (2.0 * sig)
    lhs = _DS.solve_Q(2.0 * sig * V, Z)
    rhs = _DS.solve_Q(V, Z / (2.0 * sig) ** 2)
    print(f"{sig:>10.1f} {rel(lhs, rhs):>12.3e}")

print()
print("=" * 78)
print("PART 1B -- on the REAL Sigma2/zeta6, and the B -> 2 sigma V substitution")
print("=" * 78)
N = 400
grid = ex.build_grid(n_lambda=N)
for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    print(f"\ngamma = {gam}'   (N = {N}, dlam = {grid.dlam:.4f})")
    print(f"{'sigma':>7} {'cond(Sigma2)':>13} {'Q_nom vs Q(V,Z/(2s)^2)':>24} "
          f"{'Q_smeared[B:=2sV] vs Q_nom':>28}")
    for sig in (2.0, 8.0, 32.0):
        V, A, Q_nom, rho = ex.node_stats(grid, cosg, sig, calibrate="nominal")
        S2 = V * (2.0 * sig)
        # the "written the other way" nominal
        Q_alt = np.array([_DS.solve_Q(V[k], np.zeros((6, 6, 6))) for k in range(1)])
        Zt = np.array([_DS.zeta6(cosg, float(l)) for l in grid.lam])
        Q_alt = np.array([_DS.solve_Q(V[k], Zt[k] / (2.0 * sig) ** 2)
                          for k in range(N)])
        # the smeared branch with B replaced by 2 sigma V  ==  Sigma2
        Q_sub = np.array([_DS.solve_Q(S2[k], Zt[k]) for k in range(N)])
        cnd = np.linalg.cond(S2[N // 2])
        print(f"{sig:>7.1f} {cnd:>13.4e} {rel(Q_nom, Q_alt):>24.3e} "
              f"{rel(Q_sub, Q_nom):>28.3e}")

print()
print("=" * 78)
print("PART 1C -- how B approaches Sigma2")
print("=" * 78)
print("decomposition:  B[k] = dlam sum_l rho^|k-l| A[min(k,l)]")
print("  for CONSTANT V the interior value is  Sigma2 * (r/2) coth(r/2),  r = dlam/sigma")
print("  so B/Sigma2 -> 1 needs (i) r -> 0 and (ii) sigma small vs the scale on")
print("  which Sigma2 varies (the exponential window smears lambda).")
print()
gam = 1.0
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
print(f"gamma = {gam}'")
print(f"{'N':>6} {'dlam':>8} {'sigma':>7} {'r=dl/s':>8} {'(r/2)coth':>10} "
      f"{'B/S2 mid':>10} {'|B-S2|/|S2| mid':>16} {'interior med':>13}")
for N_, sig in [(1000, 32.0), (1000, 16.0), (1000, 8.0), (1000, 4.0), (1000, 2.0),
                (1000, 1.0), (2000, 8.0), (4000, 8.0),
                (500, 16.0), (1000, 8.0), (1999, 4.0), (3997, 2.0)]:
    g = ex.build_grid(n_lambda=N_)
    V, A, Q, rho = ex.node_stats(g, cosg, sig, calibrate="nominal")
    B = ex.smeared_covariance(g, A, rho)
    S2 = V * (2.0 * sig)
    r = g.dlam / sig
    coth = (r / 2.0) / np.tanh(r / 2.0)
    k = N_ // 2
    ratio = np.abs(B[k]).max() / np.abs(S2[k]).max()
    relmid = np.abs(B[k] - S2[k]).max() / np.abs(S2[k]).max()
    # interior only (drop 5 correlation lengths at each end)
    nedge = max(1, int(5 * sig / g.dlam))
    sl = slice(nedge, N_ - nedge)
    med = np.median([np.abs(B[i] - S2[i]).max() / np.abs(S2[i]).max()
                     for i in range(nedge, N_ - nedge, max(1, (N_ - 2 * nedge) // 50))])
    print(f"{N_:>6} {g.dlam:>8.4f} {sig:>7.1f} {r:>8.4f} {coth:>10.6f} "
          f"{ratio:>10.6f} {relmid:>16.4e} {med:>13.4e}")
