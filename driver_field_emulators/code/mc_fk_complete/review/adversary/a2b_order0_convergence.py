"""A2 continued: is the 3% baseline offset of the discrete Order-0 a DISCRETISATION
artefact (must vanish as dlam -> 0) or a convention error (must not)?

Cov(k1_A, k1_B) = dlam^2 [ sum_j Wd_j^2 A_j + 2 sum_j Wd_j T_j ],
    T_j = sum_{p<j} rho^{j-p} Wd_p A_p,   T_j = rho (T_{j-1} + Wd_{j-1} A_{j-1})
-- an O(N) recursion, exact for the AR(1) chain, no truncation.
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex

DS, BG, CORE = ex._DS, ex._BG, ex._CORE


def sigma2_nodes(grid, cg):
    return np.array([CORE._assemble_6x6(grid.builder, cg, float(l)) for l in grid.lam])


def o0_disc(grid, S2, sigma, Wd=None):
    N = grid.lam.size; dlam = grid.dlam
    Wd = grid.Wd if Wd is None else Wd
    V = S2 / (2.0 * sigma)
    rho = float(np.exp(-dlam / sigma))
    A = np.empty_like(V); A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1 - rho * rho) * V[k]
    tot = np.einsum("j,jAB->AB", Wd ** 2, A, optimize=True)
    T = np.zeros((6, 6)); acc = np.zeros((6, 6))
    for j in range(1, N):
        T = rho * (T + Wd[j - 1] * A[j - 1])
        acc += Wd[j] * T
    return dlam * dlam * (tot + 2.0 * acc)


def o0_white(grid, S2, Wd=None):
    Wd = grid.Wd if Wd is None else Wd
    return np.einsum("j,jAB->AB", grid.dlam * Wd ** 2, S2, optimize=True)


ANC = ("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
       "sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz")
ad = np.load(ANC, allow_pickle=True)
m = (ad["a"] == 0) & (ad["b"] == 0) & (ad["order"] == 0)
gam_a, cos_a = [], []
for xi, yi in zip(ad["x"][m], ad["y"][m]):
    xi = np.asarray(xi, float); yi = np.asarray(yi, float)
    c = float(np.clip(np.dot(xi / np.linalg.norm(xi), yi / np.linalg.norm(yi)), -1, 1))
    gam_a.append(math.degrees(math.acos(c)) * 60.0); cos_a.append(c)
o = np.argsort(gam_a)
gam_a = np.asarray(gam_a)[o]; cos_a = np.asarray(cos_a)[o]
vs = np.asarray(ad["value"][m], float)[o]
tfin = float(ad["t_final"][0])

targets = [1.0, 5.0, 17.3]
print("=" * 96)
print("A2.4  WHITE-NOISE (sigma -> 0) discrete Order-0 / driver_stats.order0_mc  vs  N")
print("      1/N convergence => the offset is discretisation, not a wrong convention")
print("=" * 96)
print(f"{'N':>7} {'dlam':>8}" + "".join(f"{g:>12.1f}'" for g in targets))
prev = None
rows = {}
for N in (250, 500, 1000, 2000, 4000):
    grid = ex.build_grid(n_lambda=N, apply_anchor=True)
    cells = []
    for gt in targets:
        i = int(np.argmin(np.abs(gam_a - gt))); cg = float(cos_a[i])
        S2 = sigma2_nodes(grid, cg)
        o0mc = DS.order0_mc(cg, lam_f=tfin, builder=grid.builder)
        cells.append(o0_white(grid, S2)[0, 3] / o0mc)
    rows[N] = cells
    print(f"{N:>7} {grid.dlam:8.4f}" + "".join(f"{c:>13.5f}" for c in cells))
print("  successive (ratio-1) halving check  [1/N => 0.50]:")
Ns = sorted(rows)
for a, b in zip(Ns[:-1], Ns[1:]):
    print(f"    N {a}->{b}: " + "  ".join(f"{(rows[b][k]-1)/(rows[a][k]-1):.3f}"
                                          for k in range(len(targets))))
print("  Richardson (2 r(2N) - r(N)) from the last pair: " +
      "  ".join(f"{2*rows[Ns[-1]][k]-rows[Ns[-2]][k]:.5f}" for k in range(len(targets))))

print()
print("=" * 96)
print("A2.5  FINITE sigma_lambda, at FIXED sigma/dlam (so the lattice error is held fixed)")
print("      discrete AR(1) Order-0 / order0_mc.  Tests V = Sigma2/(2 sigma) + the AR(1)")
print("      recursion + the 'smearing -> Sigma2' claim, none of which the fold sees.")
print("=" * 96)
for ratio in (4.19, 8.39):
    print(f"  sigma/dlam = {ratio}")
    print(f"{'':>4}{'N':>7} {'sigma':>8}" + "".join(f"{g:>12.1f}'" for g in targets))
    for N in (500, 1000, 2000, 4000):
        grid = ex.build_grid(n_lambda=N, apply_anchor=True)
        sig = ratio * grid.dlam
        cells = []
        for gt in targets:
            i = int(np.argmin(np.abs(gam_a - gt))); cg = float(cos_a[i])
            S2 = sigma2_nodes(grid, cg)
            o0mc = DS.order0_mc(cg, lam_f=tfin, builder=grid.builder)
            cells.append(o0_disc(grid, S2, sig)[0, 3] / o0mc)
        print(f"{'':>4}{N:>7} {sig:8.3f}" + "".join(f"{c:>13.5f}" for c in cells))

print()
print("=" * 96)
print("A2.6  PRODUCTION SETTING (N = 1000): sigma sweep, ratio to order0_mc and to analysis-3")
print("=" * 96)
grid = ex.build_grid(n_lambda=1000, apply_anchor=True)
print(f"  dlam = {grid.dlam:.4f}")
for gt in targets:
    i = int(np.argmin(np.abs(gam_a - gt))); cg = float(cos_a[i])
    S2 = sigma2_nodes(grid, cg)
    o0mc = DS.order0_mc(cg, lam_f=tfin, builder=grid.builder)
    line = [f"  gamma={gam_a[i]:7.3f}'  o0mc={o0mc:.5e} an3={vs[i]:.5e} (o0mc/an3={o0mc/vs[i]:.4f})"]
    for sig in (32.0, 16.0, 8.0, 4.0, 2.0):
        r = o0_disc(grid, S2, sig)[0, 3] / o0mc
        line.append(f"    sigma={sig:5.1f} (sigma/dlam={sig/grid.dlam:5.2f})  disc/o0mc = {r:.5f}")
    line.append(f"    sigma-> 0 (white formula)                     disc/o0mc = "
                f"{o0_white(grid,S2)[0,3]/o0mc:.5f}")
    print("\n".join(line))
