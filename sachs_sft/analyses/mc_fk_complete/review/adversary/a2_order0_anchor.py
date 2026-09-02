"""ADVERSARIAL TEST A2 -- the Order-0 anchor.

THE ATTACK: "MC and exact expectation agree only because they share the author's
conventions."  The shared conventions are:  the trapezoid observable weights
w_0 = w_{N-1} = dlam/2;  the response exponent (D_j/D_k)^2;  the LOS window
Wd[j] = sum_{k>=j} w_k (D_j/D_k)^2;  the drive measure  s += f dlam;  the node
variance V = Sigma2/(2 sigma_lambda);  the AR(1) recursion; the two-ray index
layout with ray 0 = 0..2, ray 1 = 3..5, and the cross-ray block [0:3, 3:6].

EVERY ONE OF THOSE also fixes the GAUSSIAN Order-0 two-point function, which
involves no F vertex, no Q, no zeta, no placement split, and no Wick algebra of
the estimator.  Its value is fixed independently by

    O0(gamma) = int Sigma2_00(cos; la) G(la)^2 dla,   G = int_la^lam_f (D_la/D_t)^2 dt

(driver_stats.order0_mc, evaluated by 48-node Gauss-Legendre on a continuum grid)
and, one step further out, by the sft-wick analysis-3 workflow
sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz, which shares NOTHING with
this folder except the corr_op propagator itself.

So: run the SAME discrete conventions on the Order-0 observable and see whether
they land on the independent number.  Then break each convention and check the
comparison notices.
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex

DS, BG, CORE = ex._DS, ex._BG, ex._CORE
LAMF = float(BG.LAM_SOURCE_BASELINE)


def sigma2_nodes(grid, cos_gamma):
    return np.array([CORE._assemble_6x6(grid.builder, cos_gamma, float(l))
                     for l in grid.lam])


def wd_from(D, w):
    tail = np.cumsum((w / D ** 2)[::-1])[::-1]
    return D ** 2 * tail


def order0_discrete(grid, S2, sigma, Wd=None, expo=2.0, chunk=200):
    """Cov(kappa1_A, kappa1_B) on the discrete grid with the AR(1) chain.

    = sum_{j,p} dlam^2 Wd_j Wd_p rho^|j-p| A[min(j,p)]
    """
    N = grid.lam.size
    dlam = grid.dlam
    if Wd is None:
        Wd = grid.Wd
    V = S2 / (2.0 * sigma)
    rho = float(np.exp(-dlam / sigma))
    A = np.empty_like(V); A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1 - rho * rho) * V[k]
    idx = np.arange(N)
    out = np.zeros((6, 6))
    for lo in range(0, N, chunk):
        ps = np.arange(lo, min(lo + chunk, N))
        lag = np.abs(idx[:, None] - ps[None, :])
        R = (rho ** lag)[:, :, None, None] * A[np.minimum(idx[:, None], ps[None, :])]
        out += dlam * dlam * np.einsum("j,p,jpAB->AB", Wd, Wd[ps], R, optimize=True)
    return out


def order0_white(grid, S2, Wd=None):
    """sigma_lambda -> 0 limit of the same discrete object: sum_j dlam Wd_j^2 Sigma2[j]."""
    if Wd is None:
        Wd = grid.Wd
    return np.einsum("j,jAB->AB", grid.dlam * Wd ** 2, S2, optimize=True)


# ---------------------------------------------------------------- anchor data
ANC = ("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
       "sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz")
ad = np.load(ANC, allow_pickle=True)
m = (ad["a"] == 0) & (ad["b"] == 0) & (ad["order"] == 0)
xs, ys, vs = ad["x"][m], ad["y"][m], ad["value"][m]
gam_a, cos_a = [], []
for xi, yi in zip(xs, ys):
    xi = np.asarray(xi, float); yi = np.asarray(yi, float)
    c = float(np.clip(np.dot(xi / np.linalg.norm(xi), yi / np.linalg.norm(yi)), -1, 1))
    gam_a.append(math.degrees(math.acos(c)) * 60.0); cos_a.append(c)
o = np.argsort(gam_a)
gam_a = np.asarray(gam_a)[o]; cos_a = np.asarray(cos_a)[o]; vs = np.asarray(vs, float)[o]
tfin = float(ad["t_final"][0])
print(f"analysis-3 O0 anchor: {len(gam_a)} gammas, t_final={tfin}")

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
grid = ex.build_grid(n_lambda=N, apply_anchor=True)     # match order0_mc's builder
gridR = ex.build_grid(n_lambda=N, apply_anchor=False)   # what the FK code uses
print(f"grid N={N}  dlam={grid.dlam:.4f}  lam[0]={grid.lam[0]}  lam[-1]={grid.lam[-1]}")
print(f"order0_mc builder apply_c0 default -> ANCHOR_C0={DS.ANCHOR_C0:.6f}")

targets = [0.5, 1.0, 2.0, 5.0, 17.3]
print()
print("=" * 100)
print("A2.1  discrete (FK-code conventions) Order-0 cross-ray kk  vs  the two independent targets")
print("=" * 100)
hdr = (f"{'gamma':>8} {'sig':>5} {'discrete[0,3]':>14} {'order0_mc':>13} {'disc/o0mc':>10}"
       f" {'analysis3':>13} {'disc/an3':>9}")
print(hdr)
store = {}
for gt in targets:
    i = int(np.argmin(np.abs(gam_a - gt)))
    cg = float(cos_a[i]); g = float(gam_a[i])
    S2 = sigma2_nodes(grid, cg)
    o0mc = DS.order0_mc(cg, lam_f=tfin, builder=grid.builder)
    row = []
    for sig in (8.0, 4.0, 1.0):
        d = order0_discrete(grid, S2, sig)
        row.append(d[0, 3])
        print(f"{g:8.3f} {sig:5.1f} {d[0,3]:14.6e} {o0mc:13.6e} {d[0,3]/o0mc:10.5f}"
              f" {vs[i]:13.6e} {d[0,3]/vs[i]:9.5f}")
    dw = order0_white(grid, S2)
    print(f"{g:8.3f} {'0':>5} {dw[0,3]:14.6e} {o0mc:13.6e} {dw[0,3]/o0mc:10.5f}"
          f" {vs[i]:13.6e} {dw[0,3]/vs[i]:9.5f}   <- white-noise limit")
    store[gt] = dict(cos=cg, gam=g, S2=S2, o0mc=o0mc, an3=vs[i], white=dw[0, 3])
    print()

print("=" * 100)
print("A2.2  DOES IT HAVE TEETH?  break one shared convention at a time (white-noise limit,")
print("      ratio to driver_stats.order0_mc; a convention that does not matter shows 1.000)")
print("=" * 100)
D = grid.D; dlam = grid.dlam
w_trap = np.full(N, dlam); w_trap[0] = w_trap[-1] = 0.5 * dlam
w_flat = np.full(N, dlam)
variants = {}
variants["baseline (as coded)"] = grid.Wd
variants["drop trapezoid end weights"] = wd_from(D, w_flat)
tail1 = np.cumsum((w_trap / D ** 1)[::-1])[::-1]; variants["response exponent 2 -> 1"] = D ** 1 * tail1
tail4 = np.cumsum((w_trap / D ** 4)[::-1])[::-1]; variants["response exponent 2 -> 4"] = D ** 4 * tail4
variants["window reversed (sum_{k<=j})"] = D ** 2 * np.cumsum(w_trap / D ** 2)
variants["Wd -> Wd (no window at all, =1)"] = np.ones(N)
print(f"{'variant':<38}" + "".join(f"{g:>12.1f}'" for g in targets))
for name, Wd in variants.items():
    cells = []
    for gt in targets:
        st = store[gt]
        cells.append(order0_white(grid, st["S2"], Wd=Wd)[0, 3] / st["o0mc"])
    print(f"{name:<38}" + "".join(f"{c:>13.4f}" for c in cells))

print()
print("      also: the drive/variance normalisation (finite sigma, gamma=1'):")
st = store[1.0]
for name, fac in (("V = Sigma2/(2 sigma)  [as coded]", 1.0),
                  ("V = Sigma2            [x 2 sigma]", None),
                  ("V = Sigma2/sigma      [x 2]", 2.0)):
    for sig in (8.0, 4.0):
        f = (2.0 * sig) if fac is None else fac
        d = order0_discrete(grid, st["S2"] * f, sig)
        print(f"        {name:<36} sigma={sig:4.1f}  ratio to order0_mc = {d[0,3]/st['o0mc']:.4f}")

print()
print("=" * 100)
print("A2.3  index layout: full 6x6 Order-0 (white limit) at gamma = 1', normalised by [0,0]")
print("=" * 100)
st = store[1.0]
M = order0_white(grid, st["S2"])
np.set_printoptions(precision=4, suppress=False, linewidth=160)
print(M / M[0, 0])
print(f"  [0,3] (cross-ray kk) = {M[0,3]:.6e}   [0,0] (auto kk) = {M[0,0]:.6e}")
print(f"  order0_mc (cross-ray kk)                = {st['o0mc']:.6e}")
print(f"  order0_mc at cos=1 (auto, for reference) = {DS.order0_mc(1.0, lam_f=tfin, builder=grid.builder):.6e}")
