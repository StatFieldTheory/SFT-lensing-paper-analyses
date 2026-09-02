"""ADVERSARIAL TEST A5 -- the chain-level cumulant, by sampling, with NO F vertex.

<kappa1_a kappa1_b kappa1_c>_c  must equal  sum_j dlam Wd_j^3 zeta_j  in the
white-noise limit.  This exercises, all at once and against the TABULATED zeta:
    the drive measure  s += f dlam        the trapezoid observable weights
    the window Wd and its (D_j/D_k)^2     V = Sigma2/(2 sigma) and the AR(1) chain
    the two-ray index layout              the 0.5, solve_Q, and the calibration matrix
It uses NO F vertex, NO Q-placement split and NO Wick contraction of the estimator,
so it cannot share a bug with the FK assembly under test.

Variance control: kappa1[f] = kappa1[z] + eps kappa1[g] exactly (kappa1 is linear in
the drive), so the O(eps) third moment is estimated as
    < Kg_a Kz_b Kz_c + Kz_a Kg_b Kz_c + Kz_a Kz_b Kg_c >
which is a plain sample mean of a quartic in Gaussians -- no cancellation of a
large O(eps^0) piece.
"""
import sys, os, math, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex
from sachs_mc_core import _cholesky_psd

DS = ex._DS
N = int(os.environ.get("NLAM", 1000))
SIG = float(os.environ.get("SIG", 8.0))
GAM = 1.01546
M_BATCH = 40_000
NBATCH = 6
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))

grid = ex.build_grid(n_lambda=N)
V, A, Q, rho = ex.node_stats(grid, cosg, SIG, calibrate="smeared")
Vn, An, Qn, _ = ex.node_stats(grid, cosg, SIG, calibrate="nominal")
zt = np.array([DS.zeta6(cosg, float(l)) for l in grid.lam])
zs = (zt + zt.transpose(0, 1, 3, 2) + zt.transpose(0, 2, 1, 3)
      + zt.transpose(0, 2, 3, 1) + zt.transpose(0, 3, 1, 2)
      + zt.transpose(0, 3, 2, 1)) / 6.0
target = np.einsum("j,jabc->abc", grid.dlam * grid.Wd ** 3, zs, optimize=True)
print(f"N={N} dlam={grid.dlam:.4f} sigma={SIG} gamma={GAM}'  rho={rho:.6f}")
print(f"target sum_j dlam Wd^3 zeta_sym : (000) {target[0,0,0]:.6e}  "
      f"(003) {target[0,0,3]:.6e}  (033) {target[0,3,3]:.6e}  (011) {target[0,1,1]:.6e}")


def run(Quse, Csub, seed, m=M_BATCH):
    rng = np.random.default_rng(seed)
    Lf = np.array([_cholesky_psd(V[k]) for k in range(N)])
    Qf = Quse.reshape(N, 6, 36)
    s1mr = math.sqrt(max(1 - rho * rho, 0.0))
    w = np.full(N, grid.dlam); w[0] = w[-1] = 0.5 * grid.dlam
    z = np.zeros((6, m)); sz = np.zeros((6, m)); sg = np.zeros((6, m))
    Kz = np.zeros((6, m)); Kg = np.zeros((6, m))
    for k in range(N):
        innov = Lf[k] @ rng.standard_normal((6, m))
        z = innov if k == 0 else rho * z + s1mr * innov
        ww = np.einsum("im,jm->ijm", z, z, optimize=True) - Csub[k][:, :, None]
        g = 0.5 * (Qf[k] @ ww.reshape(36, m))
        sz = grid.resp[k] * sz + z * grid.dlam
        sg = grid.resp[k] * sg + g * grid.dlam
        Kz += w[k] * sz; Kg += w[k] * sg
    return Kz, Kg


def third(Kz, Kg, idx):
    a, b, c = idx
    per = Kg[a]*Kz[b]*Kz[c] + Kz[a]*Kg[b]*Kz[c] + Kz[a]*Kz[b]*Kg[c]
    return float(per.mean()), float(per.std()) / math.sqrt(per.size)


IDX = ((0, 0, 0), (0, 0, 3), (0, 3, 3), (0, 1, 1))
for label, Quse, Csub in (("smeared calibration, C = A (as coded)", Q, A),
                          ("nominal calibration, C = A", Qn, A),
                          ("smeared calibration, C = V (wrong constant)", Q, V)):
    acc = {i: [] for i in IDX}
    t0 = time.time()
    for b in range(NBATCH):
        Kz, Kg = run(Quse, Csub, 900_000 + b)
        for i in IDX:
            acc[i].append(third(Kz, Kg, i)[0])
    print(f"\n  {label}   [{NBATCH} x {M_BATCH:,} real, {time.time()-t0:.0f}s]")
    for i in IDX:
        v = np.array(acc[i]); mu = v.mean(); se = v.std(ddof=1)/math.sqrt(len(v))
        print(f"    {str(i):>10}  measured {mu:14.6e} +- {se:8.2e}   "
              f"ratio to target = {mu/target[i]:7.4f} +- {se/abs(target[i]):.4f}")
