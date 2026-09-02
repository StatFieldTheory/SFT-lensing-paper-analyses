"""A5b.  Two things.

(i) HOW DIFFERENT ARE THE TWO FIELDS?  fk_complete_core simulates with
    Lf = _cholesky_psd(V[k]), which FLOORS the negative eigenvalues of the tabulated
    Sigma2; fk_expect_exact propagates V[k] as given.  So "MC vs exact" compares two
    slightly different stochastic models.  Quantify the floor, and re-run the exact
    expectation on the FLOORED V to see how much of the FK it is worth.

(ii) Redo the chain-level <kappa1^3> test of A5 with the centering constant that
     matches what is actually simulated (the variance recursion of the FLOORED V),
     so that <g> = 0 exactly and no tadpole leaks in.
"""
import sys, os, math, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex
from sachs_mc_core import _cholesky_psd

DS = ex._DS
N, SIG, GAM = 1000, 8.0, 1.01546
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
grid = ex.build_grid(n_lambda=N)
V, A, Q, rho = ex.node_stats(grid, cosg, SIG, calibrate="smeared")

Lf = np.array([_cholesky_psd(V[k]) for k in range(N)])
Vf = np.einsum("kij,klj->kil", Lf, Lf, optimize=True)
dv = np.max(np.abs(Vf - V), axis=(1, 2)) / np.max(np.abs(V), axis=(1, 2))
nbad = int(np.sum(dv > 1e-12))
print("=" * 92)
print("A5b(i)  the MC and the exact expectation do NOT simulate the same field")
print("=" * 92)
print(f"  nodes where _cholesky_psd changes V: {nbad} / {N}")
print(f"  max |V_floored - V| / |V| over nodes = {dv.max():.3e}   (median {np.median(dv):.1e})")
if nbad:
    idx = np.argsort(dv)[-6:]
    for k in idx[::-1]:
        w = np.linalg.eigvalsh(V[k])
        print(f"    k={k:4d} lam={grid.lam[k]:8.1f}  rel change {dv[k]:.3e}  "
              f"eig_min/eig_max = {w.min()/w.max():+.3e}")
Af = np.empty_like(Vf); Af[0] = Vf[0]
for k in range(1, N):
    Af[k] = rho * rho * Af[k - 1] + (1 - rho * rho) * Vf[k]
print(f"  max |A_floored - A| / |A| = "
      f"{np.max(np.max(np.abs(Af-A),axis=(1,2))/np.max(np.abs(A),axis=(1,2))):.3e}")

T1, T2 = ex.expectation(grid, cosg, SIG, V=V, A=A, Q=Q, rho=rho)
base = float(((T1+T2)+(T1+T2).T)[0, 3])
T1f, T2f = ex.expectation(grid, cosg, SIG, V=Vf, A=Af, Q=Q, rho=rho)
flo = float(((T1f+T2f)+(T1f+T2f).T)[0, 3])
print(f"  exact FK with V as given   = {base:.6e}")
print(f"  exact FK with V floored    = {flo:.6e}   ratio = {flo/base:.6f}")
print("  => the model mismatch between the MC and its 'exact expectation' is worth "
      f"{100*abs(flo/base-1):.3f}% of the FK.")

print()
print("=" * 92)
print("A5b(ii)  chain <kappa1^3> with the centering constant that matches the simulation")
print("=" * 92)
zt = np.array([DS.zeta6(cosg, float(l)) for l in grid.lam])
zs = (zt + zt.transpose(0, 1, 3, 2) + zt.transpose(0, 2, 1, 3)
      + zt.transpose(0, 2, 3, 1) + zt.transpose(0, 3, 1, 2)
      + zt.transpose(0, 3, 2, 1)) / 6.0
target = np.einsum("j,jabc->abc", grid.dlam * grid.Wd ** 3, zs, optimize=True)
Qf = Q.reshape(N, 6, 36)
s1mr = math.sqrt(max(1 - rho * rho, 0.0))
w = np.full(N, grid.dlam); w[0] = w[-1] = 0.5 * grid.dlam
IDX = ((0, 0, 0), (0, 0, 3), (0, 3, 3), (0, 1, 1))


def run(Csub, seed, m):
    rng = np.random.default_rng(seed)
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


for label, Csub in (("C = A (from V as given; what the exact side assumes)", A),
                    ("C = A_floored (matches the simulated chain)", Af)):
    acc = {i: [] for i in IDX}
    t0 = time.time()
    NB, MB = 12, 40_000
    for b in range(NB):
        Kz, Kg = run(Csub, 5_500_000 + b, MB)
        for i in IDX:
            a_, b_, c_ = i
            per = Kg[a_]*Kz[b_]*Kz[c_] + Kz[a_]*Kg[b_]*Kz[c_] + Kz[a_]*Kz[b_]*Kg[c_]
            acc[i].append(float(per.mean()))
    print(f"\n  {label}   [{NB} x {MB:,}, {time.time()-t0:.0f}s]")
    for i in IDX:
        v = np.array(acc[i]); mu = v.mean(); se = v.std(ddof=1)/math.sqrt(len(v))
        print(f"    {str(i):>10}  ratio to  sum dlam Wd^3 zeta_sym  = "
              f"{mu/target[i]:7.4f} +- {se/abs(target[i]):.4f}")
print("\n  (the finite-sigma_lambda deficit of this observable is a few per mil at "
      "sigma=8, so the white-noise target is the right yardstick at this precision)")
