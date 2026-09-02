"""ADVERSARIAL TEST A6 -- gate 4 redone at EXACT fold nodes (no gamma interpolation),
plus the discretisation conventions that the sigma->0 limit still carries.

run_gate2.analytic_fk interpolates the fold LINEARLY in gamma between table nodes that
are 27% apart in gamma, on a curve that is convex in log-log.  Evaluate at the nodes
themselves so that choice cannot flatter (or spoil) the answer.
"""
import sys, os, math, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex
import _bootstrap

DS = ex._DS
d = np.load(_bootstrap.FOLD, allow_pickle=True)
m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
gs, vsf = [], np.asarray(d["value"], float)[m]
for x, y in zip(d["x"][m], d["y"][m]):
    x = np.asarray(x, float); y = np.asarray(y, float)
    gs.append(math.degrees(math.acos(float(np.clip(
        np.dot(x/np.linalg.norm(x), y/np.linalg.norm(y)), -1, 1)))) * 60)
gs = np.asarray(gs); o = np.argsort(gs); gs = gs[o]; vsf = vsf[o]

N = 1000
grid = ex.build_grid(n_lambda=N)
PICK = [0, 1, 3, 5, 7, 9, 11, 13, 15]
print(f"N={N}  dlam={grid.dlam:.4f}   sigma-> 0 extrapolation from sigma = 8 and 4 "
      f"(linear, as in NOTES)")
print("=" * 116)
print(f"{'gamma[node]':>12} {'fold':>13} {'ref0_inj':>12} {'ref0/fold':>10}"
      f" {'ex(s=8)/fold':>13} {'ex(s=4)/fold':>13} {'extrap/fold':>12} {'T1share':>8}")
print("=" * 116)
rows = []
for i in PICK:
    g = float(gs[i]); fold = float(vsf[i])
    cosg = float(np.cos(np.deg2rad(g / 60.0)))
    out = {}
    for sig in (8.0, 4.0):
        V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate="smeared")
        M = ex.calib_matrix(grid, V, A, rho, sig, "smeared")
        zi = ex.zeta_injected(M, Q)
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
        tot = (T1 + T2) + (T1 + T2).T
        out[sig] = (float(tot[0, 3]), float(ex.fk_reference(grid, zi)[0, 3]),
                    float((T1 + T1.T)[0, 3]) / float(tot[0, 3]))
    e8, r0, sh = out[8.0]
    e4 = out[4.0][0]
    extrap = 2 * e4 - e8
    rows.append((g, fold, r0, e8, e4, extrap, sh))
    print(f"{g:12.5f} {fold:13.6e} {r0:12.6e} {r0/fold:10.5f}"
          f" {e8/fold:13.5f} {e4/fold:13.5f} {extrap/fold:12.5f} {sh:8.4f}")
print("=" * 116)
arr = np.array([[r[5]/r[1], r[2]/r[1]] for r in rows])
print(f"  sigma->0 extrapolated exact / fold : median {np.median(arr[:,0]):.5f}  "
      f"range [{arr[:,0].min():.5f}, {arr[:,0].max():.5f}]")
print(f"  local reference (injected) / fold  : median {np.median(arr[:,1]):.5f}  "
      f"range [{arr[:,1].min():.5f}, {arr[:,1].max():.5f}]")

print()
print("=" * 116)
print("A6.2  discretisation conventions the sigma->0 limit still carries "
      "(ratio of fk_reference variants to the baseline, gamma = 1.01546')")
print("=" * 116)
g = float(gs[3]); cosg = float(np.cos(np.deg2rad(g / 60.0)))
zt = np.array([DS.zeta6(cosg, float(l)) for l in grid.lam])
zs = (zt + zt.transpose(0,1,3,2) + zt.transpose(0,2,1,3) + zt.transpose(0,2,3,1)
      + zt.transpose(0,3,1,2) + zt.transpose(0,3,2,1)) / 6.0


def ref_kernel(D, w, hd_shift):
    n = D.size
    tail = np.cumsum((w / D**2)[::-1])[::-1]
    Wd = D**2 * tail
    dl = grid.dlam
    if hd_shift == "old":       # F at the OLD state: D_{j-1}
        gj = dl * Wd[1:] / D[:-1]**4
        t4 = np.zeros(n); t4[:-1] = np.cumsum(gj[::-1])[::-1]
    else:                       # F at the NEW state: D_j
        gj = dl * Wd[1:] / D[1:]**4
        t4 = np.zeros(n); t4[:-1] = np.cumsum(gj[::-1])[::-1]
    Hd = D**4 * t4
    kern = dl * Wd * Hd
    M = np.einsum("m,Abc,mbcB->AB", kern, ex.F6, zs, optimize=True)
    return float((M + M.T)[0, 3])


fold = float(vsf[3])
for nn in (500, 1000, 2000, 4000, 8000):
    gr = ex.build_grid(n_lambda=nn)
    ztn = np.array([DS.zeta6(cosg, float(l)) for l in gr.lam])
    zsn = (ztn + ztn.transpose(0,1,3,2) + ztn.transpose(0,2,1,3) + ztn.transpose(0,2,3,1)
           + ztn.transpose(0,3,1,2) + ztn.transpose(0,3,2,1)) / 6.0
    w = np.full(nn, gr.dlam); w[0] = w[-1] = 0.5*gr.dlam
    wf = np.full(nn, gr.dlam)
    D = gr.D; dl = gr.dlam

    def kern_of(D, w, shift):
        tail = np.cumsum((w/D**2)[::-1])[::-1]; Wd = D**2*tail
        if shift == "old":
            gj = dl*Wd[1:]/D[:-1]**4
        else:
            gj = dl*Wd[1:]/D[1:]**4
        t4 = np.zeros(nn); t4[:-1] = np.cumsum(gj[::-1])[::-1]
        return dl*Wd*(D**4*t4)

    vals = {}
    for nm, (ww, sh) in (("baseline", (w, "old")), ("F at NEW state", (w, "new")),
                         ("no trapezoid ends", (wf, "old"))):
        k = kern_of(D, ww, sh)
        M = np.einsum("m,Abc,mbcB->AB", k, ex.F6, zsn, optimize=True)
        vals[nm] = float((M+M.T)[0, 3])
    print(f"  N={nn:5d} dlam={gr.dlam:7.4f}  baseline/fold={vals['baseline']/fold:.5f}"
          f"   Fnew/base={vals['F at NEW state']/vals['baseline']:.5f}"
          f"   noTrap/base={vals['no trapezoid ends']/vals['baseline']:.5f}")
