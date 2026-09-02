"""ADJUDICATOR's own check: exact expectation vs sigma_lambda at gamma=1', N=1000.

Verifies (i) NOTES gate-2 exact numbers, (ii) the adv_sigma agent's claim that
the regulator can be switched OFF on a fixed lattice so no extrapolation model
is needed, and (iii) the size of the published 2-point extrapolation's bias.
"""
import sys, time, math
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import fk_expect_exact as ex
import _bootstrap

def fold_nodes():
    d = np.load(_bootstrap.FOLD, allow_pickle=True)
    m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x/np.linalg.norm(x), y/np.linalg.norm(y)), -1, 1))))*60)
    g = np.asarray(g); v = np.asarray(d["value"], float)[m]
    o = np.argsort(g)
    return g[o], v[o]

GAM = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
N = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
SIGMAS = [float(s) for s in (sys.argv[3].split(",") if len(sys.argv) > 3
                             else "8,4,2,1,0.5,0.25".split(","))]
cosg = float(np.cos(np.deg2rad(GAM/60.0)))
gnode, vnode = fold_nodes()
ana = float(np.interp(GAM, gnode, vnode))
# nearest exact fold node
i = int(np.argmin(np.abs(gnode - GAM)))
print(f"gamma={GAM}' N={N}")
print(f"fold interpolated at {GAM}': {ana:.7e}")
print(f"nearest fold NODE {gnode[i]:.5f}': {vnode[i]:.7e}")

grid = ex.build_grid(n_lambda=N)
print(f"dlam={grid.dlam:.4f} Mpc")
zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
ref_true = ex.fk_reference(grid, zt)[0, 3]
print(f"ref_true (tabulated zeta, discrete local fold): {ref_true:.7e}  /fold={ref_true/ana:.6f}")

rows = []
for sig in SIGMAS:
    t0 = time.time()
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate="smeared")
    M = ex.calib_matrix(grid, V, A, rho, sig, "smeared")
    zi = ex.zeta_injected(M, Q)
    ref_inj = ex.fk_reference(grid, zi)[0, 3]
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    exact = ((T1+T2) + (T1+T2).T)[0, 3]
    t1 = T1[0, 3] + T1[3, 0]
    fid = np.abs(zi - 0.5*0 - zi).max()  # placeholder
    # injection fidelity vs symmetrised tabulated zeta
    zs = (zt + zt.transpose(0,2,1,3) + zt.transpose(0,3,1,2)
          + zt.transpose(0,1,3,2) + zt.transpose(0,2,3,1) + zt.transpose(0,3,2,1))/6.0
    fid = np.abs(zi - zs).max()/np.abs(zs).max()
    rows.append((sig, rho, exact, ref_inj, fid, t1))
    print(f"sigma={sig:7.3f} rho={rho:.6e} exact={exact:.7e} "
          f"exact/fold={exact/ana:.7f} exact/ref_inj={exact/ref_inj:.7f} "
          f"ref_inj/fold={ref_inj/ana:.6f} injfid={fid:.3e} T1share={t1/exact:.4f} "
          f"[{time.time()-t0:.0f}s]", flush=True)

d = dict(zip([r[0] for r in rows], [r[2]/ana for r in rows]))
if 8.0 in d and 4.0 in d:
    print(f"\npublished 2-pt linear extrapolation 2*r(4)-r(8) = {2*d[4.0]-d[8.0]:.6f}")
if 4.0 in d and 2.0 in d:
    print(f"            2*r(2)-r(4) = {2*d[2.0]-d[4.0]:.6f}")
if 2.0 in d and 1.0 in d:
    print(f"            2*r(1)-r(2) = {2*d[1.0]-d[2.0]:.6f}")
smallest = min(d)
print(f"DIRECT value at smallest sigma={smallest}: {d[smallest]:.6f}")
np.savez("review/adjudication/_v1_g%.4g_N%d.npz" % (GAM, N),
         sig=[r[0] for r in rows], exact=[r[2] for r in rows],
         ref_inj=[r[3] for r in rows], fid=[r[4] for r in rows],
         ana=ana, ref_true=ref_true, node_gamma=gnode[i], node_val=vnode[i])
