"""ADJUDICATOR check 2: decompose exact/fold into its sigma-dependent and
sigma-INdependent factors, and follow the sigma-independent factor to dlam->0.

  exact/fold = (exact/ref_inj) x (ref_inj/fold)

(exact/ref_inj) -> 1 as sigma->0 (verified in v1 to 7e-6).
(ref_inj/fold)  is sigma-independent while the injection is exact; it is the
                discrete local fold of the SYMMETRISED tabulated cumulant
                against the paper's analytic fold.

Also tests: is ref_inj really the symmetrised reference (the leg-asymmetry
claim), and how much of ref_true/fold = 1.0084 is discretisation.
"""
import sys, math
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

def sym6(z):
    return (z + z.transpose(0,1,3,2) + z.transpose(0,2,1,3) + z.transpose(0,2,3,1)
            + z.transpose(0,3,1,2) + z.transpose(0,3,2,1))/6.0

gnode, vnode = fold_nodes()
GAMS = [float(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [1.0]
NS = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [500,1000,2000,4000]

# (0) confirm ref_inj == symmetrised reference at a working sigma, N=1000
g = 1.0; cosg = float(np.cos(np.deg2rad(g/60.0)))
grid = ex.build_grid(n_lambda=1000)
zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
V, A, Q, rho = ex.node_stats(grid, cosg, 8.0, calibrate="smeared")
M = ex.calib_matrix(grid, V, A, rho, 8.0, "smeared")
zi = ex.zeta_injected(M, Q)
r_inj = ex.fk_reference(grid, zi)[0,3]
r_sym = ex.fk_reference(grid, sym6(zt))[0,3]
r_raw = ex.fk_reference(grid, zt)[0,3]
print("CHECK 0 (N=1000, gamma=1', sigma=8):")
print(f"  ref_inj  = {r_inj:.8e}")
print(f"  ref_sym  = {r_sym:.8e}   ref_inj/ref_sym = {r_inj/r_sym:.9f}")
print(f"  ref_raw  = {r_raw:.8e}   ref_raw/ref_sym = {r_raw/r_sym:.6f}  <- leg asymmetry")
print()

print(f"{'gamma':>9} {'N':>6} {'dlam':>8} {'ref_sym':>13} {'ref_raw':>13} "
      f"{'sym/fold':>9} {'raw/fold':>9}")
out = {}
for g in GAMS:
    cosg = float(np.cos(np.deg2rad(g/60.0)))
    i = int(np.argmin(np.abs(gnode - g)))
    # use the EXACT fold node value when g is a node, else interpolate
    at_node = abs(gnode[i]-g) < 1e-6
    fold = vnode[i] if at_node else float(np.interp(g, gnode, vnode))
    for N in NS:
        gr = ex.build_grid(n_lambda=N)
        z = np.array([ex._DS.zeta6(cosg, float(l)) for l in gr.lam])
        rs = ex.fk_reference(gr, sym6(z))[0,3]
        rr = ex.fk_reference(gr, z)[0,3]
        out[(g,N)] = (rs, rr, fold)
        print(f"{g:>9.5f} {N:>6d} {gr.dlam:>8.4f} {rs:>13.6e} {rr:>13.6e} "
              f"{rs/fold:>9.5f} {rr/fold:>9.5f}", flush=True)
    # Richardson in 1/N on the two finest
    a, b = out[(g,NS[-2])][0], out[(g,NS[-1])][0]
    rich = b + (b-a)*(NS[-2]/(NS[-1]-NS[-2]))*1.0
    # proper 1/N Richardson: f(N)=f_inf + c/N  ->  f_inf=(N2*f2-N1*f1)/(N2-N1)
    f_inf = (NS[-1]*b - NS[-2]*a)/(NS[-1]-NS[-2])
    print(f"   -> 1/N Richardson ref_sym = {f_inf:.6e}   /fold = {f_inf/out[(g,NS[-1])][2]:.5f}")
    ar, br = out[(g,NS[-2])][1], out[(g,NS[-1])][1]
    f_infr = (NS[-1]*br - NS[-2]*ar)/(NS[-1]-NS[-2])
    print(f"   -> 1/N Richardson ref_raw = {f_infr:.6e}   /fold = {f_infr/out[(g,NS[-1])][2]:.5f}")
