"""V9: the proposed order0_mc patch, verbatim, incl. window n_gauss."""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT="/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0,ROOT); sys.path.insert(0,ROOT+"/psd_fix")
import _bootstrap
ds,bgm,core=_bootstrap.wire()
import sigma2_repaired as R
import corr_op_C_callable as _corr_op
_TABLE_LAM=np.load(_corr_op.TABLE_PATH,allow_pickle=True)["lambda_grid"].astype(float)

def order0_mc_patched(cos_gamma, lam_f=None, builder=None, n_gauss=48):
    b = builder if builder is not None else ds._builder()
    lf = float(lam_f) if lam_f is not None else bgm.LAM_SOURCE_BASELINE
    edges = _TABLE_LAM[(_TABLE_LAM > b.lam_lo + 1e-9) & (_TABLE_LAM < lf - 1e-9)]
    edges = np.concatenate([[b.lam_lo], edges, [lf]])
    nodes, weights = leggauss(int(n_gauss))
    la = np.concatenate([0.5*(hi-lo)*nodes + 0.5*(hi+lo)
                         for lo, hi in zip(edges[:-1], edges[1:])])
    jw = np.concatenate([0.5*(hi-lo)*weights for lo, hi in zip(edges[:-1], edges[1:])])
    sig00 = np.array([b.matrix(cos_gamma, float(l))[0, 0] for l in la])
    G = ds.order0_window(b.bg, la, lf, n_gauss=n_gauss)
    return float(np.sum(jw * sig00 * G ** 2))

bg=bgm.Background(lam_min=120.0,lam_max=bgm.LAM_SOURCE_BASELINE+5.0)
d=np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
               "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"),allow_pickle=True)
g_a3,o0_a3=np.asarray(d["g"],float),np.asarray(d["o0"],float)
GAM=np.array([0.5,1.0,2.0,5.0,10.0,17.3,44.4,114.3,293.9,596.9])
B={"stock":ds.Sigma2Builder(background=bg,apply_c0=False),
   "R1":R.R1Knots(background=bg,apply_c0=False)}
print(f"{'builder':>7} {'n_gauss':>8} {'median ratio':>13}")
for k,b in B.items():
    for ng in (24,48,96):
        m=[order0_mc_patched(float(np.cos(np.deg2rad(g/60.0))),lam_f=bgm.LAM_SOURCE_BASELINE,
                             builder=b,n_gauss=ng)/float(np.interp(g,g_a3,o0_a3)) for g in GAM]
        print(f"{k:>7} {ng:>8} {np.median(m):>13.6f}",flush=True)
    m48=[ds.order0_mc(float(np.cos(np.deg2rad(g/60.0))),lam_f=bgm.LAM_SOURCE_BASELINE,
                      builder=b,n_gauss=48)/float(np.interp(g,g_a3,o0_a3)) for g in GAM]
    print(f"{k:>7} {'CURRENT48':>8} {np.median(m48):>13.6f}")
