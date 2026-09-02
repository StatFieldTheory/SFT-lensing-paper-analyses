"""How badly does the HARD-CODED Gauss-Legendre rule in order0_mc fail on the
repaired density?  n_gauss=48 is order0_mc's default (used by driver_stats'
own __main__ ANCHOR self-check); n_gauss=256 is what sachs_mc_core.simulate
uses for the `o0_analytic` it reports next to every MC number.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R
from numpy.polynomial.legendre import leggauss

LAM_F = bgm.LAM_SOURCE_BASELINE; LAM_LO = 406.0
bg = bgm.Background(lam_min=120.0, lam_max=LAM_F + 5.0)
BREAKS = R._cells(LAM_LO, LAM_F)

def panel_o0(b, cosg, npt=64):
    nd, wt = leggauss(npt); tot = 0.0
    for lo, hi in zip(BREAKS[:-1], BREAKS[1:]):
        la = 0.5*(hi-lo)*nd + 0.5*(hi+lo); jw = 0.5*(hi-lo)*wt
        s = np.array([b.matrix(cosg, float(l))[0,0] for l in la])
        G = ds.order0_window(bg, la, LAM_F, n_gauss=64)
        tot += float(np.sum(jw*s*G**2))
    return tot

NG = [48, 96, 128, 256, 512, 1024, 2048, 4096]
GAM = [0.5, 1.0, 5.0, 17.3, 44.4, 114.3]
for name, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
    b = cls(background=bg, apply_c0=False)
    print(f"\n=== {name}: order0_mc(n_gauss) / panelised truth ===")
    print(f"{'gam':>7} {'truth':>12} " + " ".join(f"{f'n={n}':>9}" for n in NG))
    for gam in GAM:
        c = float(np.cos(np.deg2rad(gam/60.0)))
        tru = panel_o0(b, c)
        vals = [ds.order0_mc(c, lam_f=LAM_F, builder=b, n_gauss=n)/tru for n in NG]
        print(f"{gam:>7.1f} {tru:>12.5e} " + " ".join(f"{v:>9.5f}" for v in vals))

# --- driver_stats' OWN __main__ anchor self-check, replayed with each builder --
print("\n=== driver_stats.__main__ [Part A] ANCHOR self-check "
      "(apply_c0=True, n_gauss=48 default) ===")
import math
anc = np.load("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
              "/sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz",
              allow_pickle=True)
m = (anc["a"]==0)&(anc["b"]==0)&(anc["order"]==0)
xs, ys, vs = anc["x"][m], anc["y"][m], anc["value"][m]
gam, cos = [], []
for xi, yi in zip(xs, ys):
    xi=np.asarray(xi,float); yi=np.asarray(yi,float)
    cc=float(np.clip(np.dot(xi/np.linalg.norm(xi), yi/np.linalg.norm(yi)),-1,1))
    gam.append(math.degrees(math.acos(cc))*60.0); cos.append(cc)
o=np.argsort(gam); gam=np.asarray(gam)[o]; cos=np.asarray(cos)[o]; vs=np.asarray(vs,float)[o]
tf = float(anc["t_final"][0])
for name, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
    b = cls(background=bg, apply_c0=True)
    o0 = np.array([ds.order0_mc(c, lam_f=tf, builder=b) for c in cos])
    ratio = o0/vs; pos = (vs>0)&(gam<200.0)
    med, sd = float(np.median(ratio[pos])), float(np.std(ratio[pos]))
    ok = abs(med-1.0)<0.10 and sd<0.05
    print(f"  {name:>6}: median={med:.4f} std={sd:.4f} "
          f"range=[{ratio[pos].min():.4f},{ratio[pos].max():.4f}] "
          f"-> ANCHOR {'PASS' if ok else 'FAIL'}   (n={pos.sum()} gammas, t_final={tf:.3f})")
    b0 = cls(background=bg, apply_c0=False)
    o0r = np.array([ds.order0_mc(c, lam_f=tf, builder=b0) for c in cos])
    print(f"          raw median (would-be new ANCHOR_RATIO) = {np.median((o0r/vs)[pos]):.5f}")
