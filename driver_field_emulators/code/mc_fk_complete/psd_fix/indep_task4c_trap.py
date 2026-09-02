"""Does a uniform Riemann/trapezoid rule (what the MC's own uniform lambda grid
realises) converge to the panelised value?"""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, _bg, _ds
bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
ev = IndepDensity(background=bgd, apply_c0=False)
cg = float(np.cos(np.deg2rad(17.3/60.0)))

edges = [406.0] + [float(x) for x in NODES if 406.0 < x < LAMF] + [LAMF]
nod, wt = leggauss(64); LA=[]; W=[]
for a, c in zip(edges[:-1], edges[1:]):
    LA.append(0.5*(c-a)*nod+0.5*(c+a)); W.append(0.5*(c-a)*wt)
LA=np.concatenate(LA); W=np.concatenate(W)
G2 = _ds.order0_window(bgd, LA, LAMF, n_gauss=64)**2
s = np.array([ev.matrix(cg, float(l))[0,0] for l in LA])
ref = float(np.sum(W*s*G2))
print(f"panelised (64/panel) reference: {ref:.10e}")
for n in (501, 1001, 2001, 4001, 8001, 16001):
    la = np.linspace(406.0, LAMF, n)
    w = np.full(n, (LAMF-406.0)/(n-1)); w[0]*=0.5; w[-1]*=0.5
    g2 = _ds.order0_window(bgd, la, LAMF, n_gauss=64)**2
    ss = np.array([ev.matrix(cg, float(l))[0,0] for l in la])
    v = float(np.sum(w*ss*g2))
    print(f"  trapezoid n={n:>6} (dlam={ (LAMF-406.0)/(n-1):6.3f} Mpc): "
          f"{v:.10e}  rel diff {v/ref-1:+.3e}")
