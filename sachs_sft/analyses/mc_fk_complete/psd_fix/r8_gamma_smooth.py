"""ATTACK (d): is the repaired Sigma2 SMOOTH in gamma, and does the xi_+/xi_-/
kappa-gamma structure survive?  Fine 60-point log gamma scan, 0.1' -> 600'.
Reports (i) pointwise Sigma2_ab at fixed lambda, (ii) the converged (panelised)
LOS integral in every component, (iii) a discrete smoothness statistic: the
second difference in log gamma, normalised by the local value.
"""
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
import sigma2_repaired as R

LAM_F = bgm.LAM_SOURCE_BASELINE; LAM_LO = 406.0
bg = bgm.Background(lam_min=120.0, lam_max=LAM_F + 5.0)
BREAKS = R._cells(LAM_LO, LAM_F)
ND, WT = leggauss(48)
_LA, _JW = [], []
for lo, hi in zip(BREAKS[:-1], BREAKS[1:]):
    _LA.append(0.5*(hi-lo)*ND + 0.5*(hi+lo)); _JW.append(0.5*(hi-lo)*WT)
_LA = np.concatenate(_LA); _JW = np.concatenate(_JW)
def G(la):
    nd, wt = leggauss(64); Dla = np.asarray(bg.D(la), float); out = np.empty_like(la)
    for i, l in enumerate(la):
        t = 0.5*(LAM_F-l)*nd + 0.5*(LAM_F+l); jw = 0.5*(LAM_F-l)*wt
        out[i] = float(np.sum(jw*(Dla[i]/np.asarray(bg.D(t)))**2))
    return out
_G2 = G(_LA)**2
LAMP = [500.0, 1000.0, 1308.9, 1400.0, 2100.7, 2280.0]   # incl. the defect window
GAM = np.geomspace(0.1, 600.0, 60)

res = {}
for bname, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
    b = cls(background=bg, apply_c0=False)
    I = np.empty((GAM.size, 3, 3)); P = np.empty((GAM.size, len(LAMP), 3, 3))
    for i, g in enumerate(GAM):
        c = float(np.cos(np.deg2rad(g/60.0)))
        S = np.array([b.matrix(c, float(l)) for l in _LA])
        I[i] = np.einsum("j,jab->ab", _JW*_G2, S)
        P[i] = np.array([b.matrix(c, float(l)) for l in LAMP])
    res[bname] = (I, P)
    print(f"[{bname}] done", flush=True)
np.savez("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/"
         "mc_fk_complete/psd_fix/_r8.npz", GAM=GAM, LAMP=np.array(LAMP),
         **{f"{k}_{n}": v for k, (I, P) in res.items()
            for n, v in (("I", I), ("P", P))})

def rough(y):
    """max |second difference in log gamma| / local |y|, interior points."""
    ly = np.asarray(y, float)
    d2 = ly[2:] - 2*ly[1:-1] + ly[:-2]
    den = np.maximum(np.abs(ly[1:-1]), 1e-300)
    return float(np.max(np.abs(d2)/den))

print("\n=== gamma-SMOOTHNESS: max |2nd difference in log-gamma| / |value| ===")
print("(60 log-spaced gamma, 0.1' -> 600'; a builder that jitters in gamma "
      "shows a large number here)")
print(f"{'quantity':>28} {'stock':>10} {'R1':>10}")
for nm, f in (("O0 kk (LOS integral)", lambda I: I[:,0,0]),
              ("O0 xi+ = 11+22",       lambda I: I[:,1,1]+I[:,2,2]),
              ("O0 xi- = 11-22",       lambda I: I[:,1,1]-I[:,2,2]),
              ("O0 kappa-gamma = 01",  lambda I: I[:,0,1])):
    print(f"{nm:>28} {rough(f(res['stock'][0])):>10.2e} {rough(f(res['R1'][0])):>10.2e}")
for j, l in enumerate(LAMP):
    for ab, nm in (((0,0),"Sigma2_00"), ((1,1),"Sigma2_11")):
        print(f"{nm+f' at lam={l:.1f}':>28} "
              f"{rough(res['stock'][1][:,j,ab[0],ab[1]]):>10.2e} "
              f"{rough(res['R1'][1][:,j,ab[0],ab[1]]):>10.2e}")

print("\n=== small-gamma limit: cross block -> within block ? ===")
Iw = {}
for bname, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
    b = cls(background=bg, apply_c0=False)
    S = np.array([b.matrix(1.0, float(l)) for l in _LA])
    Iw[bname] = np.einsum("j,jab->ab", _JW*_G2, S)
for i in (0, 5, 10):
    print(f"  gam={GAM[i]:6.3f}': stock cross/within kk = "
          f"{res['stock'][0][i,0,0]/Iw['stock'][0,0]:.6f}   R1 = "
          f"{res['R1'][0][i,0,0]/Iw['R1'][0,0]:.6f}")

print("\n=== R1/stock on the CONVERGED LOS integral, fine gamma scan ===")
I_s, I_r = res["stock"][0], res["R1"][0]
for nm, f in (("kk", lambda I: I[:,0,0]), ("xi+", lambda I: I[:,1,1]+I[:,2,2]),
              ("xi-", lambda I: I[:,1,1]-I[:,2,2]), ("kg", lambda I: I[:,0,1])):
    r = f(I_r)/f(I_s)
    ok = np.isfinite(r) & (np.abs(f(I_s)) > 1e-12*np.abs(f(I_s)).max())
    print(f"  {nm:>3}: min {r[ok].min():.6f}  max {r[ok].max():.6f}  "
          f"spread {(r[ok].max()/r[ok].min()-1)*100:+.4f}%  (n={ok.sum()}/{GAM.size})")
