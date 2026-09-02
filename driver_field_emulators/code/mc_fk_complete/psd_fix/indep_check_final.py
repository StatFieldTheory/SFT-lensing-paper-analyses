import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, assemble6, _bg, _ds
bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
ev = IndepDensity(background=bgd, apply_c0=False)
stock = _ds.Sigma2Builder(background=bgd, apply_c0=False)
cg = float(np.cos(np.deg2rad(0.5/60.0)))
def o0(): return float(_ds.order0_mc(cg, lam_f=LAMF, builder=stock, n_gauss=256))
print("stock O0(0.5') before        :", f"{o0():.10e}")
pts = []
edges = np.unique(np.concatenate(([410.0], NODES[(NODES > 410) & (NODES < 2310)], [2310.0])))
for lo, hi in zip(edges[:-1], edges[1:]):
    pts += [(float(x), "auto") for x in np.linspace(lo, hi, 202)[1:-1]]
for x in NODES[(NODES > 410) & (NODES < 2310)]:
    pts += [(float(x), "left"), (float(x), "right")]
mn = np.array([np.linalg.eigvalsh(assemble6(ev, cg, l, side=s)).min() for l, s in pts])
print("   (section-1 scan done, non-PSD", int((mn < 0).sum()), ")")
print("stock O0(0.5') after section 1:", f"{o0():.10e}")
b2 = _ds.Sigma2Builder(background=bgd, apply_c0=False)
print("FRESH stock O0(0.5')          :",
      f"{float(_ds.order0_mc(cg, lam_f=LAMF, builder=b2, n_gauss=256)):.10e}")
print("stock cache key for 0.5':", round(cg, 12), " in cache:", round(cg,12) in stock._cache)
