"""TASK 4: the Order-0 anchor  median(order0_mc / O0_analysis3)  over 10 gamma,
with the independent evaluator as the builder.

Two integrations are reported:
  * the SAME protocol d6_anchor.py uses -- one global 256-point Gauss-Legendre
    rule over [406, lam_f] -- so the number is directly comparable;
  * a PANELISED rule that splits the integral at the 21 table nodes.  The true
    integrand Sigma00(la) G(la)^2 is DISCONTINUOUS at those nodes, so a single
    global GL rule converges only algebraically; the panelised value is the
    correct one.
"""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, _bg, _ds
import sigma2_repaired as R

bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
GAM = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3, 293.9, 596.9])

ev = IndepDensity(background=bgd, apply_c0=False)
builders = {
    "stock": _ds.Sigma2Builder(background=bgd, apply_c0=False),
    "R1Knots": R.R1Knots(background=bgd, apply_c0=False),
    "R0Reference": R.R0Reference(background=bgd, apply_c0=False),
    "indep": ev,
}

def order0_panel(b, cosg, n_per_panel=48):
    """int Sigma00 G^2 dla, split at the table nodes (jump-aware)."""
    edges = [406.0] + [float(x) for x in NODES if 406.0 < x < LAMF] + [LAMF]
    nod, wt = leggauss(int(n_per_panel))
    tot = 0.0
    for a, c in zip(edges[:-1], edges[1:]):
        la = 0.5*(c-a)*nod + 0.5*(c+a)
        jw = 0.5*(c-a)*wt
        s = np.array([b.matrix(cosg, float(l))[0, 0] for l in la])
        G = _ds.order0_window(bgd, la, LAMF, n_gauss=64)
        tot += float(np.sum(jw * s * G**2))
    return tot

print("=== order0_mc / O0_analysis3, raw (apply_c0=False) ===")
print(f"{'gamma':>8} {'O0 analysis-3':>14} " + " ".join(f"{k:>12}" for k in builders)
      + f" {'indep panel':>12}")
rows = {k: [] for k in builders}; rows["indep panel"] = []
for gam in GAM:
    cosg = float(np.cos(np.deg2rad(gam/60.0)))
    a3 = float(np.interp(gam, g_a3, o0_a3))
    vals = []
    for k, b in builders.items():
        v = float(_ds.order0_mc(cosg, lam_f=LAMF, builder=b, n_gauss=256))
        rows[k].append(v/a3); vals.append(v/a3)
    vp = order0_panel(ev, cosg)/a3
    rows["indep panel"].append(vp)
    print(f"{gam:>8.2f} {a3:>14.5e} " + " ".join(f"{v:>12.5f}" for v in vals)
          + f" {vp:>12.5f}")
print(f"{'median':>23} " + " ".join(f"{np.median(rows[k]):>12.5f}" for k in builders)
      + f" {np.median(rows['indep panel']):>12.5f}")
base = float(np.median(rows["stock"]))
print(f"{'shift vs stock':>23} " +
      " ".join(f"{np.median(rows[k])/base-1:>+11.4%} " for k in builders)
      + f"{np.median(rows['indep panel'])/base-1:>+11.4%}")

print("\n=== convergence of the GLOBAL Gauss-Legendre rule on the (discontinuous)")
print("    true integrand, gamma = 17.3', vs the panelised value ===")
cosg = float(np.cos(np.deg2rad(17.3/60.0)))
ref = order0_panel(ev, cosg, n_per_panel=96)
print(f"  panelised (96/panel) : {ref:.10e}")
for ng in (48, 96, 128, 192, 256, 384, 512, 1024, 2048):
    v = float(_ds.order0_mc(cosg, lam_f=LAMF, builder=ev, n_gauss=ng))
    print(f"  global GL n={ng:>5}   : {v:.10e}   rel diff {v/ref-1:+.3e}")
for npp in (24, 48, 96, 192):
    v = order0_panel(ev, cosg, n_per_panel=npp)
    print(f"  panelised n={npp:>5}   : {v:.10e}   rel diff {v/ref-1:+.3e}")
