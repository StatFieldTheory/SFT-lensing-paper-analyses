"""PSD on the production grid + a dense scan, and the Order-0 anchor under
both the d6 protocol (global GL-256) and a node-panelised converged rule."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import sigma2_repaired as R
import candidates as K

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
L = K.TABLE_LAM
MK = [
    ("stock",      lambda ac: ds.Sigma2Builder(background=bg, apply_c0=ac)),
    ("R0 ref",     lambda ac: R.R0Reference(background=bg, apply_c0=ac)),
    ("R1 knots",   lambda ac: R.R1Knots(background=bg, apply_c0=ac)),
    ("R2 density", lambda ac: R.R2Density(background=bg, apply_c0=ac)),
    ("R2 FIXED",   lambda ac: K.R2Fixed(background=bg, apply_c0=ac)),
    ("R7 akima",   lambda ac: K.R7Akima(background=bg, apply_c0=ac)),
    ("R8 hermite", lambda ac: K.R8HermiteNodes(background=bg, apply_c0=ac)),
    ("R5 EXACT",   lambda ac: K.R5Exact(background=bg, apply_c0=ac)),
]

# ---- PSD -------------------------------------------------------------------
print("=== PSD of the 6x6.  grid A = production N=1000 uniform on [406,2313.029];")
print("    grid B = 4000 uniform (0.48 Mpc); grid C = 400 points inside [1300,1320] ===")
gA = np.linspace(406.0, bgm.LAM_SOURCE_BASELINE, 1000)
gB = np.linspace(406.0, 2328.0, 4000)
gC = np.linspace(1300.5, 1320.0, 400)
GAMS = [1.0, 17.3]
print(f"{'scheme':<11} " + " ".join(f"{'A g=%.1f'%g:>9}" for g in GAMS)
      + " " + " ".join(f"{'B g=%.1f'%g:>9}" for g in GAMS)
      + " " + " ".join(f"{'C g=%.1f'%g:>9}" for g in GAMS)
      + f" {'worst min/max eig':>19}")
for name, mk in MK:
    b = mk(False)
    row = []; worst = np.inf
    for grid in (gA, gB, gC):
        for g in GAMS:
            c = float(np.cos(np.deg2rad(g/60.0)))
            nv = 0
            for x in grid:
                M = core._assemble_6x6(b, c, float(x))
                w = np.linalg.eigvalsh(0.5*(M+M.T))
                if w[0] < 0: nv += 1
                worst = min(worst, w[0]/max(abs(w[-1]), 1e-300))
            row.append(nv)
    print(f"{name:<11} " + " ".join(f"{v:>9d}" for v in row) + f" {worst:>19.3e}")

# ---- Order-0 anchor --------------------------------------------------------
cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
GAM = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3, 293.9, 596.9])

def order0_panelled(b, cos_gamma, lam_f, n_per_panel=48):
    """Same integral as ds.order0_mc but the quadrature is SPLIT at the corr_op
    table nodes, so the discontinuities of the true density fall on panel edges."""
    brk = np.concatenate([[b.lam_lo], L[(L > b.lam_lo+1e-9) & (L < lam_f-1e-9)], [lam_f]])
    nodes, wts = leggauss(int(n_per_panel))
    tot = 0.0
    for lo, hi in zip(brk[:-1], brk[1:]):
        la = 0.5*(hi-lo)*nodes + 0.5*(hi+lo)
        jw = 0.5*(hi-lo)*wts
        s = np.array([b.matrix(cos_gamma, float(l))[0,0] for l in la])
        G = ds.order0_window(b.bg, la, lam_f, n_gauss=96)
        tot += float(np.sum(jw*s*G**2))
    return tot

for rule in ("GL-256 (d6 protocol)", "node-panelled GL-48 (converged)"):
    print(f"\n=== Order-0 anchor, median order0_mc / O0_analysis3 -- {rule} ===")
    print(f"{'scheme':<11} " + " ".join(f"{'g=%.1f'%g:>9}" for g in GAM[:6])
          + f" {'MEDIAN(10)':>11} {'vs stock':>10}")
    med = {}
    for name, mk in MK:
        b = mk(False)
        rat = []
        for gam in GAM:
            cg = float(np.cos(np.deg2rad(gam/60.0)))
            a3 = float(np.interp(gam, g_a3, o0_a3))
            v = (ds.order0_mc(cg, lam_f=bgm.LAM_SOURCE_BASELINE, builder=b, n_gauss=256)
                 if rule.startswith("GL-256")
                 else order0_panelled(b, cg, bgm.LAM_SOURCE_BASELINE))
            rat.append(float(v)/a3)
        med[name] = float(np.median(rat))
        print(f"{name:<11} " + " ".join(f"{v:>9.5f}" for v in rat[:6])
              + f" {med[name]:>11.5f} {med[name]/med['stock']-1:>+9.4%}")

# panel-count convergence for the exact builder
print("\n=== convergence of the panelled rule (R5 EXACT, gamma=17.3') ===")
b = K.R5Exact(background=bg, apply_c0=False)
cg = float(np.cos(np.deg2rad(17.3/60.0)))
prev = None
for n in (12, 24, 48, 96):
    v = order0_panelled(b, cg, bgm.LAM_SOURCE_BASELINE, n_per_panel=n)
    print("   n_per_panel=%3d  %.12e%s" % (n, v, "" if prev is None else "   rel diff %.2e"%abs(v/prev-1)))
    prev = v
for n in (128, 256, 512, 1024, 2048):
    v = ds.order0_mc(cg, lam_f=bgm.LAM_SOURCE_BASELINE, builder=b, n_gauss=n)
    print("   GLOBAL GL n=%4d  %.12e   vs panelled %+.3e" % (n, v, v/prev-1))
