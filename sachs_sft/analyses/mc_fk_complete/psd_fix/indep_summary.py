"""One-command reproduction of every headline number in the independent audit."""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, assemble6, _bg, _ds
import sigma2_repaired as R

bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
ev = IndepDensity(background=bgd, apply_c0=False)
stock = _ds.Sigma2Builder(background=bgd, apply_c0=False)
r1 = R.R1Knots(background=bgd, apply_c0=False)
r0 = R.R0Reference(background=bgd, apply_c0=False)

print("### 1. PSD of the TRUE density, [410,2310], 200 pts/cell + both node limits")
pts = []
edges = np.unique(np.concatenate(([410.0], NODES[(NODES > 410) & (NODES < 2310)], [2310.0])))
for lo, hi in zip(edges[:-1], edges[1:]):
    pts += [(float(x), "auto") for x in np.linspace(lo, hi, 202)[1:-1]]
for x in NODES[(NODES > 410) & (NODES < 2310)]:
    pts += [(float(x), "left"), (float(x), "right")]
for gam in (0.5, 1.0, 5.0, 17.3, 60.0):
    cg = float(np.cos(np.deg2rad(gam/60.0)))
    mn = np.array([np.linalg.eigvalsh(assemble6(ev, cg, l, side=s)).min() for l, s in pts])
    print(f"   gamma={gam:5.2f}'  non-PSD {int((mn<0).sum()):d}/{len(pts)}   "
          f"min eig {mn.min():+.4e}")

print("\n### 2/3. relative error in Sigma00 vs the independent evaluator, 600 lambda")
LAM = np.linspace(410.0, 2310.0, 600)
for gam in (1.0, 17.3):
    cg = float(np.cos(np.deg2rad(gam/60.0)))
    for blk, cc in (("within", 1.0), ("cross", cg)):
        ref = np.array([ev.matrix(cc, float(l))[0, 0] for l in LAM])
        line = f"   gamma={gam:5.2f}' {blk:>6}: "
        for nm, b in (("stock", stock), ("R1", r1), ("R0", r0)):
            e = np.abs(np.array([b.matrix(cc, float(l))[0, 0] for l in LAM])/ref - 1)
            line += f"{nm} max {e.max():9.4%} med {np.median(e):.2e}  "
        print(line)

print("\n### 4. Order-0 anchor: median(order0/O0_a3) over the 10 paper gammas")
anch = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                    "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"),
               allow_pickle=True)
GA3, O0A3 = np.asarray(anch["g"], float), np.asarray(anch["o0"], float)
GAM = [0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3, 293.9, 596.9]

PANEL_EDGES = [406.0] + [float(x) for x in NODES if 406.0 < x < LAMF] + [LAMF]
_nod, _wt = leggauss(48)
_L, _W = [], []
for _a, _c in zip(PANEL_EDGES[:-1], PANEL_EDGES[1:]):
    _L.append(0.5 * (_c - _a) * _nod + 0.5 * (_c + _a))
    _W.append(0.5 * (_c - _a) * _wt)
PANEL_LA = np.concatenate(_L)
PANEL_W = np.concatenate(_W)
PANEL_G2 = _ds.order0_window(bgd, PANEL_LA, LAMF, n_gauss=64) ** 2


def anchor_gl256(bld):
    """Exactly the d6_anchor.py protocol: one global 256-point GL rule."""
    return float(np.median([
        float(_ds.order0_mc(float(np.cos(np.deg2rad(gv / 60.0))), lam_f=LAMF,
                            builder=bld, n_gauss=256))
        / float(np.interp(gv, GA3, O0A3)) for gv in GAM]))


def anchor_panel(bld):
    """Jump-aware rule: 48-point GL on each of the 20 inter-node panels."""
    out = []
    for gv in GAM:
        cgv = float(np.cos(np.deg2rad(gv / 60.0)))
        sv = np.array([bld.matrix(cgv, float(x))[0, 0] for x in PANEL_LA])
        out.append(float(np.sum(PANEL_W * sv * PANEL_G2))
                   / float(np.interp(gv, GA3, O0A3)))
    return float(np.median(out))


print(f"   {'builder':>12} {'GL256 (d6 protocol)':>21} {'panel48 (jump-aware)':>22}")
anc = {}
for nm2, bld2 in (("stock", stock), ("R1Knots", r1), ("indep", ev)):
    anc[nm2] = (anchor_gl256(bld2), anchor_panel(bld2))
    print(f"   {nm2:>12} {anc[nm2][0]:>21.5f} {anc[nm2][1]:>22.5f}")
print(f"   {'R1/stock':>12} {anc['R1Knots'][0]/anc['stock'][0]-1:>+21.4%} "
      f"{anc['R1Knots'][1]/anc['stock'][1]-1:>+22.4%}")
print(f"   {'indep/stock':>12} {anc['indep'][0]/anc['stock'][0]-1:>+21.4%} "
      f"{anc['indep'][1]/anc['stock'][1]-1:>+22.4%}")
