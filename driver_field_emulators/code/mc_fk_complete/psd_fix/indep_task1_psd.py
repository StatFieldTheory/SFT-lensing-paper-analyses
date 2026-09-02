"""TASK 1: is the TRUE density PSD at every lambda in [410,2310]?

The density is piecewise smooth with jumps at the 21 table nodes, so the scan
uses a per-cell dense grid PLUS both one-sided limits at every interior node.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, assemble6, _bg

bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
ev = IndepDensity(background=bgd, apply_c0=False)

LAM_LO, LAM_HI = 410.0, 2310.0
GAM = [0.5, 1.0, 5.0, 17.3, 60.0]

# per-cell dense grid, strictly inside each cell, + node one-sided limits
pts = []   # (lam, side)
edges = np.unique(np.clip(NODES, LAM_LO, LAM_HI))
edges = edges[(edges >= LAM_LO) & (edges <= LAM_HI)]
edges = np.unique(np.concatenate(([LAM_LO], edges, [LAM_HI])))
for lo, hi in zip(edges[:-1], edges[1:]):
    if hi - lo < 1e-9:
        continue
    g = np.linspace(lo, hi, 202)[1:-1]          # strictly interior
    pts += [(float(x), "auto") for x in g]
for x in NODES:
    if LAM_LO < x < LAM_HI:
        pts.append((float(x), "left")); pts.append((float(x), "right"))
pts.append((LAM_LO, "auto")); pts.append((LAM_HI, "auto"))
print(f"scan points per gamma: {len(pts)}")

print(f"\n{'gamma':>7} {'#non-PSD':>9} {'min eig 6x6':>14} {'min eig/median':>15} "
      f"{'min Sig00 within':>17} {'min Sig00 cross':>16} {'argmin lam':>11}")
summary = {}
for gam in GAM:
    cg = float(np.cos(np.deg2rad(gam / 60.0)))
    mn = np.empty(len(pts)); s00w = np.empty(len(pts)); s00c = np.empty(len(pts))
    for i, (l, sd) in enumerate(pts):
        S6 = assemble6(ev, cg, l, side=sd)
        mn[i] = np.linalg.eigvalsh(S6).min()
        s00w[i] = ev.matrix(1.0, l, side=sd)[0, 0]
        s00c[i] = ev.matrix(cg, l, side=sd)[0, 0]
    bad = mn < 0.0
    j = int(np.argmin(mn))
    med = float(np.median(mn[~bad])) if (~bad).any() else float("nan")
    print(f"{gam:>7.2f} {int(bad.sum()):>9d} {mn.min():>14.4e} {mn.min()/med:>15.5f} "
          f"{s00w.min():>17.4e} {s00c.min():>16.4e} {pts[j][0]:>11.2f}")
    summary[gam] = (mn, bad)
    if bad.any():
        print("    NON-PSD points (first 20):")
        for k in np.flatnonzero(bad)[:20]:
            print(f"      lam={pts[k][0]:9.3f} side={pts[k][1]:5s} min eig={mn[k]:+.4e}")

# where is the minimum eigenvalue smallest relative to the local scale?
print("\n=== margin: min eigenvalue as a fraction of the largest, worst 10 lam ===")
cg = float(np.cos(np.deg2rad(17.3 / 60.0)))
rat = []
for (l, sd) in pts:
    e = np.linalg.eigvalsh(assemble6(ev, cg, l, side=sd))
    rat.append((e.min() / e.max(), l, sd))
rat.sort()
for r, l, sd in rat[:10]:
    print(f"   lam={l:9.3f} side={sd:5s}  min/max eig = {r:.4e}")
