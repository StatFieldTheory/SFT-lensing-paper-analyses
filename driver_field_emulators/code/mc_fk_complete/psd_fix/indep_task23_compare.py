"""TASK 2 + 3: R1Knots and the STOCK builder measured against the independent
cell-local evaluator, plus an independent localisation of the ringing windows."""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, assemble6, _bg, _ds
import sigma2_repaired as R          # used ONLY as the object under test

bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
ev = IndepDensity(background=bgd, apply_c0=False)
stock = _ds.Sigma2Builder(background=bgd, apply_c0=False)
r1 = R.R1Knots(background=bgd, apply_c0=False)

LAM = np.linspace(410.0, 2310.0, 600)
print("600 lambda in [410,2310]; min distance to a table node: %.4f Mpc"
      % np.min(np.abs(LAM[:, None] - NODES[None, :])))

GAM = [1.0, 17.3]
store = {}
for gam in GAM:
    cg = float(np.cos(np.deg2rad(gam / 60.0)))
    ref_w = np.array([ev.matrix(1.0, l) for l in LAM])
    ref_c = np.array([ev.matrix(cg, l) for l in LAM])
    out = {}
    for name, b in (("stock", stock), ("R1Knots", r1)):
        w = np.array([b.matrix(1.0, l) for l in LAM])
        c = np.array([b.matrix(cg, l) for l in LAM])
        out[name] = (w, c)
    store[gam] = (cg, ref_w, ref_c, out)

def relerr(a, r, mode):
    if mode == "00":
        return np.abs(a[:, 0, 0] - r[:, 0, 0]) / np.abs(r[:, 0, 0])
    # max over the 9 entries, normalised by the reference matrix max-norm
    return np.abs(a - r).reshape(len(a), -1).max(1) / np.abs(r).reshape(len(r), -1).max(1)

print("\n=== TASK 2/3: relative error vs the independent evaluator ===")
for mode, lab in (("00", "entry (0,0)"), ("mat", "max over 9 entries / ||ref||_max")):
    print(f"\n-- metric: {lab} --")
    print(f"{'gamma':>7} {'block':>7} {'builder':>9} {'max err':>11} {'median err':>12} "
          f"{'p95':>11} {'argmax lam':>11}")
    for gam in GAM:
        cg, ref_w, ref_c, out = store[gam]
        for bl, ref, idx in (("within", ref_w, 0), ("cross", ref_c, 1)):
            for name in ("stock", "R1Knots"):
                e = relerr(out[name][idx], ref, mode)
                print(f"{gam:>7.2f} {bl:>7} {name:>9} {e.max():>11.4%} "
                      f"{np.median(e):>12.3e} {np.percentile(e,95):>11.3e} "
                      f"{LAM[int(np.argmax(e))]:>11.2f}")

print("\n=== pooled over both gammas and both blocks (2400 matrices) ===")
for mode, lab in (("00", "entry (0,0)"), ("mat", "max-entry / ||ref||_max")):
    for name in ("stock", "R1Knots"):
        e = []
        for gam in GAM:
            cg, ref_w, ref_c, out = store[gam]
            e.append(relerr(out[name][0], ref_w, mode))
            e.append(relerr(out[name][1], ref_c, mode))
        e = np.concatenate(e)
        print(f"  {lab:26s} {name:>9}: max {e.max():.4%}  median {np.median(e):.3e}")

# ---------------------------------------------------------------- task 3 ----
print("\n=== TASK 3: independent localisation of the stock builder's ringing ===")
gam = 17.3
cg, ref_w, ref_c, out = store[gam]
e = relerr(out["stock"][0], ref_w, "00")          # within block, entry 00
THR = 0.05                                        # 5% relative error
bad = e > THR
# group consecutive flagged lambdas into windows
idx = np.flatnonzero(bad)
wins = []
if idx.size:
    start = idx[0]; prev = idx[0]
    for k in idx[1:]:
        if k == prev + 1:
            prev = k; continue
        wins.append((start, prev)); start = k; prev = k
    wins.append((start, prev))
print(f"  windows with >|{THR:.0%}| error in Sigma00(within), gamma={gam}': {len(wins)}")
print(f"  {'#':>3} {'lam range [Mpc]':>22} {'width':>7} {'peak err':>10} "
      f"{'peak lam':>9} {'prev node':>10} {'offset':>8}")
for n, (i0, i1) in enumerate(wins, 1):
    lo, hi = LAM[i0], LAM[i1]
    k = i0 + int(np.argmax(e[i0:i1 + 1]))
    prev_node = NODES[np.searchsorted(NODES, LAM[k], side="right") - 1]
    print(f"  {n:>3} [{lo:8.1f},{hi:8.1f}] {hi-lo:>7.1f} {e[k]:>10.2%} "
          f"{LAM[k]:>9.1f} {prev_node:>10.0f} {LAM[k]-prev_node:>8.1f}")

print("\n  signed error profile around node 1300 (stock/indep - 1), 4 Mpc steps:")
cgx = cg
for l in np.arange(1284.0, 1345.0, 4.0):
    a = stock.matrix(1.0, float(l))[0, 0]; r = ev.matrix(1.0, float(l))[0, 0]
    print(f"    lam={l:7.1f}  stock={a:.5e}  indep={r:.5e}  err={a/r-1:+8.2%}")
