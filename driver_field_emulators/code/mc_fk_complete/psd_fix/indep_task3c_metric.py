"""Can the brief's 'max 42.5% / 0.20%' be reproduced?  And how many windows
does the stock builder actually drive negative?"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, assemble6, _bg, _ds
import sigma2_repaired as R

bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
ev = IndepDensity(background=bgd, apply_c0=False)
stock = _ds.Sigma2Builder(background=bgd, apply_c0=False)
r0 = R.R0Reference(background=bgd, apply_c0=False)
r1 = R.R1Knots(background=bgd, apply_c0=False)
G = [1.0, 17.3]

def run(LAM, tag, ref_builder):
    out = {}
    for name, b in (("stock", stock), ("R1", r1)):
        es = []
        for gam in G:
            cgx = float(np.cos(np.deg2rad(gam/60.0)))
            for cc in (1.0, cgx):
                r = np.array([ref_builder.matrix(cc, float(l))[0, 0] for l in LAM])
                v = np.array([b.matrix(cc, float(l))[0, 0] for l in LAM])
                es.append(np.abs(v/r - 1.0))
        e = np.concatenate(es)
        out[name] = (e.max(), np.median(e))
    print(f"  {tag:38s} stock max {out['stock'][0]:8.2%} med {out['stock'][1]:.2e} | "
          f"R1 max {out['R1'][0]:7.3%} med {out['R1'][1]:.2e}")

print("=== grids, reference = R0Reference (as the brief used) ===")
run(np.linspace(410, 2310, 600), "linspace(410,2310,600)", r0)
run(np.linspace(410, 2310, 300), "linspace(410,2310,300)", r0)
run(np.linspace(406, 2328, 600), "linspace(406,2328,600)", r0)
run(np.linspace(406, 2313.029, 600), "linspace(406,2313.029,600)", r0)
run(np.linspace(420, 2300, 600), "linspace(420,2300,600)", r0)
run(np.geomspace(410, 2310, 600), "geomspace(410,2310,600)", r0)
print("=== same grids, reference = independent evaluator ===")
run(np.linspace(410, 2310, 600), "linspace(410,2310,600)", ev)
run(np.linspace(406, 2328, 600), "linspace(406,2328,600)", ev)
run(np.geomspace(410, 2310, 600), "geomspace(410,2310,600)", ev)

print("\n=== how many contiguous windows does the STOCK builder drive NEGATIVE? ===")
FINE = np.arange(410.0, 2310.0 + 1e-9, 0.25)
for gam in G:
    cgx = float(np.cos(np.deg2rad(gam/60.0)))
    s00 = np.array([stock.matrix(1.0, float(l))[0, 0] for l in FINE])
    mn6 = np.array([np.linalg.eigvalsh(np.asarray(
        __import__("sachs_mc_core")._assemble_6x6(stock, cgx, float(l)))).min()
        for l in FINE])
    for arr, lab in ((s00, "Sigma00(within) < 0"), (mn6, "min eig 6x6 < 0")):
        bad = arr < 0
        idx = np.flatnonzero(bad); wins = []
        if idx.size:
            s = p = idx[0]
            for k in idx[1:]:
                if k == p+1: p = k; continue
                wins.append((s, p)); s = p = k
            wins.append((s, p))
        seg = [f"[{FINE[a]:.1f},{FINE[b]:.1f}]" for a, b in wins]
        print(f"  gamma={gam:5.2f}'  {lab:22s}: {len(wins):2d} windows, "
              f"{bad.sum()*0.25:6.2f} Mpc total")
        print(f"      {' '.join(seg)}")
    # and the independent reference for contrast
    mn6r = np.array([np.linalg.eigvalsh(assemble6(ev, cgx, float(l))).min() for l in FINE])
    print(f"      independent evaluator on the same fine grid: "
          f"{(mn6r < 0).sum()} of {FINE.size} non-PSD, min eig {mn6r.min():+.3e}")
