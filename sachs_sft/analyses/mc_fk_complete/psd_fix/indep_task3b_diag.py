"""Why do my max errors disagree with the brief?  Test R0Reference itself,
validate the independent evaluator by the INTEGRAL identity, and take a fine
census of the ringing windows."""
import sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, _bg, _ds, _corr_op, _N1, _n2_of_cos
import sigma2_repaired as R

bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
ev = IndepDensity(background=bgd, apply_c0=False)
stock = _ds.Sigma2Builder(background=bgd, apply_c0=False)
r0 = R.R0Reference(background=bgd, apply_c0=False)
r1 = R.R1Knots(background=bgd, apply_c0=False)

cg = float(np.cos(np.deg2rad(17.3 / 60.0)))

# ---------------------------------------------------------------- 0 ---------
print("=== 0. INTEGRAL IDENTITY:  int_{l0}^{t} D^4 Sigma dl + D^4C|_{l0} == D^4 C|_t ===")
print("    (an end-to-end check of the independent derivative, using only C values)")
def Cd(t):
    return float(_corr_op.C_fn(_N1, float(t), _n2_of_cos(cg), float(t))[0, 0])
l0 = 410.0
g0 = float(bgd.D(l0))**4 * Cd(l0)
print(f"  {'t':>8} {'D^4 C(t) direct':>18} {'g0 + int D^4 Sigma':>20} {'rel diff':>11}")
for t in (600.0, 1000.0, 1300.0, 1301.0, 1800.0, 2200.0, 2310.0):
    # integrate cell by cell so the quadrature never straddles a jump
    edges = [l0] + [float(x) for x in NODES if l0 < x < t] + [t]
    tot = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        f = lambda l: float(bgd.D(l))**4 * ev.matrix(cg, l)[0, 0]
        tot += quad(f, a, b, limit=200, epsabs=0, epsrel=1e-11)[0]
    direct = float(bgd.D(t))**4 * Cd(t)
    print(f"  {t:>8.1f} {direct:>18.10e} {g0+tot:>20.10e} {abs((g0+tot)/direct-1):>11.2e}")

# ---------------------------------------------------------------- 1 ---------
print("\n=== 1. R0Reference measured against the independent evaluator ===")
LAM = np.linspace(410.0, 2310.0, 600)
ref = np.array([ev.matrix(cg, float(l))[0, 0] for l in LAM])
tab = {}
for name, b in (("R0Reference", r0), ("R1Knots", r1), ("stock", stock)):
    v = np.array([b.matrix(cg, float(l))[0, 0] for l in LAM])
    tab[name] = v
    e = np.abs(v/ref - 1.0)
    j = int(np.argmax(e))
    node = NODES[np.searchsorted(NODES, LAM[j], side="right")-1]
    print(f"  {name:>12}: max {e.max():.4%} at lam={LAM[j]:.2f} "
          f"(node {node:.0f} + {LAM[j]-node:.2f} Mpc), median {np.median(e):.3e}")

print("\n  -- how R0 behaves in the 10 Mpc after node 2037 (where stock peaks) --")
print(f"  {'lam':>9} {'indep':>14} {'R0':>14} {'R1':>14} {'stock':>14}")
for l in (2036.0, 2036.9, 2037.2, 2038.0, 2040.0, 2043.0, 2049.0, 2060.0):
    print(f"  {l:>9.1f} {ev.matrix(cg,l)[0,0]:>14.5e} {r0.matrix(cg,l)[0,0]:>14.5e} "
          f"{r1.matrix(cg,l)[0,0]:>14.5e} {stock.matrix(cg,l)[0,0]:>14.5e}")

# ---------------------------------------------------------------- 2 ---------
print("\n=== 2. errors EXCLUDING a guard band around every table node ===")
for guard in (0.0, 1.0, 2.0, 4.0, 8.0, 12.09, 24.0):
    m = np.min(np.abs(LAM[:, None] - NODES[None, :]), axis=1) > guard
    line = f"  guard {guard:5.2f} Mpc ({m.sum():3d}/600 kept): "
    for name in ("stock", "R1Knots", "R0Reference"):
        e = np.abs(tab[name][m]/ref[m] - 1.0)
        line += f"{name} max {e.max():8.4%} med {np.median(e):.2e}   "
    print(line)

# ---------------------------------------------------------------- 3 ---------
print("\n=== 3. fine-grid census of the stock builder's ringing windows ===")
FINE = np.arange(410.0, 2310.0 + 1e-9, 0.25)
rf = np.array([ev.matrix(1.0, float(l))[0, 0] for l in FINE])
sf = np.array([stock.matrix(1.0, float(l))[0, 0] for l in FINE])
ef = np.abs(sf/rf - 1.0)
for THR in (0.10, 0.25, 0.50, 1.00, 2.00):
    bad = ef > THR
    idx = np.flatnonzero(bad)
    wins = []
    if idx.size:
        s = p = idx[0]
        for k in idx[1:]:
            if k == p + 1:
                p = k; continue
            wins.append((s, p)); s = p = k
        wins.append((s, p))
    widths = [FINE[b]-FINE[a] for a, b in wins]
    print(f"  threshold {THR:5.0%}: {len(wins):3d} windows, "
          f"total flagged length {sum(widths):7.1f} Mpc, "
          f"widths {min(widths) if widths else 0:.2f}-{max(widths) if widths else 0:.2f} Mpc")
print()
THR = 1.00
bad = ef > THR; idx = np.flatnonzero(bad); wins = []
s = p = idx[0]
for k in idx[1:]:
    if k == p + 1:
        p = k; continue
    wins.append((s, p)); s = p = k
wins.append((s, p))
print(f"  windows at threshold {THR:.0%}  ({len(wins)} of them):")
print(f"  {'#':>3} {'start':>8} {'end':>8} {'width':>7} {'peak err':>10} "
      f"{'prev node':>10} {'start-node':>11} {'end-node':>9}")
for n, (a, b) in enumerate(wins, 1):
    lo, hi = FINE[a], FINE[b]
    k = a + int(np.argmax(ef[a:b+1]))
    nd = NODES[np.searchsorted(NODES, lo, side="right") - 1]
    print(f"  {n:>3} {lo:>8.2f} {hi:>8.2f} {hi-lo:>7.2f} {ef[k]:>10.1%} "
          f"{nd:>10.0f} {lo-nd:>11.2f} {hi-nd:>9.2f}")
print(f"\n  interior table nodes in [410,2310]: "
      f"{int(((NODES>410)&(NODES<2310)).sum())}")
