"""A3: the attack's own tests.
  (1) vary R1's per-cell sampling density over a wide range
  (2) vary R0's FD step over a wide range
  (3) vary the stock builder's n_nodes (incl. the 'already PSD' 80 and 120)
Each measured under the d6 rule (global GL-256) AND a jump-aware panelised rule.
If the Order-0 shift is a real property of the tabulated C it must be
independent of (1) and (2)."""
import numpy as np, sys, time
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix/adv")
import advlib as A
ds = A.ds; bg = A.BG

GL   = lambda lo,hi: A.rule_gl(lo,hi,256)
PAN  = lambda lo,hi: A.rule_panel(lo,hi,48)
PAN96= lambda lo,hi: A.rule_panel(lo,hi,96)

E = A.Exact(background=bg, apply_c0=False)
med = lambda b, r: float(np.median(A.anchor(b, r)))
mE_gl, mE_pan = med(E, GL), med(E, PAN)
b_cur = ds.Sigma2Builder(background=bg, apply_c0=False)
mC_gl, mC_pan = med(b_cur, GL), med(b_cur, PAN)
print(f"EXACT   GL256 {mE_gl:.6f}   PANEL48 {mE_pan:.6f}   PANEL96 {med(E,PAN96):.6f}")
print(f"stock160 GL256 {mC_gl:.6f}   PANEL48 {mC_pan:.6f}   PANEL96 {med(b_cur,PAN96):.6f}")
print()

print("(1) R1 per-cell sampling spacing sweep  [production dx = 12 Mpc]")
print(f"{'dx[Mpc]':>9} {'#samples':>9} {'GL256 med':>11} {'shift vs stockGL':>17} "
      f"{'PANEL med':>11} {'shift vs stockPAN':>18} {'maxErr vs EXACT':>16}")
lamq = np.linspace(410, 2310, 601)
cos1 = float(np.cos(np.deg2rad(1.0/60)))
ex = np.array([E.matrix(cos1, float(l))[0,0] for l in lamq])
for dx in (48.0, 24.0, 12.0, 6.0, 3.0, 1.5, 0.75):
    b = A.R1Var(background=bg, apply_c0=False, dx=dx)
    b._density_splines(cos1)
    ns = sum(len(sp[0][0].x) for _,_,sp in b._cache[round(cos1,12)])
    v = np.array([b.matrix(cos1, float(l))[0,0] for l in lamq])
    e1, e2 = med(b, GL), med(b, PAN)
    print(f"{dx:9.2f} {ns:9d} {e1:11.6f} {e1/mC_gl-1:+16.4%} {e2:11.6f} "
          f"{e2/mC_pan-1:+17.4%} {np.max(np.abs(v/ex-1)):16.3e}")
print()

print("(2) R0 finite-difference step sweep  [default fd_h = 0.5 Mpc]")
print(f"{'fd_h[Mpc]':>10} {'GL256 med':>11} {'shift vs stockGL':>17} {'PANEL med':>11} {'shift vs stockPAN':>18}")
for h in (8.0, 4.0, 2.0, 1.0, 0.5, 0.2, 0.05, 0.01):
    b = A.R.R0Reference(background=bg, apply_c0=False, fd_h=h)
    e1, e2 = med(b, GL), med(b, PAN)
    print(f"{h:10.3f} {e1:11.6f} {e1/mC_gl-1:+16.4%} {e2:11.6f} {e2/mC_pan-1:+17.4%}")
print()

print("(3) stock builder n_nodes sweep  [production n_nodes = 160]")
print(f"{'n_nodes':>8} {'dlam':>7} {'GL256 med':>11} {'shift vs 160GL':>15} "
      f"{'PANEL med':>11} {'vs EXACT(PAN)':>14} {'nonPSD/600':>11} {'minEig/maxEig':>14}")
for n in (60, 80, 120, 160, 240, 320, 480, 640, 1280, 2560):
    b = ds.Sigma2Builder(background=bg, apply_c0=False, n_nodes=n)
    e1, e2 = med(b, GL), med(b, PAN)
    # PSD scan on 600 lambda at gamma=1'
    bad = 0; worst = 1e9
    for l in lamq:
        w, c = ds.sigma2_blocks_from(b, cos1, float(l)) if hasattr(ds,'sigma2_blocks_from') else (None,None)
        W = b.matrix(1.0, float(l)); X = b.matrix(cos1, float(l))
        M = np.zeros((6,6)); M[:3,:3]=W; M[3:,3:]=W; M[:3,3:]=X; M[3:,:3]=X.T
        ev = np.linalg.eigvalsh(M)
        if ev.min() < 0: bad += 1
        worst = min(worst, ev.min()/max(abs(ev.max()),1e-300))
    print(f"{n:8d} {(2328-406)/(n-1):7.2f} {e1:11.6f} {e1/mC_gl-1:+14.4%} {e2:11.6f} "
          f"{e2/mE_pan-1:+13.4%} {bad:11d} {worst:14.3e}")
