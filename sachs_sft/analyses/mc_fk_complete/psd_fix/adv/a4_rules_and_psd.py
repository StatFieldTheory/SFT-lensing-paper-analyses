"""A4: (a) is the panelised rule the right limit?  cross-check with a uniform
trapezoid, which is what the MC's own uniform lambda grid realises.
(b) density error of stock vs the parameter-free EXACT.
(c) PSD of the EXACT density in the critical window, both one-sided limits."""
import numpy as np, sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix/adv")
import advlib as A
ds = A.ds; bg = A.BG
E = A.Exact(background=bg, apply_c0=False)
b_cur = ds.Sigma2Builder(background=bg, apply_c0=False)
b_r1 = A.R.R1Knots(background=bg, apply_c0=False)
med = lambda b,r: float(np.median(A.anchor(b,r)))

print("(a) rule convergence, median anchor")
print(f"{'rule':>22} {'stock160':>11} {'R1':>11} {'EXACT':>11} {'R1/stock-1':>12} {'EXACT/stock-1':>14}")
rules = [("GL 128",        lambda lo,hi: A.rule_gl(lo,hi,128)),
         ("GL 256 (d6)",   lambda lo,hi: A.rule_gl(lo,hi,256)),
         ("GL 512",        lambda lo,hi: A.rule_gl(lo,hi,512)),
         ("GL 2048",       lambda lo,hi: A.rule_gl(lo,hi,2048)),
         ("PANEL GL24",    lambda lo,hi: A.rule_panel(lo,hi,24)),
         ("PANEL GL48",    lambda lo,hi: A.rule_panel(lo,hi,48)),
         ("PANEL GL96",    lambda lo,hi: A.rule_panel(lo,hi,96)),
         ("TRAP n=501",    lambda lo,hi: A.rule_trap(lo,hi,501)),
         ("TRAP n=1001",   lambda lo,hi: A.rule_trap(lo,hi,1001)),
         ("TRAP n=4001",   lambda lo,hi: A.rule_trap(lo,hi,4001)),
         ("TRAP n=16001",  lambda lo,hi: A.rule_trap(lo,hi,16001))]
for nm, r in rules:
    a,b,c = med(b_cur,r), med(b_r1,r), med(E,r)
    print(f"{nm:>22} {a:11.6f} {b:11.6f} {c:11.6f} {b/a-1:+11.4%} {c/a-1:+13.4%}")

print("\n(b) density error vs EXACT, gamma=1', 601 lambda in [410,2310] and a 0.25 Mpc scan")
cos1 = float(np.cos(np.deg2rad(1.0/60)))
for tag, lam in (("601 log-free grid", np.linspace(410,2310,601)),
                 ("7601 @ 0.25 Mpc",   np.arange(410,2310.01,0.25))):
    ex = np.array([E.matrix(cos1,float(l))[0,0] for l in lam])
    for nm,b in (("stock160",b_cur),("R1",b_r1)):
        v = np.array([b.matrix(cos1,float(l))[0,0] for l in lam])
        rel = np.abs(v/ex-1)
        print(f"   {tag:18s} {nm:9s} median {np.median(rel):.3e}  max {np.max(rel):.3e}"
              f"  at lam={lam[np.argmax(rel)]:.2f}   n(>10%) {int((rel>0.1).sum())}")

print("\n(c) PSD of the EXACT (parameter-free) density, both one-sided limits at nodes")
for g in (0.5, 1.0, 5.0, 17.3):
    cosg = float(np.cos(np.deg2rad(g/60)))
    pts = list(np.arange(410, 2310.01, 0.5))
    bad = 0; mn = 1e9
    for l in pts:
        for side in ("left","right"):
            W = E.matrix(1.0, float(l), side=side); X = E.matrix(cosg, float(l), side=side)
            M = np.zeros((6,6)); M[:3,:3]=W; M[3:,3:]=W; M[:3,3:]=X; M[3:,:3]=X.T
            ev = np.linalg.eigvalsh(M)
            if ev.min()<0: bad += 1
            mn = min(mn, ev.min()/max(abs(ev.max()),1e-300))
    print(f"   gamma={g:6.2f}'  nonPSD {bad}/{2*len(pts)}   min(minEig/maxEig) = {mn:+.3e}")

print("\n(d) the critical node lam=1308.93: stock vs EXACT eigenvalues (gamma=1')")
for l in (1307.0, 1308.93, 1310.84, 1312.0):
    for nm,b in (("stock160",b_cur),("R1",b_r1),("EXACT",E)):
        W=b.matrix(1.0,l); X=b.matrix(cos1,l)
        M=np.zeros((6,6)); M[:3,:3]=W; M[3:,3:]=W; M[:3,3:]=X; M[3:,:3]=X.T
        ev=np.linalg.eigvalsh(M)
        print(f"   lam={l:8.2f} {nm:9s} Sigma00={b.matrix(1.0,l)[0,0]:+.5e}  "
              f"minEig={ev.min():+.4e}  minEig/maxEig={ev.min()/abs(ev.max()):+.3e}")
    print()
