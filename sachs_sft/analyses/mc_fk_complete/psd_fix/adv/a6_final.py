"""A6: (a) do the coarse 'already PSD' stock builders get the density right?
(b) the converged anchor, and how it compares to the hardcoded ANCHOR_RATIO."""
import numpy as np, sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix/adv")
import advlib as A
ds = A.ds; bg = A.BG
E = A.Exact(background=bg, apply_c0=False)
cos1 = float(np.cos(np.deg2rad(1.0/60)))
lamq = np.linspace(410, 2310, 601)
ex = np.array([E.matrix(cos1,float(l))[0,0] for l in lamq])
GL = lambda lo,hi: A.rule_gl(lo,hi,256); PAN = lambda lo,hi: A.rule_panel(lo,hi,96)
med = lambda b,r: float(np.median(A.anchor(b,r)))
mE = med(E, PAN)
print("(a) does coarse stock sampling get the DENSITY right?  gamma=1', vs EXACT")
print(f"{'builder':>14} {'nonPSD/601':>11} {'medErr':>10} {'maxErr':>10} {'n(>10%)':>8} "
      f"{'GL256 anchor':>13} {'PANEL anchor':>13} {'PANEL vs EXACT':>15}")
rows = [("stock n=60",ds.Sigma2Builder(background=bg,apply_c0=False,n_nodes=60)),
        ("stock n=80",ds.Sigma2Builder(background=bg,apply_c0=False,n_nodes=80)),
        ("stock n=120",ds.Sigma2Builder(background=bg,apply_c0=False,n_nodes=120)),
        ("stock n=160*",ds.Sigma2Builder(background=bg,apply_c0=False)),
        ("stock n=640",ds.Sigma2Builder(background=bg,apply_c0=False,n_nodes=640)),
        ("R1 knots",A.R.R1Knots(background=bg,apply_c0=False)),
        ("R1 dx=1.5",A.R1Var(background=bg,apply_c0=False,dx=1.5)),
        ("R0 h=0.01",A.R.R0Reference(background=bg,apply_c0=False,fd_h=0.01)),
        ("EXACT",E)]
for nm,b in rows:
    v = np.array([b.matrix(cos1,float(l))[0,0] for l in lamq]); rel=np.abs(v/ex-1)
    bad=0
    for l in lamq:
        W=b.matrix(1.0,float(l)); X=b.matrix(cos1,float(l))
        M=np.zeros((6,6)); M[:3,:3]=W; M[3:,3:]=W; M[:3,3:]=X; M[3:,:3]=X.T
        if np.linalg.eigvalsh(M).min()<0: bad+=1
    g,p = med(b,GL), med(b,PAN)
    print(f"{nm:>14} {bad:11d} {np.median(rel):10.3e} {np.max(rel):10.3e} "
          f"{int((rel>0.1).sum()):8d} {g:13.6f} {p:13.6f} {p/mE-1:+14.4%}")
print(f"\n  hardcoded driver_stats.ANCHOR_RATIO = {ds.ANCHOR_RATIO}")
print(f"  converged (PANEL96) stock = {med(ds.Sigma2Builder(background=bg,apply_c0=False),PAN):.6f}"
      f"   EXACT = {mE:.6f}   -> repair moves it by {mE/med(ds.Sigma2Builder(background=bg,apply_c0=False),PAN)-1:+.4%}"
      f",  vs the {mE/ds.ANCHOR_RATIO-1:+.2%} the hardcoded constant is already off")
