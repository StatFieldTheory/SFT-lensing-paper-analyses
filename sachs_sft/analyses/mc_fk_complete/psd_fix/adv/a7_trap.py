"""A7: what shift does the MC's OWN rule realise?  sachs_mc_core uses
lam = linspace(lam_min, lam_source, n_lambda) with trapezoid endpoint weights."""
import numpy as np, sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix/adv")
import advlib as A
ds = A.ds; bg = A.BG
E = A.Exact(background=bg, apply_c0=False)
cur = ds.Sigma2Builder(background=bg, apply_c0=False)
r1 = A.R.R1Knots(background=bg, apply_c0=False)
med = lambda b,r: float(np.median(A.anchor(b,r)))
truth = med(E, lambda lo,hi: A.rule_panel(lo,hi,96))
print(f"converged truth (EXACT, panelised GL96) = {truth:.6f}")
print(f"{'n_lambda':>9} {'dlam':>7} {'stock':>10} {'R1':>10} {'R1/stock-1':>12} "
      f"{'stock/truth-1':>14} {'R1/truth-1':>12}")
sh=[]
for n in (400, 600, 800, 1000, 1200, 1500, 2000, 3000):
    r = lambda lo,hi,n=n: A.rule_trap(lo,hi,n)
    a,b = med(cur,r), med(r1,r)
    sh.append(b/a-1)
    print(f"{n:9d} {(A.LAM_F-406)/(n-1):7.2f} {a:10.6f} {b:10.6f} {b/a-1:+11.4%} "
          f"{a/truth-1:+13.4%} {b/truth-1:+11.4%}")
sh=np.array(sh)
print(f"\nMC-rule realised shift R1/stock-1: mean {sh.mean():+.4%}  rms {sh.std():.4%}  "
      f"range [{sh.min():+.4%}, {sh.max():+.4%}]   -- vs the claimed +1.61%")
