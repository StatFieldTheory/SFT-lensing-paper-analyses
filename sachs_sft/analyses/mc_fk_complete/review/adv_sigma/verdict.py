"""Adversarial verdict on the sigma_lambda -> 0 extrapolation."""
import numpy as np
ANA = 1.59771599e-05
BLOCKS = [("../../_gate2_highstat.npz", "A 20260827"),
          ("../../_gate2_highstat_B.npz", "B 314159001"),
          ("../../_papergrid.npz", "C 20260828"),
          ("_mcblockD.npz", "D 777000111 (MINE)")]

print("=" * 80)
print("1.  THE LIMIT IS NOT AN EXTRAPOLATION.  Switch the regulator off on the")
print("    SAME lattice the production run uses (N = 1000, gamma = 1').")
print("=" * 80)
d = np.load("_fixN1000_g1.npz"); o = np.argsort(d["sigma"])[::-1]
dl = float(d["dlam"][0])
for i in o:
    s = d["sigma"][i]
    print(f"    sigma={s:>6.2f}  rho={np.exp(-dl/s):>9.2e}  exact/ref_inj={d['exact'][i]/d['ref_inj'][i]:.7f}"
          f"  exact/analytic={d['exact'][i]/ANA:.7f}")
TRUE = float(d["exact"][d["sigma"] == 0.25][0]) / ANA
print(f"\n    DIRECT sigma->0 value (no model): exact/analytic = {TRUE:.6f}")
print(f"    published 2-point linear extrapolation from sigma = 8,4 : 1.000514")
print(f"    method bias of the published recipe                     : "
      f"{100*(1.000514-TRUE):+.3f}%")

print()
print("=" * 80)
print("2.  MONTE-CARLO / EXACT, four independent seed blocks (112 seeds x 24000)")
print("=" * 80)
pool = {8.0: [], 4.0: []}
per = {8.0: [], 4.0: []}
EX = {}
for f, lbl in BLOCKS:
    dd = np.load(f, allow_pickle=True)
    for sig in (8.0, 4.0):
        m = (dd["gamma"] == 1.0) & (dd["sigma"] == sig)
        if not m.any():
            continue
        i = int(np.where(m)[0][0])
        s = np.asarray(dd["seeds"][i], float); E = float(dd["exact"][i]); EX[sig] = E
        se = s.std(ddof=1) / np.sqrt(s.size)
        pool[sig].append(s); per[sig].append((lbl, s.mean() / E, se / E, s.size))
        print(f"    sigma={sig:>4.1f}  {lbl:>20}  n={s.size:>3d}  MC/exact = "
              f"{s.mean()/E:.4f} +/- {se/E:.4f}")
for sig in (8.0, 4.0):
    y = np.array([p[1] for p in per[sig]]); w = 1.0 / np.array([p[2] for p in per[sig]])**2
    mu = (w*y).sum()/w.sum(); chi2 = (w*(y-mu)**2).sum(); dof = y.size-1
    s = np.concatenate(pool[sig]); se = s.std(ddof=1)/np.sqrt(s.size)
    print(f"    POOLED sigma={sig:>4.1f} n={s.size:>3d}  MC/exact = {s.mean()/EX[sig]:.5f}"
          f" +/- {se/EX[sig]:.5f}   block chi2={chi2:.2f}/{dof} -> inflate error x"
          f"{max(1.0,np.sqrt(chi2/dof)):.2f}  => +/- {se/EX[sig]*max(1.0,np.sqrt(chi2/dof)):.5f}")
n = min(len(np.concatenate(pool[8.0])), len(np.concatenate(pool[4.0])))
a8 = np.concatenate(pool[8.0])[:n]; a4 = np.concatenate(pool[4.0])[:n]
r8 = a8/EX[8.0]; r4 = a4/EX[4.0]; dif = r4-r8
print(f"\n    PAIRED test for a sigma-dependent MC bias between the two anchors:")
print(f"      (MC/exact)|_4 - (MC/exact)|_8 = {dif.mean():+.5f} +/- "
      f"{dif.std(ddof=1)/np.sqrt(n):.5f}   ({dif.mean()/(dif.std(ddof=1)/np.sqrt(n)):+.2f} sigma), n={n}")

print()
print("=" * 80)
print("3.  THE MC *DOES* DEVELOP A SYSTEMATIC AS THE REGULATOR IS REMOVED --")
print("    but only BELOW the anchors.  (MC floors the non-PSD tabulated V;")
print("    the exact evaluator does not.  Same Q, different field.)")
print("=" * 80)
print("      sigma      A neg-eig nodes    exact(MC's field)/exact(as evaluated)")
import re
for line in open("_psd.log"):
    if line.startswith("sig="):
        sg = re.search(r"sig=\s*([\d.]+)", line).group(1)
        na = re.search(r"A neg-nodes=\s*(\d+)", line).group(1)
        rt = re.search(r"ratio=([\d.]+)", line).group(1)
        print(f"      {sg:>6}   {na:>10}                    {rt:>10}   "
              f"({100*(float(rt)-1):+.3f}%)")
print("      -> +0.014% at sigma = 8 and +0.034% at sigma = 4 (the anchors),")
print("         but +2.87% at sigma = 2.  The observed MC/exact = 1.045 at")
print("         sigma = 2 is therefore mostly SYSTEMATIC, not seed scatter.")

print()
print("=" * 80)
print("4.  THE CORRECTED NUMBER")
print("=" * 80)
s8 = np.concatenate(pool[8.0]); se8 = s8.std(ddof=1)/np.sqrt(s8.size)/EX[8.0]
y = np.array([p[1] for p in per[8.0]]); w = 1.0/np.array([p[2] for p in per[8.0]])**2
mu = (w*y).sum()/w.sum(); sc = np.sqrt(((w*(y-mu)**2).sum())/(y.size-1))
r = s8.mean()/EX[8.0]
print(f"    MC/exact (sigma-independent, tested above) = {r:.5f} +/- {se8*max(1,sc):.5f}")
print(f"    exact(sigma->0)/analytic, DIRECT, N=1000   = {TRUE:.5f}")
print(f"    => MC(sigma->0)/analytic at N=1000         = {r*TRUE:.4f} +/- {se8*max(1,sc)*TRUE:.4f}")
print(f"       (paper quotes 1.002 +/- 0.006 in two places and 0.997 +/- 0.006 in one)")
dN = np.load("_refN.npz"); m = dN["sigma"] == 0.1
Ns = dN["N"][m]; ri = dN["ref_inj"][m]/ANA; x = np.argsort(Ns)
rich = 2*ri[x][-1]-ri[x][-2]
print(f"\n    grid dependence of that limit (regulator off, ref_inj/analytic):")
for i in x:
    print(f"       N={int(Ns[i]):>5d}  {ri[i]:.6f}")
print(f"       1/N Richardson (dlam -> 0)     {rich:.6f}")
print(f"    => in the JOINT limit the exact estimator sits {100*(rich-1):+.2f}% from the fold,")
print(f"       so MC/analytic would read {r*rich:.4f}, not {r*TRUE:.4f}.")
