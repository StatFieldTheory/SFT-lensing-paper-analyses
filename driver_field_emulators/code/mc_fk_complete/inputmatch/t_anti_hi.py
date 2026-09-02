"""High-statistics agreement test: production estimator vs antithetic.

Pushes the PRODUCTION simulate_ff_crn to enough realisations that its own error
is a few tenths of a percent of the analytic FF, so the antithetic wrapper is
tested for BIAS, not merely consistency.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _bootstrap
ds, bgm, mc = _bootstrap.wire()
from ff_anti import ff_moment_anti

GAM, NL, NSEED, NR = 1.0, 1000, 8, 300_000
prod, anti = [], []
pse, ase_ = [], []
t0 = time.time()
for s in range(NSEED):
    sd = 101 + 13 * s
    c = mc.MCConfig(n_real=NR, batch_size=3000, n_lambda=NL,
                    use_f_vertex=True, apply_anchor=False, seed=sd)
    r = mc.simulate_ff_crn(c, GAM)
    prod.append(r.ff_moment[0, 0]); pse.append(r.ff_moment_err[0, 0])
    a, ase, npair, nb = ff_moment_anti(mc, c, GAM)
    anti.append(a); ase_.append(ase)
    print(f"  seed {sd}: prod {prod[-1]:.5e} +/-{pse[-1]:.1e} | "
          f"anti {a:.5e} +/-{ase:.1e}", flush=True)
p = np.array(prod); a = np.array(anti)
pm, am = p.mean(), a.mean()
pe = p.std(ddof=1) / np.sqrt(p.size); ae = a.std(ddof=1) / np.sqrt(a.size)
print(f"\nn_real total = {NSEED*NR:,}  ({time.time()-t0:.0f} s)")
print(f"production  {pm:.6e} +/- {pe:.2e}  (seed SEM)   internal SE/sqrt(n) "
      f"{np.mean(pse)/np.sqrt(NSEED):.2e}")
print(f"antithetic  {am:.6e} +/- {ae:.2e}  (seed SEM)   internal SE/sqrt(n) "
      f"{np.mean(ase_)/np.sqrt(NSEED):.2e}")
d = pm - am
print(f"difference  {d:+.3e}  = {d/am*100:+.2f}%  of antithetic; "
      f"combined sigma {np.hypot(pe, ae):.2e} -> {d/np.hypot(pe, ae):+.2f} sigma")
np.savez(Path(__file__).resolve().parent / "out" / "t_anti_hi.npz",
         prod=p, anti=a, pse=np.array(pse), ase=np.array(ase_))
