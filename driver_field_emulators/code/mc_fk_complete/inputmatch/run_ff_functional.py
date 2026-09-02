import sys, time
import numpy as np
sys.path.insert(0, "..")
from ff_functional import ff_terms, o0_from_kernel
from pathlib import Path

cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
g_a3 = np.asarray(d["g"], float)
o0_a3 = np.asarray(d["o0"], float); ff_a3 = np.asarray(d["ff_mom"], float)


def interp(g, y, x):
    return float(np.exp(np.interp(np.log(x), np.log(g), np.log(np.abs(y))))
                 * np.sign(np.interp(np.log(x), np.log(g), y)))


print("# grid convergence at gamma = 1'   (FF moment [0,0], kind=full)")
for n in (81, 121, 161, 241, 321):
    t0 = time.time()
    r = ff_terms(1.0, n=n, kind="full")
    print(f"  n={n:4d}  moment={r['moment'][0,0]:.6e}  conn22={r['conn22'][0,0]:.4e} "
          f" disc={r['disc'][0,0]:.4e}  k1k3+k3k1={(r['k1k3']+r['k3k1'])[0,0]:.4e}"
          f"   ({time.time()-t0:.0f}s)")

print("\n# O0 from the same machinery (sanity vs o0_two_ways)")
for kind in ("full", "col", "mc"):
    print(f"  {kind:5s} O0(1') = {o0_from_kernel(1.0, n=321, kind=kind):.6e}")

print("\n# FF: full vs collapsed vs mc-realised, and the analysis-3 reference")
gam = [0.5, 1.0, 2.6, 5.0, 6.7, 17.3, 44.4]
print(f"  {'gamma':>7} {'FF_full':>12} {'FF_col':>12} {'FF_mc':>12} "
      f"{'FF_a3':>12} {'full/a3':>8} {'col/full':>9} {'mc/full':>8}")
rows = []
for g in gam:
    rf = ff_terms(g, n=241, kind="full")["moment"][0, 0]
    rc = ff_terms(g, n=241, kind="col")["moment"][0, 0]
    rm = ff_terms(g, n=241, kind="mc")["moment"][0, 0]
    a3 = interp(g_a3, ff_a3, g)
    rows.append((g, rf, rc, rm, a3))
    print(f"  {g:7.2f} {rf:12.5e} {rc:12.5e} {rm:12.5e} {a3:12.5e} "
          f"{rf/a3:8.4f} {rc/rf:9.5f} {rm/rf:8.5f}")
np.savez("ff_functional.npz", rows=np.array(rows))
