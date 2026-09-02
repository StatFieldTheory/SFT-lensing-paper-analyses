"""FF Monte-Carlo: lambda-step convergence and the anchor test.

Usage:  run_nl_scan.py <gamma_arcmin> <tag> [nl,nl,...] [--anchor]

The FF vertex enters the Euler recursion as a DRIFT, so its discretisation error
is O(dlam); the noise term is exact by construction (``resp`` telescopes).  This
sweeps n_lambda at fixed statistics and fits the continuum limit.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _bootstrap
ds, bgm, mc = _bootstrap.wire()
from ff_anti import ff_moment_anti

GAM = float(sys.argv[1])
TAG = sys.argv[2]
NLS = [int(x) for x in sys.argv[3].split(",")]
ANCHOR = "--anchor" in sys.argv
NSEED, NR = 8, 96_000

cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
ana = float(np.interp(GAM, np.asarray(d["g"], float), np.asarray(d["ff_mom"], float)))

print(f"gamma = {GAM}'  anchor={'ON' if ANCHOR else 'OFF'}  "
      f"analytic FF = {ana:.6e}   ({NSEED} seeds x {NR:,} real)")
print(f"{'n_lambda':>9} {'dlam[Mpc]':>10} {'FF_mc':>13} {'seedSEM':>10} "
      f"{'intSE':>10} {'/analytic':>10} {'+/-':>8} {'sec':>6}")
rows = []
for nl in NLS:
    t0 = time.time(); vals = []; ses = []
    for s in range(NSEED):
        cfg = mc.MCConfig(n_real=NR, batch_size=3000, n_lambda=nl,
                          use_f_vertex=True, apply_anchor=ANCHOR, seed=11 + 7 * s)
        v, se, npair, nb = ff_moment_anti(mc, cfg, GAM)
        vals.append(v); ses.append(se)
    v = np.array(vals); m = v.mean()
    sem = v.std(ddof=1) / np.sqrt(v.size)
    ise = float(np.mean(ses)) / np.sqrt(NSEED)
    dl = (bgm.LAM_SOURCE_BASELINE - 406.0) / (nl - 1)
    print(f"{nl:>9} {dl:>10.4f} {m:>13.6e} {sem:>10.2e} {ise:>10.2e} "
          f"{m/ana:>10.4f} {sem/ana:>8.4f} {time.time()-t0:>6.0f}", flush=True)
    rows.append((nl, dl, m, sem, ise, np.array(vals)))

np.savez(Path(__file__).resolve().parent / "out" / f"nl_{TAG}.npz",
         gamma=GAM, anchor=ANCHOR, analytic=ana,
         nl=np.array([r[0] for r in rows]), dlam=np.array([r[1] for r in rows]),
         ff=np.array([r[2] for r in rows]), sem=np.array([r[3] for r in rows]),
         intse=np.array([r[4] for r in rows]),
         per_seed=np.array([r[5] for r in rows]), n_real=NSEED * NR)
