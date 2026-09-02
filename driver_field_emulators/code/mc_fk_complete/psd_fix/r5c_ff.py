"""ATTACK (c): does the FF channel move when the repaired builder is swapped in?

sachs_mc_core._precompute builds its OWN `_ds.Sigma2Builder`, so the builder is
swapped by binding the class onto the driver_stats MODULE OBJECT at runtime.
The original files are untouched.  A counting wrapper on `_cholesky_psd` records
how often the eigen-floor fires in each configuration.
"""
import sys, time
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R
from dataclasses import replace

_STOCK = ds.Sigma2Builder
_orig_chol = core._cholesky_psd
COUNT = {"n": 0}
def _counting(M):
    try:
        np.linalg.cholesky(M)
    except np.linalg.LinAlgError:
        COUNT["n"] += 1
    return _orig_chol(M)
core._cholesky_psd = _counting

GAMS = [1.0, 10.0]
SEEDS = [12345, 20260828, 314159, 777001, 555, 90210, 424242, 1000003, 11, 22, 33, 44, 55, 66, 77, 88, 101, 202, 303, 404, 505, 606, 707, 808, 909, 1111, 2222, 3333, 4444, 5555, 6666, 7777]
NL = int(sys.argv[1]) if len(sys.argv) > 1 else 600
NR = int(sys.argv[2]) if len(sys.argv) > 2 else 12000
res = {}
for bname, cls in (("stock", _STOCK), ("R1", R.R1Knots)):
    ds.Sigma2Builder = cls
    for gam in GAMS:
        for sd in SEEDS:
            cfg = core.MCConfig(n_real=NR, batch_size=2000, n_lambda=NL,
                                use_f_vertex=True, apply_anchor=False, seed=sd)
            COUNT["n"] = 0
            t0 = time.time()
            r = core.simulate_ff_crn(cfg, gam)
            res[(bname, gam, sd)] = (float(r.ff_moment[0,0]), float(r.ff_conn[0,0]),
                                     float(r.ff_disc[0,0]), float(r.ff_moment_err[0,0]),
                                     COUNT["n"], r.n_good)
            print(f"{bname:>6} gam={gam:>5.1f} seed={sd:<9} "
                  f"ffmom={r.ff_moment[0,0]:+.6e} conn={r.ff_conn[0,0]:+.6e} "
                  f"disc={r.ff_disc[0,0]:+.6e} se={r.ff_moment_err[0,0]:.2e} "
                  f"chol_fallback={COUNT['n']} ({time.time()-t0:.0f}s)", flush=True)
ds.Sigma2Builder = _STOCK
np.savez(f"/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/"
         f"mc_fk_complete/psd_fix/_r5c_ff_N{NL}.npz",
         keys=np.array([f"{a}|{b}|{c}" for a,b,c in res], dtype=object),
         vals=np.array(list(res.values()), dtype=float))

print(f"\n=== FF channel, paired by seed (n_lambda={NL}, n_real={NR}) ===")
for gam in GAMS:
    for j, lbl in ((0,"ff_mom"), (1,"ff_conn"), (2,"ff_disc")):
        s = np.array([res[("stock",gam,sd)][j] for sd in SEEDS])
        r = np.array([res[("R1",gam,sd)][j] for sd in SEEDS])
        d = r - s
        print(f"gam={gam:>5.1f} {lbl:>8}: stock {s.mean():+.6e} +/- {s.std(ddof=1)/np.sqrt(len(SEEDS)):.2e} | "
              f"R1 {r.mean():+.6e} +/- {r.std(ddof=1)/np.sqrt(len(SEEDS)):.2e} | "
              f"paired diff {d.mean():+.4e} +/- {d.std(ddof=1)/np.sqrt(len(SEEDS)):.2e}  "
              f"({d.mean()/abs(s.mean())*100:+.4f}% of stock)")
