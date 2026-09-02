"""Validate the antithetic estimator against the production simulate_ff_crn."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import _bootstrap
ds, bgm, mc = _bootstrap.wire()
from ff_anti import ff_moment_anti

print(f"{'gamma':>6} {'nl':>5} {'n_real':>8} | {'production':>11} {'SE':>9} | "
      f"{'antithetic':>11} {'SE':>9} | {'var ratio':>9} {'t_prod':>7} {'t_anti':>7}")
for gam in (1.0, 6.7):
    for nl in (1000, 4000):
        nr = 48000
        cfg = mc.MCConfig(n_real=nr, batch_size=3000, n_lambda=nl,
                          use_f_vertex=True, apply_anchor=False, seed=11)
        t0 = time.time(); r = mc.simulate_ff_crn(cfg, gam); t1 = time.time()
        cfg2 = mc.MCConfig(n_real=nr, batch_size=3000, n_lambda=nl,
                           use_f_vertex=True, apply_anchor=False, seed=11)
        a, ase, npair, nb = ff_moment_anti(mc, cfg2, gam); t2 = time.time()
        # variance ratio per unit CPU: production SE^2 * t_prod vs anti SE^2 * t_anti
        vr = (r.ff_moment_err[0, 0] ** 2 * (t1 - t0)) / (ase ** 2 * (t2 - t1))
        print(f"{gam:>6.2f} {nl:>5d} {nr:>8d} | {r.ff_moment[0,0]:>11.4e} "
              f"{r.ff_moment_err[0,0]:>9.2e} | {a:>11.4e} {ase:>9.2e} | "
              f"{vr:>9.1f} {t1-t0:>7.1f} {t2-t1:>7.1f}  (nblow={nb})", flush=True)
