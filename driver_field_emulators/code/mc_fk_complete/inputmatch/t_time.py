"""Timing + internal-SE vs seed-scatter calibration for simulate_ff_crn."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
import _bootstrap
ds, bgm, mc = _bootstrap.wire()

for (nl, nr) in ((500, 6000), (1000, 6000), (2000, 6000), (4000, 6000)):
    t0 = time.time()
    cfg = mc.MCConfig(n_real=nr, batch_size=3000, n_lambda=nl,
                      use_f_vertex=True, apply_anchor=False, seed=11)
    r = mc.simulate_ff_crn(cfg, 1.0)
    dt = time.time() - t0
    v, se = r.ff_moment[0, 0], r.ff_moment_err[0, 0]
    print(f"nl={nl:5d} nr={nr:6d}  ff={v:.5e} intSE={se:.2e} ({se/abs(v)*100:.2f}%)"
          f"  t={dt:.1f}s  -> {dt/(nr*nl)*1e6:.4f} us per (real*node)", flush=True)
