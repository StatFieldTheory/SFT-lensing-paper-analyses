"""How big is the per-shell change in Sigma2, not just in its integral?

The MC consumes Sigma2(cos, lam) shell by shell, so the pointwise change matters
even where the integral telescopes.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import sigma2_repaired as R  # noqa: E402
from t12_diag_aware_builder import R4DiagAware  # noqa: E402

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
lam = np.linspace(410.0, 2320.0, 400)
r1 = R.R1Knots(background=bg, apply_c0=False)
r4 = R4DiagAware(background=bg, apply_c0=False)
cur = ds.Sigma2Builder(background=bg, apply_c0=False)
for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    a = np.array([cur.matrix(cosg, float(t))[0, 0] for t in lam])
    b = np.array([r1.matrix(cosg, float(t))[0, 0] for t in lam])
    c = np.array([r4.matrix(cosg, float(t))[0, 0] for t in lam])
    r = c / b - 1.0
    print(f"gamma={gam}':  Sigma2_00 R4 vs R1  median {np.median(r):+.3f}  "
          f"min {r.min():+.3f}  max {r.max():+.3f}   "
          f"| R1 vs current median {np.median(b/a-1):+.3f}")
    q = np.percentile(np.abs(r), [50, 90, 99])
    print(f"          |R4/R1 - 1| percentiles 50/90/99: "
          f"{q[0]:.3f} / {q[1]:.3f} / {q[2]:.3f}")
