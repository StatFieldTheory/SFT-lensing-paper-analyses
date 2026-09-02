"""Why does R0Reference integrate differently from np.gradient of the same F?"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import sigma2_repaired as R  # noqa: E402
import corr_op_C_callable as _corr_op  # noqa: E402

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
cosg = float(np.cos(np.deg2rad(1.0 / 60.0)))
n2 = ds._n2_of_cos(cosg)
r0 = R.R0Reference(background=bg, apply_c0=False)
cur = ds.Sigma2Builder(background=bg, apply_c0=False)

la = np.linspace(1150.0, 1450.0, 61)
D = np.asarray(bg.D(la), float)
C = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0] for t in la])
F = D ** 4 * C
Fp = np.gradient(F, la, edge_order=2)
print(f"{'lam':>9} {'R0 sigma2':>13} {'grad(F)/D^4':>13} {'current':>13} "
      f"{'R0/grad':>9}")
for i, t in enumerate(la):
    a = r0.matrix(cosg, float(t))[0, 0]
    b = Fp[i] / D[i] ** 4
    c = cur.matrix(cosg, float(t))[0, 0]
    print(f"{t:9.2f} {a:13.5e} {b:13.5e} {c:13.5e} {a/b:9.3f}")
