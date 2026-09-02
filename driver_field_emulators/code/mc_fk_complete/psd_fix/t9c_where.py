"""Locate the disagreement between the R0 density and grad(F) on one grid."""
from __future__ import annotations
import sys
import numpy as np
ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import sigma2_repaired as R  # noqa: E402
import corr_op_C_callable as _corr_op  # noqa: E402

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
lam_f = bgm.LAM_SOURCE_BASELINE
cosg = float(np.cos(np.deg2rad(1.0 / 60.0)))
n2 = ds._n2_of_cos(cosg)
la = np.linspace(406.0, lam_f, 4001)
D = np.asarray(bg.D(la), float)
G = ds.order0_window(bg, la, lam_f, n_gauss=256)
W = G ** 2 / D ** 4
C = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0] for t in la])
F = D ** 4 * C
Fp = np.gradient(F, la, edge_order=2)
r0 = R.R0Reference(background=bg, apply_c0=False)
S = np.array([r0.matrix(cosg, float(t))[0, 0] for t in la])
i_r0 = np.trapezoid(S * G ** 2, la)
i_gr = np.trapezoid(W * Fp, la)
print(f"int W*Sigma2(R0)   = {i_r0:.6e}")
print(f"int W*grad(F)      = {i_gr:.6e}   ratio {i_r0/i_gr:.4f}")
d = (S * G ** 2 - W * Fp)
cum = np.concatenate([[0.0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(la))])
print("\ncumulative excess of the R0 route (fraction of the total integral):")
for x in (500, 700, 1000, 1300, 1600, 1900, 2100, 2250, 2300, lam_f):
    j = int(np.searchsorted(la, x)) - 1
    print(f"  lam <= {x:8.1f}   {cum[j]/i_gr:+8.4f}")
k = np.argsort(np.abs(d))[::-1][:8]
print("\nlargest pointwise integrand differences:")
for j in sorted(k):
    print(f"  lam={la[j]:9.3f}  R0={S[j]:11.4e}  grad={Fp[j]/D[j]**4:11.4e}")
