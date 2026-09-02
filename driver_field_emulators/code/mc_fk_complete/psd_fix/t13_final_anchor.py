"""Final anchor table: GL-256 (what order0_mc does) vs a converged trapezoid."""
from __future__ import annotations
import sys
import numpy as np
from pathlib import Path
ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import sigma2_repaired as R  # noqa: E402
from t12_diag_aware_builder import R4DiagAware  # noqa: E402

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
lam_f = bgm.LAM_SOURCE_BASELINE
d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"),
            allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
la = np.linspace(406.01, lam_f, 4001)
G = ds.order0_window(bg, la, lam_f, n_gauss=256)
B = {"current": ds.Sigma2Builder(background=bg, apply_c0=False),
     "R1 knots": R.R1Knots(background=bg, apply_c0=False),
     "R4 diag-aware": R4DiagAware(background=bg, apply_c0=False)}
GAM = [0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3]
gl = {k: [] for k in B}
tz = {k: [] for k in B}
for gam in GAM:
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    a3 = float(np.interp(gam, g_a3, o0_a3))
    for k, b in B.items():
        gl[k].append(float(ds.order0_mc(cosg, lam_f=lam_f, builder=b, n_gauss=256)) / a3)
        s = np.array([b.matrix(cosg, float(t))[0, 0] for t in la])
        tz[k].append(float(np.trapezoid(s * G ** 2, la)) / a3)
print("median order0_mc / O0_analysis3 over gamma in [0.5', 114.3']")
print(f"{'builder':>15} {'GL-256 (order0_mc)':>20} {'trapezoid n=4001':>19}")
for k in B:
    print(f"{k:>15} {np.median(gl[k]):>20.5f} {np.median(tz[k]):>19.5f}")
print()
for k in B:
    print(f"{k:>15}  shift vs current:  GL-256 {np.median(gl[k])/np.median(gl['current'])-1:+8.3%}"
          f"   trapezoid {np.median(tz[k])/np.median(tz['current'])-1:+8.3%}")
