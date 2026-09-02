"""PSD and accuracy of each candidate repair, against the R0 yardstick."""
import sys, time
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
lam = np.linspace(410.0, 2310.0, 600)
GAM = (1.0, 17.3)

builders = {
    "current (interpolating spline deriv)": ds.Sigma2Builder(background=bg, apply_c0=False),
    "R1 knots at table nodes": R.R1Knots(background=bg, apply_c0=False),
    "R2 density + PCHIP": R.R2Density(background=bg, apply_c0=False),
    "R3 smoothing spline": R.R3Smoothed(background=bg, apply_c0=False),
    "R0 reference (dense FD)": R.R0Reference(background=bg, apply_c0=False),
}
print("=== PSD of the assembled 6x6, and accuracy vs R0 ===")
ref = {}
b0 = builders["R0 reference (dense FD)"]
for gam in GAM:
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    ref[gam] = np.array([core._assemble_6x6(b0, cosg, float(l)) for l in lam])
for name, b in builders.items():
    t0 = time.time()
    line = [f"{name:<38}"]
    for gam in GAM:
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        M = np.array([core._assemble_6x6(b, cosg, float(l)) for l in lam])
        e = np.array([np.linalg.eigvalsh(m) for m in M])
        nneg = int((e.min(axis=1) < 0).sum())
        worst = e.min() / np.median(e[e > 0])
        rel = np.abs(M - ref[gam]).max() / np.abs(ref[gam]).max()
        med = np.median(np.abs(M - ref[gam]).reshape(len(lam), -1).max(axis=1)
                        / np.abs(ref[gam]).reshape(len(lam), -1).max(axis=1))
        line.append(f"g={gam:<5} neg={nneg:>3}/600 worst={worst:+7.4f} "
                    f"maxerr={rel:8.2e} mederr={med:8.2e}")
    print("  " + "  |  ".join(line) + f"  [{time.time()-t0:.0f}s]")
