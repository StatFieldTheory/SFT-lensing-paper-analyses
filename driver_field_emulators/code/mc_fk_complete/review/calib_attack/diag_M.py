"""Do the alternative calibration matrices actually tend to Sigma2 as sigma->0?

Prints, per calibration, the kernel-weighted mean and max of
|M[k]-Sigma2[k]|_max / |Sigma2[k]|_max at fixed sigma/dlam, plus the edge nodes.
Cheap: no expectation, no Q.
"""
import sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/review/calib_attack")
import numpy as np, altcal
import fk_expect_exact as ex

gam = float(sys.argv[1]); pairs = sys.argv[2]; names = sys.argv[3].split(",")
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
for item in pairs.split(","):
    Ns, sigs = item.split(":"); N, sig = int(Ns), float(sigs)
    grid = ex.build_grid(n_lambda=N)
    V, A, Z, rho = altcal.base_stats(grid, cosg, sig)
    S2 = V * (2 * sig)
    kern = grid.dlam * grid.Wd * grid.Hd; kern = kern / kern.sum()
    s2n = np.abs(S2).reshape(N, -1).max(1)
    print(f"\nN={N} sigma={sig} dlam={grid.dlam:.4f} s/dl={sig/grid.dlam:.3f}"
          f"  gamma={gam}'")
    print(f"{'calib':>12} {'kwmean|dM|/S2':>14} {'max|dM|/S2':>11} "
          f"{'M/S2 k=0':>9} {'k=2':>7} {'k=mid':>7} {'k=N-1':>7}")
    for nm in names:
        M = altcal.calib_matrices(grid, V, A, rho, sig, nm)
        d = np.abs(M - S2).reshape(N, -1).max(1) / s2n
        r = np.abs(M).reshape(N, -1).max(1) / s2n
        print(f"{nm:>12} {float((d*kern).sum()):>14.5f} {float(d.max()):>11.4f} "
              f"{r[0]:>9.4f} {r[2]:>7.4f} {r[N//2]:>7.4f} {r[-1]:>7.4f}")
