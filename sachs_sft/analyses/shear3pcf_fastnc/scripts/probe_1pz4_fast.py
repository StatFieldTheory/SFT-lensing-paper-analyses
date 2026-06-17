#!/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
"""FAST surplus-(A) check: the per-shell (1+z)^4 radial-convention ratio.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python (canoes-free
analytic algebra + reads the PyCCL oracle npz).

Pure analytic power-count of the deployed-canoes 3-pt radial response vs the PyCCL
standard, per shell, decomposed into (1+z) powers; plus the lensing-weighted
<(1+z)^4>(gamma) cross-check vs the deployed/PyCCL ~10x oracle ratio.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/zzhang/projects/canoes")
_HERE = Path(__file__).resolve().parent
_OUT = _HERE.parent / "outputs"
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402

OMEGA_M = 0.3160919980475834
H = 0.6711
Z_SOURCE = 5.0
H0 = 100.0 / 299792.458  # h/Mpc, canoes _H0_H_PER_MPC


def main() -> int:
    zf = np.linspace(0.0, Z_SOURCE, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=OMEGA_M, h=H) for zi in zf])
    af = 1.0 / (1.0 + zf)
    from scipy.integrate import cumulative_trapezoid
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(Z_SOURCE / 40, Z_SOURCE, 40)
    chi_phys = np.interp(z, zf, chif)
    chi_s = float(chif[-1])
    a = 1.0 / (1.0 + z)
    opz = 1.0 + z

    K_geo = chi_phys * (chi_s - chi_phys) / chi_s
    A_pois = -1.5 * OMEGA_M * H0 ** 2 * opz
    # deployed canoes 3-pt radial response (the factor multiplying the shared
    # D^4/chi^4 angular machinery), folded in chi: a^8 K_geo^3 (A(1+z)^4)^3
    resp_can = (a ** 8) * (K_geo ** 3) * (A_pois * opz ** 4) ** 3
    # standard (PyCCL) radial response: g^3, g = 3/2 Om H0^2 (1+z) K_geo
    g_std = 1.5 * OMEGA_M * H0 ** 2 * opz * K_geo
    resp_std = g_std ** 3
    ratio_shell = resp_can / resp_std
    sym = -(opz ** 4)   # symbolic prediction

    print("=" * 84)
    print("SURPLUS (A): per-shell deployed-canoes / standard radial response ratio")
    print("=" * 84)
    print(f"{'shell':>5s} {'z':>7s} {'1+z':>7s} {'S_can/S_std':>13s} {'-(1+z)^4':>10s} "
          f"{'ratio/(-(1+z)^4)':>16s}")
    for i in (5, 10, 20, 30, 39):
        print(f"{i:5d} {z[i]:7.3f} {opz[i]:7.3f} {ratio_shell[i]:13.5e} {sym[i]:10.3f} "
              f"{ratio_shell[i]/sym[i]:16.6f}")
    print(f"\n  symbolic: a^8 [A(a)(1+z)^4]^3 / [3/2 Om H0^2 (1+z) K]^3")
    print(f"          = a^8 * (-1)^3 (1+z)^12 = a^8 (-(1+z)^12) = -(1+z)^4")
    print(f"  max |ratio_shell/(-(1+z)^4) - 1| = {np.max(np.abs(ratio_shell/sym - 1.0)):.3e}")

    # cross-check vs deployed/PyCCL oracle ratio (~10x) via lensing weight
    print()
    print("=" * 84)
    print("CROSS-CHECK: lensing-weighted (1+z)^4 vs deployed-canoes/PyCCL oracle ratio")
    print("=" * 84)
    orc = np.load(_OUT / "scalar_vs_pyccl_tree.npz", allow_pickle=True)
    ga = orc["gamma_arcmin"]
    realcan_std = orc["ratio_realcan_std"]
    z_eff = orc["z_eff"]
    pred = (1.0 + z_eff) ** 4
    print(f"{'gamma':>8s} {'realcan/std':>11s} {'z_eff':>7s} {'(1+ze)^4':>9s} "
          f"{'ratio/(1+ze)^4':>15s}")
    for i in range(ga.size):
        print(f"{ga[i]:8.3f} {realcan_std[i]:11.3f} {z_eff[i]:7.3f} {pred[i]:9.3f} "
              f"{realcan_std[i]/pred[i]:15.4f}")

    np.savez(_OUT / "powercount_1pz4_fast.npz",
             z=z, opz=opz, ratio_shell=ratio_shell, sym=sym,
             gamma_arcmin=ga, realcan_std=realcan_std, z_eff=z_eff, pred_1pz4=pred)
    print(f"\nsaved -> {_OUT/'powercount_1pz4_fast.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
