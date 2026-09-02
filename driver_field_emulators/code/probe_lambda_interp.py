"""Which lambda interpolation of the 16-shell table is right?

The production fold reads the vertex table through the callable's
LINEAR-in-lambda interpolation; redshift_fold_trends.py uses PCHIP of
sign * ln|zeta|. On the shell grid built for z_s = 5 the two differ by
10-30% at mid-path lambda, exactly where the FK fold carries half its
weight, so the choice is a measurable systematic of the production fold.

This probe computes the TRUTH at mid-gap lambdas with the certified
bookkeeping chain (redshift_pershell.py [A2]: pref * D^4 * raw z=0
quadrature reproduces the table shells to 0.5% at gamma <= 20'), then
scores both interpolants against it. gamma <= 20' and z(lambda) <= 1.1
keep the LOW admixture below the certification tolerance.

Run (canoes venv, from driver_field_emulators/code):
    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python probe_lambda_interp.py
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator

from b_model.cosmology import FIDUCIAL, growth
from redshift_pershell import (ARCMIN, ELL_MAX, H0_HMPC, PRODUCTS, Dc,
                               collapsed_rows, raw_vertex_z0)

LAM_MID = np.array([600.0, 1040.0, 1420.0, 1700.0, 2100.0])
GAMMAS = np.array([0.5, 5.0, 20.0])


def prefactor(z: np.ndarray, chi_hmpc: np.ndarray) -> np.ndarray:
    om, h = FIDUCIAL.Omega_m, FIDUCIAL.h
    a_leg = -1.5 * om * H0_HMPC ** 2 * (1.0 + z)
    r_ttt = (a_leg * (1.0 + z) ** 4) ** 3
    jac = (1.0 + z) ** (-4)
    return jac * (2 * np.pi) ** (-4) * r_ttt * chi_hmpc ** (-4) * h ** 4


def main() -> None:
    tree = np.load(PRODUCTS / "collapsed_dense_tree_cut1000.npz")
    lam = np.asarray(tree["lambda_shells_Mpc"], float)
    rows, g_act = collapsed_rows(tree, GAMMAS)
    zeta = np.asarray(tree["zeta_TTT"], float)[rows]

    z_mid = np.asarray([float(Dc.z_of_lambda(x)) for x in LAM_MID])
    chi_mid = np.asarray([float(Dc.chi_of_lambda(x)) for x in LAM_MID]) \
        * FIDUCIAL.h
    pref = prefactor(z_mid, chi_mid)
    d4 = growth(z_mid) ** 4

    print("truth = pref * D^4 * raw z=0 quadrature (certified to 0.5% "
          "on the stored shells at these gammas)\n")
    print(f"{'gamma':>7} {'lambda':>8} {'z':>6}   {'truth':>12} "
          f"{'linear/truth':>13} {'logpchip/truth':>15}")
    worst_lin, worst_log = 0.0, 0.0
    for i, g in enumerate(g_act):
        y = zeta[i]
        f_lin = lambda x: np.interp(x, lam, y)          # noqa: E731
        s = float(np.sign(y[0]))
        f_log = PchipInterpolator(lam, np.log(np.abs(y)))
        for j, lm in enumerate(LAM_MID):
            truth = pref[j] * d4[j] * raw_vertex_z0(g * ARCMIN, chi_mid[j],
                                                    ELL_MAX)
            r_lin = float(f_lin(lm)) / truth
            r_log = s * float(np.exp(f_log(lm))) / truth
            worst_lin = max(worst_lin, abs(np.log(abs(r_lin))))
            worst_log = max(worst_log, abs(np.log(abs(r_log))))
            print(f"{g:7.2f} {lm:8.0f} {z_mid[j]:6.3f}   {truth:+12.4e} "
                  f"{r_lin:13.3f} {r_log:15.3f}")
        print()
    print(f"worst |ln ratio| vs truth: linear {worst_lin:.3f}   "
          f"logpchip {worst_log:.3f}")


if __name__ == "__main__":
    main()
