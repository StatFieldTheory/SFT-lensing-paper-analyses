"""Order-0 by integration by parts -- a route that never differentiates C.

order0_mc = int_a^b W F' dla  with  W = G^2 / D^4,  F = D^4 C(la,la)
          = [W F]_a^b - int_a^b W' F dla
          = - W(a) F(a) - int_a^b W' F dla        (G(lam_f) = 0)

W and its derivative involve only the smooth background, so this evaluates the
SAME anchor integral using C VALUES ONLY.  If the production-vs-corrected shift
agrees with the derivative-based route (t5), the shift is not an artefact of how
the density is differentiated.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402

ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import corr_op_C_callable as _corr_op  # noqa: E402

TABLE_LAM = np.load(_corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)


def main() -> None:
    bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
    lam_f = bgm.LAM_SOURCE_BASELINE
    a = 406.0
    n = 20001
    la = np.linspace(a, lam_f, n)
    D = np.asarray(bg.D(la), float)
    G = ds.order0_window(bg, la, lam_f, n_gauss=256)
    W = G ** 2 / D ** 4
    Wp = np.gradient(W, la, edge_order=2)

    print(f"{'gamma':>8} {'O0 prod (IBP)':>15} {'O0 corr (IBP)':>15} "
          f"{'corr/prod-1 [%]':>16}")
    for gam in (0.5, 1.0, 5.0, 17.3, 44.4, 114.3):
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        n2 = ds._n2_of_cos(cosg)
        C_prod = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0]
                           for t in la])
        node = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0]
                         for t in TABLE_LAM])
        if np.all(node > 0):
            g = PchipInterpolator(TABLE_LAM, np.log(node))
            C_corr = np.exp(g(la))
            model = "pchip-log"
        else:
            g = PchipInterpolator(TABLE_LAM, node)
            C_corr = g(la)
            model = "pchip-lin"  # node values are not all positive at this gamma
        out = []
        for C in (C_prod, C_corr):
            F = D ** 4 * C
            out.append(-W[0] * F[0] - np.trapezoid(Wp * F, la))
        print(f"{gam:>8.2f} {out[0]:>15.6e} {out[1]:>15.6e} "
              f"{100*(out[1]/out[0]-1):>16.3f}   {model}")


if __name__ == "__main__":
    main()
