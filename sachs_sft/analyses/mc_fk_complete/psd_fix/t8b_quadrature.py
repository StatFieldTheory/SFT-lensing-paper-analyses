"""Is the Order-0 shift sensitive to the quadrature, or to the model?

Three evaluations of the SAME integral, per C-diagonal model:
  GL-256   : what driver_stats.order0_mc does (Gauss-Legendre on Sigma2 * G^2)
  TRAPZ    : fine uniform trapezoid of W * F', F' by central differences
  IBP      : fine uniform trapezoid of -W(a)F(a) - W' F   (no derivative of C)

Sigma2 built from the production table is a SAWTOOTH (the bilinear diagonal has
a slope discontinuity at each of the 21 nodes and dips negative just past each),
so a 256-node Gauss-Legendre rule is not obviously converged on it; the
corrected diagonal is smooth and should be easy.  Printing all three separates
"the correction changes the answer" from "the rule cannot integrate a sawtooth".
"""
from __future__ import annotations

import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator

ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402

ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import corr_op_C_callable as _corr_op  # noqa: E402
import sigma2_repaired as R  # noqa: E402
from t5_order0_impact import CorrectedDiag  # noqa: E402

TABLE_LAM = np.load(_corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)
A = 406.0


def c_diag(model: str, cosg: float, la: np.ndarray) -> np.ndarray:
    n2 = ds._n2_of_cos(cosg)
    if model == "prod":
        return np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0]
                         for t in la])
    node = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0]
                     for t in TABLE_LAM])
    g = PchipInterpolator(TABLE_LAM, np.log(node))
    return np.exp(g(la))


def main() -> None:
    bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
    lam_f = bgm.LAM_SOURCE_BASELINE
    for gam in (1.0, 17.3):
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        print(f"\n=== gamma = {gam}' ===")
        print(f"{'model':>6} {'n':>7} {'GL-256':>14} {'TRAPZ':>14} {'IBP':>14}")
        res = {}
        for model in ("prod", "corr"):
            b = (R.R0Reference(background=bg, apply_c0=False) if model == "prod"
                 else CorrectedDiag(background=bg, apply_c0=False))
            gl = float(ds.order0_mc(cosg, lam_f=lam_f, builder=b, n_gauss=256))
            for n in (20001, 80001):
                la = np.linspace(A, lam_f, n)
                D = np.asarray(bg.D(la), float)
                G = ds.order0_window(bg, la, lam_f, n_gauss=256)
                W = G ** 2 / D ** 4
                F = D ** 4 * c_diag(model, cosg, la)
                Fp = np.gradient(F, la, edge_order=2)
                trap = float(np.trapezoid(W * Fp, la))
                ibp = float(-W[0] * F[0] - np.trapezoid(np.gradient(W, la,
                                                                    edge_order=2) * F, la))
                print(f"{model:>6} {n:>7} {gl:>14.6e} {trap:>14.6e} {ibp:>14.6e}")
                res[(model, n, "t")] = trap
                res[(model, n, "i")] = ibp
            res[(model, "gl")] = gl
        for tagn, lbl in (("gl", "GL-256"), ((80001, "t"), "TRAPZ n=80001"),
                          ((80001, "i"), "IBP   n=80001")):
            k = (lambda m: res[(m, tagn)] if tagn == "gl" else res[(m,) + tagn])
            print(f"  corr/prod - 1 [{lbl}] = "
                  f"{100*(k('corr')/k('prod')-1):+.3f} %")


if __name__ == "__main__":
    main()
