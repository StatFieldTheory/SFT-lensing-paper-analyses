"""Can the TRUE diagonal C(lam,lam) be recovered from the production nodes?

The dense rebuilds give the truth at sub-cell lambda.  This scores three cheap
one-dimensional models built from the PRODUCTION node values alone:

  bilinear   : what the production callable actually returns
               B(u) = (1-u)^2 F_ii + 2u(1-u) F_ij + u^2 F_jj      (the Bezier)
  linear     : straight line through the two diagonal node values
  loglinear  : straight line through ln C  (geometric interpolation)
  pchip-log  : PCHIP through ln C over ALL 21 diagonal node values

If one of these tracks the dense truth to <~1%, the full-range true diagonal can
be reconstructed without rebuilding 20 cells, and the Order-0 impact of the
21-node grid can be evaluated directly.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, "/Users/zzhang/projects/angular_statistics/canoes/src")
from canoes.sachs.sft_input.corr_op.table import (  # noqa: E402
    load_table_2d_from_npz, make_C_fn_lambda_project,
)

HERE = Path(__file__).resolve().parent
PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
            "sachs_sft/callables/C_propagator/corr_op/"
            "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz")
N1 = np.array([0.0, 0.0, 1.0])
ARCMIN = np.pi / (180.0 * 60.0)


def n2_of_cos(c: float) -> np.ndarray:
    c = float(np.clip(c, -1.0, 1.0))
    return np.array([np.sqrt(max(0.0, 1.0 - c * c)), 0.0, c])


def main(tags: list[str]) -> None:
    tabp = load_table_2d_from_npz(str(PROD))
    lamp = np.asarray(tabp.cl_table.lambda_grid, float)
    Cp = make_C_fn_lambda_project(tabp)
    for cosg, gtag in ((1.0, "cos=1"), (float(np.cos(17.3 * ARCMIN)), "gamma=17.3'")):
        n2 = n2_of_cos(cosg)
        diag_p = np.array([Cp(N1, float(t), n2, float(t))[0, 0] for t in lamp])
        pch = PchipInterpolator(lamp, np.log(diag_p))
        print(f"\n=== {gtag} : C00 diagonal models vs dense truth ===")
        print(f"{'lambda':>9} {'truth':>13} {'bilinear':>10} {'linear':>10} "
              f"{'loglin':>10} {'pchip-log':>10}   (all columns: model/truth-1 [%])")
        for tag in tags:
            p = HERE / "dense" / f"corr_op_dense_{tag}.npz"
            if not p.exists():
                continue
            tabd = load_table_2d_from_npz(str(p))
            lam = np.asarray(tabd.cl_table.lambda_grid, float)
            Cd = make_C_fn_lambda_project(tabd)
            print(f"  -- {tag}")
            for t in lam:
                tr = Cd(N1, float(t), n2, float(t))[0, 0]
                bil = Cp(N1, float(t), n2, float(t))[0, 0]
                i = int(np.clip(np.searchsorted(lamp, t) - 1, 0, lamp.size - 2))
                a, b = lamp[i], lamp[i + 1]
                u = (t - a) / (b - a)
                lin = (1 - u) * diag_p[i] + u * diag_p[i + 1]
                log = np.exp((1 - u) * np.log(diag_p[i]) + u * np.log(diag_p[i + 1]))
                pc = float(np.exp(pch(t)))
                print(f"{t:9.2f} {tr:13.5e} {100*(bil/tr-1):10.3f} "
                      f"{100*(lin/tr-1):10.3f} {100*(log/tr-1):10.3f} "
                      f"{100*(pc/tr-1):10.3f}")


if __name__ == "__main__":
    main(sys.argv[1:] or ["calib3"])
