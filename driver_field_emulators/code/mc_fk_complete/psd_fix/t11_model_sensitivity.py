"""How much of the +6% depends on the choice of 1-D diagonal model?

Recomputes the Order-0 shift by the IBP route with four reconstructions of the
true diagonal from the 21 exact node values, and with the four dense-rebuilt
cells replaced by the MEASURED truth (interpolated inside those cells from the
5 dense nodes) so at least those cells carry no model at all.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402

ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import corr_op_C_callable as _corr_op  # noqa: E402

sys.path.insert(0, "/Users/zzhang/projects/angular_statistics/canoes/src")
from canoes.sachs.sft_input.corr_op.table import (  # noqa: E402
    load_table_2d_from_npz, make_C_fn_lambda_project,
)

HERE = Path(__file__).resolve().parent
TABLE_LAM = np.load(_corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)
DENSE = ("cellA5", "cellB5", "cellC5", "cellD5")


def main(gam: float = 1.0) -> None:
    bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
    lam_f = bgm.LAM_SOURCE_BASELINE
    la = np.linspace(406.0, lam_f, 40001)
    D = np.asarray(bg.D(la), float)
    G = ds.order0_window(bg, la, lam_f, n_gauss=256)
    W = G ** 2 / D ** 4
    Wp = np.gradient(W, la, edge_order=2)
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    n2 = ds._n2_of_cos(cosg)
    C0 = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0] for t in la])
    node = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0]
                     for t in TABLE_LAM])

    def O0(C):
        F = D ** 4 * C
        return -W[0] * F[0] - np.trapezoid(Wp * F, la)

    base = O0(C0)
    models = {
        "pchip-log": np.exp(PchipInterpolator(TABLE_LAM, np.log(node))(la)),
        "cubic-log": np.exp(CubicSpline(TABLE_LAM, np.log(node))(la)),
        "pchip-lin": PchipInterpolator(TABLE_LAM, node)(la),
        "loglinear": np.exp(np.interp(la, TABLE_LAM, np.log(node))),
    }
    # dense-measured truth inside the four rebuilt cells, pchip-log elsewhere
    hyb = models["pchip-log"].copy()
    for tag in DENSE:
        p = HERE / "dense" / f"corr_op_dense_{tag}.npz"
        if not p.exists():
            continue
        tab = load_table_2d_from_npz(str(p))
        lam = np.asarray(tab.cl_table.lambda_grid, float)
        Cd = make_C_fn_lambda_project(tab)
        v = np.array([Cd(ds._N1, float(t), n2, float(t))[0, 0] for t in lam])
        m = (la >= lam[0]) & (la <= lam[-1])
        hyb[m] = np.exp(PchipInterpolator(lam, np.log(v))(la[m]))
    models["hybrid (dense cells measured)"] = hyb

    print(f"gamma = {gam}'   O0(production diagonal) = {base:.6e}")
    print(f"{'diagonal model':>30} {'O0':>14} {'shift [%]':>11}")
    for k, C in models.items():
        v = O0(C)
        print(f"{k:>30} {v:>14.6e} {100*(v/base-1):>11.3f}")


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 1.0)
