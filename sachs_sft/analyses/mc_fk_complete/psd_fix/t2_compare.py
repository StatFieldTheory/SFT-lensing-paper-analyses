"""Measure the production table's own lambda-interpolation error.

TRUTH  = a corr_op table rebuilt with the IDENTICAL production settings but a
         denser lambda grid inside one cell; evaluated AT its own nodes, so no
         lambda interpolation is involved.
PROD   = the deployed 21-node table's bilinear interpolant at the same lambda.

CONTROL: the dense grid shares its two end nodes with the production grid.  Any
disagreement THERE is not interpolation, it is the source-adaptive shared chi
grid changing when the lambda grid changes.  The control has to be small for
the mid-cell number to mean anything.

Also reports the second-order prediction  delta(u) = h^2 u(1-u) C_xx  which is
the exact leading behaviour of a bilinear interpolant restricted to the cell
diagonal (mixed term cancels; f_xx = f_yy by symmetry).

Usage: python t2_compare.py <tag> [tag ...]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

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
    Cp = make_C_fn_lambda_project(load_table_2d_from_npz(str(PROD)))
    gammas = [(1.0, "cos=1 (within-ray)"),
              (float(np.cos(1.0 * ARCMIN)), "gamma=1'"),
              (float(np.cos(17.3 * ARCMIN)), "gamma=17.3'")]
    for tag in tags:
        path = HERE / "dense" / f"corr_op_dense_{tag}.npz"
        tab = load_table_2d_from_npz(str(path))
        lam = np.asarray(tab.cl_table.lambda_grid, float)
        Cd = make_C_fn_lambda_project(tab)
        print(f"\n################ dense table {path.name}")
        print(f"  dense nodes: {lam.tolist()}")
        for cosg, tag2 in gammas:
            n2 = n2_of_cos(cosg)
            print(f"\n  --- {tag2} ---")
            for comp, cname in (((0, 0), "C00 <Phi00 Phi00>"),
                                ((1, 1), "C11 <RePsi0 RePsi0>")):
                print(f"   {cname}")
                print(f"   {'lambda':>9} {'u':>5} {'TRUTH(dense node)':>19} "
                      f"{'PROD(bilinear)':>16} {'prod/truth-1 [%]':>17}")
                for t in lam:
                    u = (t - lam[0]) / (lam[-1] - lam[0])
                    tr = Cd(N1, float(t), n2, float(t))[comp]
                    pr = Cp(N1, float(t), n2, float(t))[comp]
                    flag = "   <-- CONTROL (shared node)" if (
                        np.isclose(u, 0.0) or np.isclose(u, 1.0)) else ""
                    print(f"   {t:9.2f} {u:5.2f} {tr:19.6e} {pr:16.6e} "
                          f"{100*(pr/tr-1):17.3f}{flag}")


if __name__ == "__main__":
    main(sys.argv[1:] or ["calib3"])
