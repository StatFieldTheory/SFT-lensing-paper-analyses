"""(a) how the mid-cell error shrinks when the cell is halved, measured directly
   (b) which lambda cells carry the Order-0 shift.

(a) uses ONLY the dense tables: their own nodes bracket the cell midpoint at half
    the production spacing, so
        B_half = 0.25 C(a',a') + 0.5 C(a',b') + 0.25 C(b',b')
    with a', b' the dense nodes either side of the midpoint, is exactly what the
    production callable would return if the table had twice the resolution.
    Comparing B_half and B_full to the same truth C(m,m) gives the empirical
    h-scaling of the defect (a smooth function would give 4x, a diagonal ridge
    only 2x).

(b) uses the exact IBP decomposition
        dO0 = - int W' (F_corr - F_prod) dlam
    restricted to one cell at a time; the cell sum is the total by construction.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
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
N1 = np.array([0.0, 0.0, 1.0])


def part_a() -> None:
    print("=== (a) mid-cell bilinear error vs cell width (cos = 1, C00) ===")
    print(f"{'cell':>16} {'h_full':>7} {'err(h)':>9} {'h_half':>7} {'err(h/2)':>9} "
          f"{'ratio':>7}")
    Cp = make_C_fn_lambda_project(load_table_2d_from_npz(str(_corr_op.TABLE_PATH)))
    for tag in ("cellA5", "cellB5", "cellC5", "cellD5"):
        p = HERE / "dense" / f"corr_op_dense_{tag}.npz"
        if not p.exists():
            continue
        tab = load_table_2d_from_npz(str(p))
        lam = np.asarray(tab.cl_table.lambda_grid, float)
        Cd = make_C_fn_lambda_project(tab)
        a, b = float(lam[0]), float(lam[-1])
        m = float(lam[len(lam) // 2])
        ap, bp = float(lam[len(lam) // 2 - 1]), float(lam[len(lam) // 2 + 1])
        truth = Cd(N1, m, N1, m)[0, 0]
        # full-width production cell (bilinear from the production table)
        full = Cp(N1, m, N1, m)[0, 0]
        # half-width cell, built from the dense table's own exact entries
        half = (0.25 * Cd(N1, ap, N1, ap)[0, 0]
                + 0.5 * Cd(N1, ap, N1, bp)[0, 0]
                + 0.25 * Cd(N1, bp, N1, bp)[0, 0])
        ef, eh = 100 * (full / truth - 1), 100 * (half / truth - 1)
        print(f"{a:7.1f}-{b:7.1f} {b-a:7.1f} {ef:8.2f}% {bp-ap:7.1f} {eh:8.2f}% "
              f"{ef/eh:7.2f}")


def part_b(gam: float) -> None:
    bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
    lam_f = bgm.LAM_SOURCE_BASELINE
    la = np.linspace(406.0, lam_f, 40001)
    D = np.asarray(bg.D(la), float)
    G = ds.order0_window(bg, la, lam_f, n_gauss=256)
    W = G ** 2 / D ** 4
    Wp = np.gradient(W, la, edge_order=2)
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    n2 = ds._n2_of_cos(cosg)
    C = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0] for t in la])
    node = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0]
                     for t in TABLE_LAM])
    g = PchipInterpolator(TABLE_LAM, np.log(node))
    F0, F1 = D ** 4 * C, D ** 4 * np.exp(g(la))
    O0 = -W[0] * F0[0] - np.trapezoid(Wp * F0, la)
    print(f"\n=== (b) per-cell Order-0 attribution at gamma = {gam}' "
          f"(O0_prod = {O0:.6e}) ===")
    print(f"{'cell':>16} {'h':>6} {'dO0/O0 [%]':>11}")
    tot = 0.0
    for lo, hi in zip(TABLE_LAM[:-1], TABLE_LAM[1:]):
        m = (la >= lo) & (la <= hi)
        if m.sum() < 2:
            continue
        d = -np.trapezoid(Wp[m] * (F1[m] - F0[m]), la[m]) / O0 * 100
        tot += d
        print(f"{lo:7.1f}-{hi:7.1f} {hi-lo:6.1f} {d:11.3f}")
    print(f"{'total':>23} {tot:11.3f}")


if __name__ == "__main__":
    part_a()
    part_b(1.0)
