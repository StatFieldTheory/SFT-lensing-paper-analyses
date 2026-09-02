"""Is C(lam1, lam2) smooth across the diagonal, or does it have a ridge?

Prints the raw (lam1, lam2) matrix of the angular-summed C00 at cos=1 for a
dense table and, on the same nodes, the production table's value.  Entries that
are NODES of BOTH tables are exact in both -- they are controls.  Entries that
are dense-only nodes expose what the production bilinear cell is missing.
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


def main(tag: str) -> None:
    tabd = load_table_2d_from_npz(str(HERE / "dense" / f"corr_op_dense_{tag}.npz"))
    lam = np.asarray(tabd.cl_table.lambda_grid, float)
    Cd = make_C_fn_lambda_project(tabd)
    tabp = load_table_2d_from_npz(str(PROD))
    lamp = np.asarray(tabp.cl_table.lambda_grid, float)
    Cp = make_C_fn_lambda_project(tabp)
    n2 = N1.copy()
    print(f"dense nodes {lam.tolist()}")
    print("\nDENSE C00(lam1,lam2) [exact at every entry]")
    hdr = "        " + "".join(f"{x:13.1f}" for x in lam)
    print(hdr)
    Md = np.array([[Cd(N1, float(a), n2, float(b))[0, 0] for b in lam] for a in lam])
    Mp = np.array([[Cp(N1, float(a), n2, float(b))[0, 0] for b in lam] for a in lam])
    for i, a in enumerate(lam):
        print(f"{a:8.1f}" + "".join(f"{v:13.5e}" for v in Md[i]))
    print("\nPRODUCTION C00 on the same (lam1,lam2), bilinear where not a prod node")
    print(hdr)
    for i, a in enumerate(lam):
        print(f"{a:8.1f}" + "".join(f"{v:13.5e}" for v in Mp[i]))
    print("\nPROD/DENSE - 1 [%]   (entries that are prod nodes in BOTH indices "
          "are controls)")
    print(hdr)
    for i, a in enumerate(lam):
        row = ""
        for j, b in enumerate(lam):
            ctrl = (np.min(np.abs(lamp - a)) < 1e-9) and (np.min(np.abs(lamp - b)) < 1e-9)
            row += f"{100*(Mp[i,j]/Md[i,j]-1):11.3f}{'*' if ctrl else ' '} "
        print(f"{a:8.1f}" + row)
    print("* = both indices are production nodes -> chi-grid control, "
          "no interpolation involved")

    # Off-diagonal decay right next to the diagonal (ridge diagnostic), from the
    # dense table only.
    print("\nridge diagnostic (dense table): C00(lam1, lam1+d)/C00(lam1,lam1)")
    for i, a in enumerate(lam[:-1]):
        for j in range(i, lam.size):
            d = lam[j] - a
            print(f"  lam1={a:7.1f}  d={d:7.1f}  ratio={Md[i, j]/Md[i, i]:8.5f}")
        break


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "calib3")
