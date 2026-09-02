"""Build a corr_op table on a DENSER lambda grid over a restricted range.

Purpose: measure the lambda-interpolation error of the production 21-node
`corr_op_table_zs5_21pt_...npz` inside one of its wide cells, by computing the
underlying C at lambda values that are table NODES of the dense build (hence
exact) and comparing against the production table's bilinear interpolant.

Every build argument except `--lambda-grid-custom` is copied verbatim from
`callables/C_propagator/corr_op/build_corr_op_table.py` (the validated
production setup), so the only difference between the two tables is the
lambda grid (and, as a consequence, the source-adaptive shared chi grid --
which is what the control comparison at the shared nodes measures).

Usage:
    python t1_build_dense.py <tag> <comma-separated lambda nodes> [n_chi]
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

CANOES_ROOT = Path("/Users/zzhang/projects/angular_statistics/canoes")
PY = sys.executable
HERE = Path(__file__).resolve().parent
OUT = HERE / "dense"

PK = "examples/data/PCAMB_pyccl_stf_fid_z0.txt"
COSMO = ["--omega-m", "0.3160919980475834", "--h", "0.6711", "--n-s", "0.97"]


def build(tag: str, nodes: str, n_chi: str = "256") -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"corr_op_dense_{tag}.npz"
    args = [
        PY, "-m", "canoes.sachs.sft_input.corr_op.build",
        "--pk", PK, *COSMO,
        "--ell-max", "5000",
        "--dense-threshold", "500",
        "--n-log-per-decade", "16",
        "--lambda-grid", "custom",
        "--lambda-grid-custom", nodes,
        "--source-coord", "lambda_project_mpc",
        "--chi-min", "50.0",
        "--n-chi-per-pair", n_chi,
        "--accuracy", "default",
        "--ell-hi-threshold", "800",
        "--output", str(out),
    ]
    print("[t1] " + " ".join(args), flush=True)
    env = {**os.environ,
           "CANOES_SUPPRESS_METAL_WARNING": "1",
           "PYTHONPATH": str(CANOES_ROOT / "src")}
    t0 = time.perf_counter()
    rc = subprocess.run(args, cwd=str(CANOES_ROOT), env=env).returncode
    print(f"[t1] tag={tag} rc={rc} wall={time.perf_counter()-t0:.1f}s", flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(build(sys.argv[1], sys.argv[2],
                           sys.argv[3] if len(sys.argv) > 3 else "256"))
