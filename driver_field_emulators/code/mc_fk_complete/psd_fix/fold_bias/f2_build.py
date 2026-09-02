"""Build a corr_op table on an arbitrary lambda grid, with a settable chi_min.

Same as psd_fix/t1_build_dense.py except that --chi-min is exposed, because a
node-matched rebuild for the Order-0 fold needs lambda nodes far below the
production chi_min = 50 Mpc floor (the first Gauss-Legendre node is 5.57 Mpc,
and `_build_shared_chi_grid` requires chi_min < min(lambda_grid)).

Usage:  python f2_build.py <tag> <comma-separated nodes> [n_chi] [chi_min]
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

CANOES_ROOT = Path("/Users/zzhang/projects/angular_statistics/canoes")
PY = sys.executable
OUT = Path(__file__).resolve().parent / "tables"
PK = "examples/data/PCAMB_pyccl_stf_fid_z0.txt"
COSMO = ["--omega-m", "0.3160919980475834", "--h", "0.6711", "--n-s", "0.97"]


def build(tag: str, nodes: str, n_chi: str = "256", chi_min: str = "50.0") -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"corr_op_{tag}.npz"
    args = [
        PY, "-m", "canoes.sachs.sft_input.corr_op.build",
        "--pk", PK, *COSMO,
        "--ell-max", "5000",
        "--dense-threshold", "500",
        "--n-log-per-decade", "16",
        "--lambda-grid", "custom",
        "--lambda-grid-custom", nodes,
        "--source-coord", "lambda_project_mpc",
        "--chi-min", chi_min,
        "--n-chi-per-pair", n_chi,
        "--accuracy", "default",
        "--ell-hi-threshold", "800",
        "--lambda-chunk-size", "4",
        "--output", str(out),
    ]
    print("[f2] " + " ".join(args), flush=True)
    env = {**os.environ, "CANOES_SUPPRESS_METAL_WARNING": "1",
           "PYTHONPATH": str(CANOES_ROOT / "src")}
    t0 = time.perf_counter()
    rc = subprocess.run(args, cwd=str(CANOES_ROOT), env=env).returncode
    print(f"[f2] tag={tag} rc={rc} wall={time.perf_counter()-t0:.1f}s", flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(build(*sys.argv[1:]))
