"""Reproducible build of the corr_op C_propagator table (zs5-21pt, STF fiducial).

Runs the canoes corr_op build CLI with the EXACT validated setup (the args that
produced the previously-validated `corr_op_prod_ellmax5000_zs5_21pt_...` table:
PyCCL cross-check at 5-10% on kk/xi_+ Order-0). The rebuilt table cross-checks
bit-for-bit (max abs diff 2.9e-15) against that archived table. The output is
consumed by `corr_op_C_callable.py` in this folder.

Usage:
    python build_corr_op_table.py            # full production build (~1.2 h)
    python build_corr_op_table.py --smoke    # tiny grid wiring check (~30 s)

Convention (read back from the table by the callable, never hand-typed):
    source_coord = lambda_project_mpc   (corr_op's natural source coordinate is
                                         the affine λ; χ is only the inner
                                         line-of-sight integration grid)
    ell_hi_threshold = 800              (S3: Limber for ℓ>800, t-form for ℓ≤800)
    cl_tensor_convention = bare ; ell_taper_frac = 0.2  (CorrOpTable2DConfig defaults)
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
PY = "/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python"
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz"

# z=0 matter P(k), PyCCL STF fiducial. (canoes-relative path; _load_inputs applies
# apply_poisson_to_phi in the physical-Mpc convention.)
PK = "examples/data/PCAMB_pyccl_stf_fid_z0.txt"

COSMO = ["--omega-m", "0.3160919980475834", "--h", "0.6711", "--n-s", "0.97"]
BUILD_ARGS = [
    "--pk", PK,
    *COSMO,
    "--ell-max", "5000",
    "--dense-threshold", "500",
    "--n-log-per-decade", "16",
    "--lambda-grid", "zs5-21pt",          # → source_coord=lambda_project_mpc
    "--chi-min", "50.0",
    "--n-chi-per-pair", "256",
    "--accuracy", "default",              # FFTlog Nmax=256
    "--ell-hi-threshold", "800",          # S3 high-ℓ Limber split
]


def main() -> int:
    smoke = "--smoke" in sys.argv
    args = [PY, "-m", "canoes.sachs.sft_input.corr_op.build"]
    if smoke:
        args += ["--pk", PK, *COSMO, "--smoke-test",
                 "--output", str(HERE / "_smoke_corr_op_table.npz")]
    else:
        args += [*BUILD_ARGS, "--output", str(OUTPUT)]
    print("[build_corr_op_table] running:\n  " + " ".join(args))
    env = {**os.environ, "CANOES_SUPPRESS_METAL_WARNING": "1"}
    return subprocess.run(args, cwd=str(CANOES_ROOT), env=env).returncode


if __name__ == "__main__":
    raise SystemExit(main())
