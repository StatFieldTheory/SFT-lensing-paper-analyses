"""Measure whether the FF run's coarser propagator time grid matters.

The three deployed runs share every geometric setting except one: the FF
run tabulates its propagators with dt = 2.0 while Order-0 and FK use
dt = 1.0. The paper sums the three curves, so an inconsistency in the
propagator grid is a difference between terms that are meant to be
directly comparable. This is independent of the kappa3 vertex, since the
FF sweep uses vertex_types [F] and never queries the kappa3 table.

This script re-runs the FF sweep with dt = 1.0, everything else identical,
into an isolated directory, so the two can be differenced. It never writes
into the deployed run directory.

Usage::

    SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
    $SFT run_ff_dt_check.py --dt 1.0 --n-jobs 6
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import yaml

_REPO = Path(__file__).resolve().parents[5]
SACHS_SFT = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft"
FF_CONFIG = (SACHS_SFT / "sftwick_outputs" / "2PCF"
             / "C_corr_op_K_limber_FF" / "config_L2.yaml")
PRODUCTS = Path(__file__).resolve().parent / "products"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--n-jobs", type=int, default=6,
                    help="cap sft-wick parallelism; the machine is shared")
    ap.add_argument("--n-gauss", type=int, default=24)
    args = ap.parse_args()

    config = yaml.safe_load(FF_CONFIG.read_text())
    work = PRODUCTS / f"_ff_dt{args.dt:g}"
    work.mkdir(parents=True, exist_ok=True)

    config["propagators"]["dt"] = float(args.dt)
    config["propagators"]["n_jobs"] = int(args.n_jobs)
    config["sweep"]["n_jobs"] = int(args.n_jobs)
    for entry in config.get("output", []):
        entry["path"] = str(work / Path(entry["path"]).name)
    for section in ("expand", "propagators"):
        cache = config.get(section, {}).get("cache_path")
        if cache:
            config[section]["cache_path"] = str(work / Path(cache).name)

    production_dir = FF_CONFIG.parent.resolve()
    for entry in config.get("output", []):
        if Path(entry["path"]).resolve().parent == production_dir:
            raise SystemExit("variant output still points at the production run")

    config_path = work / f"config_ff_dt{args.dt:g}.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    output = PRODUCTS / f"xi_FF_dt{args.dt:g}.npz"
    env = dict(os.environ)
    env.update(SFT_WICK_CONFIG_YAML=str(config_path),
               SFT_WICK_EXPAND_ORDERS="0,2",
               SFT_WICK_SWEEP_N_GAUSS=str(args.n_gauss))
    command = [sys.executable, str(SACHS_SFT / "scripts" / "run_FF_single.py"),
               str(SACHS_SFT / "scripts" / "inputs" / "kappa2_grid_lcut0_low_empty.npz"),
               str(output)]
    print("[ff] " + " ".join(command), flush=True)
    completed = subprocess.run(command, env=env, cwd=str(work))
    if completed.returncode != 0:
        sys.exit(f"run_FF_single failed with code {completed.returncode}")
    print(f"[ff] wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
