"""Run the production FK two-point fold against an alternative kappa3 table.

The production callable pins its table by an explicit ``TABLE_PATH``
constant, deliberately: sachs_sft resolves products by path, never by
environment variable, so that any result is traceable to one file. To
compare table variants without weakening that rule, this script MATERIALISES
one callable module per variant (a verbatim copy of the production callable
with TABLE_PATH rewritten) and one config per variant, so every run still
points at a concrete, auditable pair of files.

Usage (from ``sachs_sft/callables/kappa3_vertex/rebuild``)::

    SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
    $SFT run_fk_variant.py products/<table>.npz

Writes ``<table stem>_xi.npz`` next to the table, plus the generated
callable and config in a ``_variant_<stem>`` folder beside it.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

SACHS_SFT = Path(__file__).resolve().parents[3]
CALLABLE_DIR = SACHS_SFT / "callables" / "kappa3_vertex" / "equal_time_limber"
PRODUCTION_CALLABLE = CALLABLE_DIR / "equal_time_limber_kappa3_callable.py"
PRODUCTION_CONFIG = (
    SACHS_SFT / "sftwick_outputs" / "2PCF" / "C_corr_op_K_limber_FK"
    / "config_L2.yaml"
)


def materialise_variant(table: Path | None, work_dir: Path,
                        t_final: list[float] | None = None,
                        n_jobs: int | None = None,
                        base_config: Path | None = None,
                        orders: str = "0,2",
                        callable_src: Path | None = None) -> tuple[Path | None, Path]:
    """Write a callable + config pair bound to ``table``.

    ``t_final`` replaces the sweep's source-distance grid. The propagator
    tables are built to ``t_max`` independently of it, so a list of source
    distances costs one sweep rather than one run each. This is what makes
    the source-redshift dependence affordable to measure.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    base = base_config or PRODUCTION_CONFIG

    # An Order-0 config has no three-point vertex, so there is no table to
    # substitute. Driving that config is how the Order-0 reference is obtained
    # at source distances other than the single production one.
    variant_callable = None
    if table is not None:
        # `callable_src` swaps the vertex callable itself, which is how the
        # permutation-aware replacement is exercised through the production
        # runner rather than through a re-implementation of the fold.
        src = callable_src or PRODUCTION_CALLABLE
        source = src.read_text()
        markers = [
            ('TABLE_PATH = _HERE / '
             '"equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz"'),
            'TABLE_PATH = _HERE / "table_permclosed.npz"',
        ]
        marker = next((m for m in markers if m in source), None)
        if marker is None:
            raise RuntimeError(f"{src.name}: no recognised TABLE_PATH line")
        variant_callable = work_dir / f"callable_{work_dir.name}.py"
        variant_callable.write_text(
            source.replace(marker, f'TABLE_PATH = Path(r"{table}")'))

    config = yaml.safe_load(base.read_text())
    # The run folders write their configs with paths relative to the YAML folder
    # (since 2026-09-02); the variant config lives elsewhere, so every inherited
    # module and data path is made absolute here before anything is rewritten.
    def _abs(p):
        return p if Path(p).is_absolute() else str((base.parent / p).resolve())
    sysd = config["system"]
    sysd["linear"]["R_time_module"] = _abs(sysd["linear"]["R_time_module"])
    sysd["noise"]["kappa2"]["module"] = _abs(sysd["noise"]["kappa2"]["module"])
    for v in sysd.get("vertices", []):
        if "coupling_path" in v:
            v["coupling_path"] = _abs(v["coupling_path"])
    for v in sysd.get("nonlocal_vertices", []):
        if "coupling_module" in v:
            v["coupling_module"] = _abs(v["coupling_module"])
    if config.get("propagators", {}).get("c_closed_form_module"):
        config["propagators"]["c_closed_form_module"] = _abs(config["propagators"]["c_closed_form_module"])
    if variant_callable is not None:
        for vertex in config["system"].get("nonlocal_vertices", []):
            if vertex.get("name") == "K":
                vertex["coupling_module"] = str(variant_callable)

    # Redirect every write target into the variant folder. This is not a
    # tidiness measure: the production config's output.path is ABSOLUTE, and
    # run_FF_single.py unlinks that path before the run and moves the result
    # out of it afterwards. Inheriting it would delete the deployed FK 2PCF
    # (it did once, 2026-08-25; recovered from the _PROD_REF copy). The
    # expand and propagator caches are likewise cleared by the runner, so
    # they must not point at the production directory either.
    for entry in config.get("output", []):
        entry["path"] = str(work_dir / Path(entry["path"]).name)
    for section in ("expand", "propagators"):
        cache = config.get(section, {}).get("cache_path")
        if cache:
            config[section]["cache_path"] = str(work_dir / Path(cache).name)

    if t_final is not None:
        config["sweep"]["t_final_grid"] = [float(v) for v in t_final]

    # Worker count is a MEMORY control, not just a speed knob. sft-wick forks
    # with loky and every worker reloads each input NPZ, so the budget is
    # per-worker footprint times n_jobs. The production config asks for -1,
    # i.e. every core, which is only safe when nothing else is running. An
    # earlier session lost a machine this way.
    if n_jobs is not None:
        for section in ("propagators", "sweep"):
            if section in config and "n_jobs" in config[section]:
                config[section]["n_jobs"] = int(n_jobs)

    variant_config = work_dir / f"config_{work_dir.name}.yaml"
    variant_config.write_text(yaml.safe_dump(config, sort_keys=False))

    production_dir = base.resolve().parent
    for entry in config.get("output", []):
        if Path(entry["path"]).resolve().parent == production_dir:
            raise RuntimeError("variant output still points at the production run")
    return variant_callable, variant_config


def run(config: Path, output: Path, n_gauss: int,
        orders: str = "0,2") -> None:
    """Drive the production runner, not a re-implementation of it.

    ``scripts/run_FF_single.py`` is the script that produced the deployed FK
    2PCF, including its module-cache hygiene (it is meant to be run as a
    subprocess so the callables' module-level table caches start clean).
    Reusing it keeps a variant run byte-comparable with production apart
    from the table under test.
    """
    runner = SACHS_SFT / "scripts" / "run_FF_single.py"
    kappa2_npz = SACHS_SFT / "scripts" / "inputs" / "kappa2_grid_lcut0_low_empty.npz"
    env = dict(os.environ)
    env.update(
        SFT_WICK_CONFIG_YAML=str(config),
        SFT_WICK_EXPAND_ORDERS=orders,
        SFT_WICK_SWEEP_N_GAUSS=str(n_gauss),
    )
    command = [sys.executable, str(runner), str(kappa2_npz), str(output)]
    print("[fk] " + " ".join(command))
    completed = subprocess.run(command, env=env, cwd=str(config.parent))
    if completed.returncode != 0:
        sys.exit(f"run_FF_single failed with code {completed.returncode}")
    print(f"[fk] wrote {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path, nargs="?", default=None,
                        help="kappa3 table; omit when driving an Order-0 config")
    parser.add_argument("--n-gauss", type=int, default=24)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--t-final", type=float, nargs="+", default=None,
                        help="source affine distances lambda(z_s); default "
                             "keeps the production z_s = 5 value")
    parser.add_argument("--tag", default=None,
                        help="suffix distinguishing runs on the same table")
    parser.add_argument("--config", type=Path, default=None,
                        help="base config to copy; defaults to the production "
                             "FK run's config_L2.yaml")
    parser.add_argument("--orders", default="0,2")
    parser.add_argument("--callable", dest="callable_src", type=Path, default=None,
                        help="vertex callable to use instead of the production "
                             "one; its TABLE_PATH is rewritten to `table`")
    parser.add_argument("--work-dir", type=Path, default=None)
    parser.add_argument("--n-jobs", type=int, default=8,
                        help="loky workers for the propagator and sweep stages. "
                             "Bounds memory: the production -1 takes every core "
                             "and every worker reloads the inputs.")
    args = parser.parse_args()

    table = args.table.resolve() if args.table else None
    if table is not None and not table.exists():
        sys.exit(f"table not found: {table}")
    stem = (table.stem if table else "order0") + (f"_{args.tag}" if args.tag else "")
    default_dir = (table.parent if table else Path.cwd())
    work_dir = args.work_dir or default_dir / f"_variant_{stem}"
    _, config = materialise_variant(table, work_dir, args.t_final, args.n_jobs,
                                    args.config, args.orders,
                                    args.callable_src)
    output = args.out or (default_dir / f"{stem}_xi.npz")
    run(config, output, args.n_gauss, args.orders)


if __name__ == "__main__":
    main()
