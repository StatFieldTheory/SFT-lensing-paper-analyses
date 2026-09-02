"""Helper: run sft-wick FF once for a given kappa2_grid.npz.

Args:
    kappa2_npz : path to the kappa2 grid file
    output_npz : path to write the FF result
    sigma2_npz : optional path to the sigma2 high-ell table

Run as a subprocess to ensure ``kappa2_callable._GRID`` module-level cache is
fresh each time.
"""

from __future__ import annotations

import os
import shutil
import sys
import time
from dataclasses import replace
from pathlib import Path


def _resolve_here_path(path_like: str | os.PathLike | None, here: Path) -> Path | None:
    if path_like is None:
        return None
    p = Path(path_like)
    return p if p.is_absolute() else (here / p).resolve()


def _clear_cache(path_like: str | os.PathLike | None, here: Path) -> None:
    p = _resolve_here_path(path_like, here)
    if p is not None and p.exists():
        shutil.rmtree(p)


def main() -> None:
    if len(sys.argv) not in (3, 4):
        print(
            "usage: run_FF_single.py <kappa2_npz> <output_npz> [sigma2_npz]",
            file=sys.stderr,
        )
        sys.exit(2)
    kappa2_path = Path(sys.argv[1]).resolve()
    output_path = Path(sys.argv[2]).resolve()
    sigma2_path = Path(sys.argv[3]).resolve() if len(sys.argv) == 4 else None

    here = Path(__file__).resolve().parent
    pkg = here.parent
    repo_root = pkg.parent.parent

    # Keep the canoes-pipeline script directory first so all helper modules
    # used by config_FF_canoes.yaml resolve from this local pipeline.
    for d in (repo_root, here):
        d_str = str(d)
        if d_str in sys.path:
            sys.path.remove(d_str)
        sys.path.insert(0, d_str)

    os.environ["SFT_WICK_KAPPA2_NPZ"] = str(kappa2_path)
    if sigma2_path is None:
        os.environ.pop("SFT_WICK_SIGMA2_NPZ", None)
        os.environ["SFT_WICK_GL_INCLUDE_SIGMA2"] = "0"
    else:
        os.environ["SFT_WICK_SIGMA2_NPZ"] = str(sigma2_path)
        os.environ["SFT_WICK_GL_INCLUDE_SIGMA2"] = "1"

    # Pre-resolve SFT_WICK_KAPPA3_NPZ relative to the caller's CWD before the
    # chdir below.  kappa3_callable runs inside the joblib workers whose CWD
    # has been switched to ``here`` (= scripts/canoes_pipeline/scripts/), so a
    # relative ``inputs/foo.npz`` from the user (CWD = canoes_pipeline/) would
    # otherwise be mis-resolved to canoes_pipeline/scripts/inputs/foo.npz.
    kappa3_env = os.environ.get("SFT_WICK_KAPPA3_NPZ")
    if kappa3_env:
        os.environ["SFT_WICK_KAPPA3_NPZ"] = str(
            Path(kappa3_env).expanduser().resolve()
        )

    from sft_wick.workflow.config import (  # type: ignore[import-not-found]
        load_workflow_config,
        run_workflow,
    )

    # YAML config selectable via SFT_WICK_CONFIG_YAML env var (default
    # config_FF_canoes.yaml for the dense kappa2 path). For the chi_delta
    # path set SFT_WICK_CONFIG_YAML=config_FF_chi_delta.yaml.
    cfg_yaml_name = os.environ.get(
        "SFT_WICK_CONFIG_YAML", "config_FF_canoes.yaml"
    )
    cfg_path = (here / cfg_yaml_name).resolve()
    # Everything relative in the YAML resolves against the YAML's own folder:
    # sft-wick resolves module and data paths that way at load time, and its
    # output and cache paths are CWD-relative, so the CWD is set to that folder.
    base = cfg_path.parent
    os.chdir(base)
    cfg = load_workflow_config(cfg_path)
    _clear_cache(cfg.expand.cache_path, base)
    _clear_cache(cfg.propagators.cache_path, base)
    sweep_n_gauss = os.environ.get("SFT_WICK_SWEEP_N_GAUSS")
    if sweep_n_gauss:
        cfg = replace(
            cfg,
            sweep=replace(cfg.sweep, n_gauss=int(sweep_n_gauss)),
        )
    expand_orders = os.environ.get("SFT_WICK_EXPAND_ORDERS")
    if expand_orders:
        orders = tuple(int(part.strip()) for part in expand_orders.split(",") if part.strip())
        cfg = replace(
            cfg,
            expand=replace(cfg.expand, orders=orders),
        )
    # Resolve the actual output path from the loaded YAML so this helper works
    # for any config (config_FF_canoes.yaml, config_FF_full_iso_F.yaml, etc.)
    # rather than only the dense-canoes default. YAML output paths are
    # relative to the YAML's parent dir (we already chdir'd to ``here``).
    yaml_out = Path(cfg.output[0].path)
    default_out = yaml_out if yaml_out.is_absolute() else (base / yaml_out).resolve()
    if default_out.exists():
        # A run folder carrying a PRODUCTION marker holds a deployed result:
        # refuse to unlink it unless the caller says so explicitly. A variant
        # config that inherited a production output path destroyed the
        # deployed FK sweep on 2026-08-25.
        if (default_out.parent / "PRODUCTION").exists() and not os.environ.get("SFT_WICK_FORCE_OVERWRITE"):
            sys.exit(f"[run_FF_single] refusing to overwrite {default_out}: its folder is marked "
                     "PRODUCTION; set SFT_WICK_FORCE_OVERWRITE=1 to override")
        default_out.unlink()

    print(f"[run_FF_single] kappa2 = {kappa2_path.name}")
    if sigma2_path is None:
        print("[run_FF_single] sigma2 = disabled")
    else:
        print(f"[run_FF_single] sigma2 = {sigma2_path.name}")
    print(f"[run_FF_single] sweep.n_gauss = {cfg.sweep.n_gauss}")
    print(f"[run_FF_single] expand.orders = {cfg.expand.orders}")
    print(f"[run_FF_single] -> {output_path}")
    t0 = time.perf_counter()
    run_workflow(cfg)
    elapsed = time.perf_counter() - t0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(default_out), str(output_path))
    print(f"[run_FF_single] done in {elapsed:.1f} s")


if __name__ == "__main__":
    main()
