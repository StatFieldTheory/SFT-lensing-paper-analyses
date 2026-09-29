"""Run isolated SFT-Wick Order-0 folds at the C table's actual source plane.

The production covariance table is read without modification. Its angular
window is ell=2..5000 with the recorded 20 percent cosine taper. The source
affine coordinate is obtained from the same canoes lightcone mapping that
built the table. Each invocation uses a fresh output directory and one worker.
"""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
import yaml
from scipy.optimize import brentq

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SACHS = HERE.parents[2]
PRODUCTION = SACHS / "sftwick_outputs/2PCF/C_corr_op_O0/config_L2.yaml"
C_MODULE = SACHS / "callables/C_propagator/corr_op/corr_op_C_callable.py"
C_TABLE = C_MODULE.parent / "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz"
CANOES = Path(os.environ.get("CANOES_ROOT", "~/projects/angular_statistics/canoes")).expanduser().resolve()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _commit(path: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def _tracked_state(path: Path) -> dict:
    """Fingerprint tracked working changes without archiving unrelated local files."""
    status = subprocess.run(
        ["git", "-C", str(path), "status", "--porcelain", "--untracked-files=no"],
        capture_output=True, check=True,
    ).stdout
    diff = subprocess.run(
        ["git", "-C", str(path), "diff", "HEAD", "--binary"],
        capture_output=True, check=True,
    ).stdout
    return {"tracked_dirty": bool(status.strip()),
            "tracked_status": status.decode(),
            "tracked_diff_sha256": hashlib.sha256(diff).hexdigest()}


def source_geometry(redshift: float) -> dict:
    """Invert the production table's own monotone lightcone mapping."""
    sys.path.insert(0, str(CANOES / "src"))
    from canoes.cosmo.lightcone import chi_of_lambda, z_of_chi

    with np.load(C_TABLE) as saved:
        config = json.loads(saved["config_json"].item())
        table_metadata = json.loads(saved["cl_table_metadata_json"].item())
        nodes = np.asarray(saved["lambda_grid"], dtype=float)
        ell = np.asarray(saved["ell_assembly_grid"], dtype=int)
    cosmology = config["cosmology"]
    cosmo = SimpleNamespace(**cosmology)

    def mapped_redshift(affine: float) -> float:
        chi = chi_of_lambda(
            cosmo, np.array([affine]), convention="project",
            distance_unit="Mpc", z_max=1200.0,
        )
        return float(z_of_chi(cosmo, chi, distance_unit="Mpc", z_max=1200.0)[0])

    affine = brentq(
        lambda value: mapped_redshift(value) - redshift,
        float(nodes[0]), float(nodes[-1]), xtol=1e-10,
    )
    chi = float(chi_of_lambda(
        cosmo, np.array([affine]), convention="project",
        distance_unit="Mpc", z_max=1200.0,
    )[0])
    bracket_index = int(np.searchsorted(nodes, affine))
    return {
        "z_source_requested": redshift,
        "z_source_recovered": mapped_redshift(affine),
        "lambda_source_mpc": affine,
        "chi_source_mpc": chi,
        "lambda_source_bracket_mpc": nodes[bracket_index - 1:bracket_index + 1].tolist(),
        "lambda_mapping": "canoes project-affine mapping, physical Mpc, z_max=1200, E0=1",
        "table_config": config,
        "table_metadata": table_metadata,
        "ell_min": int(ell.min()),
        "ell_max": int(ell.max()),
        "cosmology": cosmology,
        "low_affine_policy": "The actual make_C_fn_sftwick factory returns zero if either affine coordinate is below the first table node at 406 Mpc. Integrate only the supported interval without bridging the jump at that node.",
        "outer_integral": "SFT-Wick integrate_over=all integrates both external Sachs fields with unit affine weights. No R propagator enters the zero-vertex contraction. R is already folded into C.",
    }


def _make_config(directory: Path, n_gauss: int, source_affine: float, limit_angles: int | None) -> Path:
    config = yaml.safe_load(PRODUCTION.read_text())
    config["system"]["linear"]["R_time_module"] = str(SACHS / "scripts/D_callable.py")
    config["system"]["noise"]["kappa2"]["module"] = str(SACHS / "scripts/kappa2_callable.py")
    config["system"]["vertices"][0]["coupling_path"] = str(SACHS / "scripts/inputs/F_tensor.npy")
    config["expand"]["orders"] = [0]
    config["expand"]["n_jobs"] = 1
    config["expand"]["cache_path"] = str(directory / "cache_expand")
    propagators = config["propagators"]
    propagators["n_jobs"] = 1
    propagators["cache_path"] = str(directory / "cache_propagators")
    propagators["c_closed_form_module"] = str(C_MODULE)
    propagators["c_closed_form_attr"] = "C_fn"
    propagators["c_closed_form_vectorized"] = False
    config["sweep"]["t_final_grid"] = [source_affine]
    config["sweep"]["n_gauss"] = n_gauss
    config["sweep"]["n_jobs"] = 1
    if limit_angles is not None:
        config["sweep"]["positions_grid"]["y"] = config["sweep"]["positions_grid"]["y"][:limit_angles]
    config["output"] = [{"type": "npz", "path": str(directory / "order0_raw.npz")}]
    path = directory / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


def _save_observables(raw_path: Path, output_path: Path) -> None:
    # The object arrays are produced locally by this exact framework run.
    with np.load(raw_path, allow_pickle=True) as raw:
        a, b = np.asarray(raw["a"], int), np.asarray(raw["b"], int)
        order = np.asarray(raw["order"], int)
        values = np.real_if_close(np.asarray(raw["value"]))
        if np.iscomplexobj(values) or not np.all(np.isfinite(values)):
            raise ValueError("Order-0 output is not finite and real")
        x = np.stack(raw["x"]).astype(float)
        y = np.stack(raw["y"]).astype(float)
        x /= np.linalg.norm(x, axis=1)[:, None]
        y /= np.linalg.norm(y, axis=1)[:, None]
        angles = np.rad2deg(np.arccos(np.clip(np.sum(x * y, axis=1), -1, 1))) * 60.0
        channels = {}
        gamma = None
        for first, second in [(0, 0), (0, 1), (1, 1), (2, 2)]:
            mask = (a == first) & (b == second) & (order == 0)
            sort = np.argsort(angles[mask])
            current_gamma = angles[mask][sort]
            if gamma is None:
                gamma = current_gamma
            elif not np.array_equal(current_gamma, gamma):
                raise ValueError("Component angle grids disagree")
            channels[first, second] = values[mask][sort]
        np.savez_compressed(
            output_path,
            gamma_arcmin=gamma,
            theta_arcmin=gamma,
            xi_kappa=channels[0, 0],
            xi_plus=channels[1, 1] + channels[2, 2],
            xi_minus=channels[1, 1] - channels[2, 2],
            xi_kappa_gamma=-channels[0, 1],
            kk=channels[0, 0],
            xip=channels[1, 1] + channels[2, 2],
            xim=channels[1, 1] - channels[2, 2],
            kgt=-channels[0, 1],
            lambda_source_mpc=np.asarray(raw["t_final"], float)[0],
        )


def _integrate_fixed_table(config_path: Path, raw_path: Path) -> dict:
    """Independently integrate the same zero-vertex covariance contraction.

    A product trapezoid rule containing every affine interpolation knot is
    exact for the fixed table's bilinear interpolant. This checks the outer
    GL rule only, without claiming convergence of the underlying covariance.
    """
    spec = importlib.util.spec_from_file_location("corr_op_C_callable", C_MODULE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {C_MODULE}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    config = yaml.safe_load(config_path.read_text())
    final = float(config["sweep"]["t_final_grid"][0])
    with np.load(C_TABLE) as saved:
        knots = np.asarray(saved["lambda_grid"], dtype=float)
    knots = np.r_[knots[knots < final], final]
    widths = np.diff(knots)
    weights = np.r_[widths[0] / 2, (widths[:-1] + widths[1:]) / 2, widths[-1] / 2]
    first, second = np.meshgrid(knots, knots, indexing="ij")
    rows = {key: [] for key in ("x", "y", "t_final", "a", "b", "order", "value")}
    for x in config["sweep"]["positions_grid"]["x"]:
        for y in config["sweep"]["positions_grid"]["y"]:
            # Use the scalar factory, whose angular channels retain the exact
            # input direction. The current batch factory rounds cos(theta)
            # before constructing those channels, changing small angles.
            covariance = np.array([
                module.C_fn(np.asarray(x), float(t1), np.asarray(y), float(t2))
                for t1, t2 in zip(first.ravel(), second.ravel(), strict=True)
            ]).reshape(len(knots), len(knots), 3, 3)
            integrated = np.einsum("i,ijab,j->ab", weights, covariance, weights)
            for a, b in config["sweep"]["component_pairs"]:
                for key, value in zip(
                    rows, (x, y, final, a, b, 0, float(integrated[a, b])), strict=True,
                ):
                    rows[key].append(value)
    np.savez_compressed(raw_path, **{key: np.asarray(value) for key, value in rows.items()})
    return {
        "integration_method": "Product trapezoid over every bilinear interpolation cell",
        "integration_knots_mpc": knots.tolist(),
        "n_affine_knots": len(knots),
        "validation_scope": "Exact integration of this fixed piecewise-bilinear C callable above its 406 Mpc support floor. The unsupported interval contributes exactly zero. Does not validate the table's physical or numerical accuracy.",
        "angular_geometry_policy": "Scalar C_fn uses unrounded input directions. Batch factory rounding is bypassed.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z-source", type=float, default=5.0)
    parser.add_argument("--n-gauss", type=int, default=24)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--limit-angles", type=int, help="Small smoke run only")
    parser.add_argument("--describe", action="store_true", help="Print geometry without running the fold")
    parser.add_argument("--table-integral-check", action="store_true", help="Independent exact integral of the fixed bilinear C interpolant")
    args = parser.parse_args()
    if args.n_gauss < 2 or not np.isfinite(args.z_source) or args.z_source <= 0:
        parser.error("Require n_gauss >= 2 and a finite positive source redshift")
    if args.limit_angles is not None and not 1 <= args.limit_angles <= 40:
        parser.error("limit-angles must be between 1 and 40")
    geometry = source_geometry(args.z_source)
    if args.describe:
        print(json.dumps(geometry, indent=2))
        return
    default_name = "order0_table_integral" if args.table_integral_check else f"order0_gl{args.n_gauss}"
    directory = (args.output_dir or HERE / "outputs" / default_name).resolve()
    if not directory.is_relative_to(HERE):
        parser.error(f"Output must stay inside {HERE}")
    directory.mkdir(parents=True, exist_ok=False)
    config_path = _make_config(directory, args.n_gauss, geometry["lambda_source_mpc"], args.limit_angles)
    initial_hash = _sha256(C_TABLE)
    metadata = {
        **geometry,
        "n_gauss": args.n_gauss,
        "n_jobs": 1,
        "limit_angles": args.limit_angles,
        "production_config": str(PRODUCTION),
        "production_config_sha256": _sha256(PRODUCTION),
        "covariance_table": str(C_TABLE),
        "covariance_table_sha256": initial_hash,
        "covariance_callable_sha256": _sha256(C_MODULE),
        "runner_sha256": _sha256(Path(__file__)),
        "canoes_commit": _commit(CANOES),
        "canoes_working_state": _tracked_state(CANOES),
        "python": sys.executable,
        "status": "running",
        "integration_method": "SFT-Wick external-time Gauss-Legendre rule",
        "angular_geometry_policy": "Scalar C_fn uses unrounded input directions. Batch factory rounding is bypassed.",
    }
    metadata_path = directory / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    from sft_wick.workflow.config import load_workflow_config, run_workflow
    import sft_wick

    metadata["sft_wick_source"] = str(Path(sft_wick.__file__).resolve())
    metadata["sft_wick_commit"] = _commit(Path(sft_wick.__file__).resolve().parents[2])
    metadata["sft_wick_version"] = version("sft-wick")
    metadata["sft_wick_working_state"] = _tracked_state(Path(sft_wick.__file__).resolve().parents[2])
    start = time.perf_counter()
    if args.table_integral_check:
        metadata.update(_integrate_fixed_table(config_path, directory / "order0_raw.npz"))
    else:
        run_workflow(load_workflow_config(config_path), progress=True)
    _save_observables(directory / "order0_raw.npz", directory / "order0_observables.npz")
    if _sha256(C_TABLE) != initial_hash:
        raise RuntimeError("Production covariance table changed during the fold")
    metadata.update(status="complete", elapsed_seconds=time.perf_counter() - start)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"output_dir": str(directory), "elapsed_seconds": metadata["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
