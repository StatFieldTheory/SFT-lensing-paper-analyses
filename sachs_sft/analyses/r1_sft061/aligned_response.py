"""Cached Sachs response on the exact lightcone map used by the C table.

All distances are physical Mpc. The project affine convention has E0=1.
The Jacobi amplitude D=a*chi is the manuscript's verified FLRW solution.
Native canoes interpolation knots are retained, avoiding a second resampling
grid or the older D_callable's different z_max=10 tabulation.
"""

from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
SACHS = HERE.parents[1]
CANOES_ROOT = Path(os.environ.get("CANOES_ROOT", "~/projects/angular_statistics/canoes")).expanduser().resolve()
sys.path.insert(0, str(CANOES_ROOT / "src"))

from canoes.cosmo import lightcone as _lightcone  # noqa: E402

OMEGA_M = 0.3160919980475834
H_REDUCED = 0.6711
BACKGROUND_Z_MAX = 1200.0
BACKGROUND_N_Z = 5001
OBSERVER_E0 = 1.0
COSMOLOGY = SimpleNamespace(Omega_m=OMEGA_M, h=H_REDUCED)
SOURCE_REDSHIFTS = np.array([1.0, 1.7, 2.5, 3.2, 4.0, 5.0])

# The private table accessor is used only to retain the exact breakpoints of
# the public maps. Coordinate conversion itself uses their public API, and
# validate() checks the cached implementation against those maps directly.
_Z, _CHI, _, _ = _lightcone._zchi_growth_tables(
    COSMOLOGY, z_max=BACKGROUND_Z_MAX, n_z=BACKGROUND_N_Z,
    distance_unit="Mpc",
)
_LAMBDA = _lightcone.lambda_of_chi(
    COSMOLOGY, _CHI, convention="project", z_max=BACKGROUND_Z_MAX,
    distance_unit="Mpc",
)
MAX_AFFINE_MPC = float(_LAMBDA[-1])
for _table in (_Z, _CHI, _LAMBDA):
    _table.setflags(write=False)


def _validated_affine(value: np.ndarray | float) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(array)) or np.any(array < 0) or np.any(array > MAX_AFFINE_MPC):
        raise ValueError(f"Affine coordinates must be finite and within [0, {MAX_AFFINE_MPC}] Mpc")
    return array


def chi_of_lambda(value: np.ndarray | float) -> np.ndarray | float:
    """Exact piecewise-linear public canoes affine-to-chi map, cached."""
    array = _validated_affine(value)
    result = np.interp(array, _LAMBDA, _CHI)
    return float(result) if array.ndim == 0 else result


def z_of_lambda(value: np.ndarray | float) -> np.ndarray | float:
    """Use the same two interpolation stages as canoes's public maps."""
    chi = np.asarray(chi_of_lambda(value))
    result = np.interp(chi, _CHI, _Z)
    return float(result) if chi.ndim == 0 else result


def D_at(value: np.ndarray | float) -> np.ndarray | float:
    """Background Jacobi amplitude D=a*chi in physical Mpc."""
    chi = np.asarray(chi_of_lambda(value))
    redshift = np.interp(chi, _CHI, _Z)
    result = chi / (1.0 + redshift)
    return float(result) if chi.ndim == 0 else result


@lru_cache(maxsize=16384)
def _D_scalar(value: float) -> float:
    if not np.isfinite(value) or not 0.0 <= value <= MAX_AFFINE_MPC:
        raise ValueError(f"Affine coordinate {value} is outside the cached lightcone")
    chi = float(np.interp(value, _LAMBDA, _CHI))
    return chi / (1.0 + float(np.interp(chi, _CHI, _Z)))


def response_R(t: np.ndarray | float, lambda_prime: np.ndarray | float) -> np.ndarray | float:
    """Retarded response Theta(t-lambda_prime) [D(lambda_prime)/D(t)]^2.

    The coincident positive-time value is one, matching the existing response
    callable. The observer value is zero. Acausal or pre-observer pairs return
    zero without evaluating D outside its physical domain.
    """
    if np.ndim(t) == 0 and np.ndim(lambda_prime) == 0:
        final, source = float(t), float(lambda_prime)
        if not np.isfinite(final) or not np.isfinite(source):
            raise ValueError("Response times must be finite")
        if source < 0.0 or final <= 0.0 or final < source:
            return 0.0
        denominator = _D_scalar(final)
        return (_D_scalar(source) / denominator) ** 2
    final, source = np.broadcast_arrays(np.asarray(t, dtype=float), np.asarray(lambda_prime, dtype=float))
    if np.any(~np.isfinite(final)) or np.any(~np.isfinite(source)):
        raise ValueError("Response times must be finite")
    result = np.zeros_like(final)
    causal = (source >= 0.0) & (final > 0.0) & (final >= source)
    if np.any(causal):
        result[causal] = (np.asarray(D_at(source[causal])) / np.asarray(D_at(final[causal]))) ** 2
    return result


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate() -> dict:
    """Check public mapping parity and the actual C-builder response factor."""
    affine_sources = np.interp(SOURCE_REDSHIFTS, _Z, _LAMBDA)
    probes = np.r_[0.0, 50.0, 406.0, 550.0, 1000.0, 1700.0, 2300.0, affine_sources]
    public_chi = _lightcone.chi_of_lambda(
        COSMOLOGY, probes, convention="project", z_max=BACKGROUND_Z_MAX,
        distance_unit="Mpc",
    )
    public_z = _lightcone.z_of_chi(
        COSMOLOGY, public_chi, z_max=BACKGROUND_Z_MAX, distance_unit="Mpc",
    )
    np.testing.assert_allclose(chi_of_lambda(probes), public_chi, rtol=0, atol=0)
    np.testing.assert_allclose(z_of_lambda(probes), public_z, rtol=0, atol=0)
    np.testing.assert_allclose(z_of_lambda(affine_sources), SOURCE_REDSHIFTS, rtol=0, atol=5e-13)
    np.testing.assert_allclose(D_at(probes), public_chi / (1 + public_z), rtol=0, atol=0)

    from canoes.sachs.sft_input.corr_op._limber import _id_basis_windows
    from canoes.nuell.transfer._sachs_decomp import get_sachs_decomposition
    from canoes.nuell.transfer.sachs_driver import COEFS_PSI0

    inner_affine = np.array([50.0, 406.0, 1000.0, 1700.0, 2200.0, 2320.0])
    inner_chi = np.asarray(chi_of_lambda(inner_affine))
    source_chi = np.asarray(chi_of_lambda(affine_sources))
    lightcone, inner_opz = _lightcone.build_lightcone(
        COSMOLOGY, inner_chi, z_max=BACKGROUND_Z_MAX, distance_unit="Mpc",
    )
    windows, _ = _id_basis_windows(
        decomp=get_sachs_decomposition(COEFS_PSI0), chi=inner_chi,
        lam=source_chi, opz_chi=inner_opz,
        opz_lam=1.0 + np.asarray(z_of_lambda(affine_sources)),
        m_id=lightcone.M[0], lightcone=lightcone, l2_arr=np.array([6.0]),
    )
    from_c_builder = windows[0] * inner_chi[None, :]**2 / (
        inner_opz[None, :]**2 * lightcone.M[0][None, :]
    )
    response = response_R(affine_sources[:, None], inner_affine[None, :])
    np.testing.assert_allclose(response, from_c_builder, rtol=3e-14, atol=1e-15)
    scalar = np.array([[response_R(float(t), float(s)) for s in inner_affine] for t in affine_sources])
    np.testing.assert_allclose(response, scalar, rtol=2e-15, atol=1e-15)
    np.testing.assert_array_equal(response_R([0, 10, -1, 10], [0, 20, 0, -1]), np.zeros(4))
    np.testing.assert_array_equal(response_R(affine_sources, affine_sources), np.ones(6))
    return {
        "status": "passed",
        "public_coordinate_maps": "Exact equality at observer, selected table nodes and all requested source planes",
        "D_identity": "Exact equality to public chi/(1+z) at all probes",
        "C_builder_response_max_abs_error": float(np.max(np.abs(response - from_c_builder))),
        "C_builder_helper": "canoes.sachs.sft_input.corr_op._limber._id_basis_windows with COEFS_PSI0, stripping its recorded growth and density factors",
        "scalar_array_max_abs_error": float(np.max(np.abs(response - scalar))),
        "causal_and_observer_zero": "passed",
        "coincident_positive_time_unity": "passed",
    }


def write_source_map(path: Path = HERE / "source_map.json") -> None:
    """Write pinned source geometry and validation without running a fold."""
    if path.exists():
        raise FileExistsError(f"Preserve the existing source map before writing {path}")
    validation = validate()
    affine = np.interp(SOURCE_REDSHIFTS, _Z, _LAMBDA)
    source_files = [Path(__file__).resolve(), Path(_lightcone.__file__).resolve(),
                    CANOES_ROOT / "src/canoes/cosmo/growth.py",
                    CANOES_ROOT / "src/canoes/sachs/sft_input/corr_op/table.py",
                    CANOES_ROOT / "src/canoes/sachs/sft_input/corr_op/_limber.py"]
    table_path = SACHS / "callables/C_propagator/corr_op/corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz"
    commit = subprocess.run(["git", "-C", str(CANOES_ROOT), "rev-parse", "HEAD"],
                            text=True, capture_output=True, check=True).stdout.strip()
    payload = {
        "z": SOURCE_REDSHIFTS.tolist(),
        "lambda_mpc": affine.tolist(),
        "chi_mpc": np.asarray(chi_of_lambda(affine)).tolist(),
        "D_mpc": np.asarray(D_at(affine)).tolist(),
        "mapping_provenance": {
            "cosmology": {"Omega_m": OMEGA_M, "h": H_REDUCED, "geometry": "flat matter-plus-Lambda, without radiation"},
            "E0": OBSERVER_E0,
            "convention": "project",
            "distance_unit": "physical Mpc",
            "background_z_max": BACKGROUND_Z_MAX,
            "background_n_z": BACKGROUND_N_Z,
            "native_knots": int(_Z.size),
            "max_affine_mpc": MAX_AFFINE_MPC,
            "source_mapping": "Same public canoes lambda_of_chi, chi_of_lambda and z_of_chi as the C-table build, cached on their native interpolation knots",
            "canoes_root": str(CANOES_ROOT),
            "canoes_commit": commit,
            "source_sha256": {str(file): _sha256(file) for file in source_files},
            "C_table_path": str(table_path),
            "C_table_sha256": _sha256(table_path),
            "jacobi_identity_source": "sections/appendix.tex, appendix D equals a chi, verified by mathematica/derive_jacobi_flat_flrw.wl",
        },
        "response_module": str(Path(__file__).resolve()),
        "response_attr": "response_R",
        "validation": validation,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"path": str(path), "z": payload["z"], "lambda_mpc": payload["lambda_mpc"], "validation": validation}, indent=2))


__all__ = ["D_at", "chi_of_lambda", "z_of_lambda", "response_R", "validate", "MAX_AFFINE_MPC"]


if __name__ == "__main__":
    write_source_map()
