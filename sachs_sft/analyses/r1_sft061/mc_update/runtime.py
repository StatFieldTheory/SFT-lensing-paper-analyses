"""Process-local adapters for the existing Sachs Monte Carlo implementations."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
R1 = HERE.parent
ANALYSES = R1.parent
SACHS = ANALYSES.parent
DEFAULT_RESPONSE = R1 / "aligned_response.py"
DEFAULT_TABLE = SACHS / "callables/kappa3_vertex/rebuild/products/pieces_nphi512/table_permclosed_np512.npz"
DEFAULT_FOLD = R1 / "true_redshift_fk_gl24_corrected_main/main/xi_main.npz"
SOURCE = 2317.9695958203943


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure(table: Path, response: Path, endpoint: float) -> SimpleNamespace:
    """Finish legacy imports, then redirect inputs in this child process only."""
    for folder in (ANALYSES / "mc_fk_complete", ANALYSES / "mc_fk_complete/inputmatch"):
        sys.path.insert(0, str(folder))
    import fk_expect_exact as ex
    import fk_complete_core as fc
    import ff_functional as ff
    import ff_anti
    import perm_aware_kappa3_callable as pa

    spec = importlib.util.spec_from_file_location("r1_mc_aligned_response", response.resolve())
    if spec is None or spec.loader is None:
        raise ImportError(response)
    aligned = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(aligned)
    if not np.isclose(float(aligned.z_of_lambda(endpoint)), 5.0, rtol=0, atol=1e-10):
        raise ValueError("This update requires the actual z_s=5 endpoint")

    class Background:
        """Use the existing aligned map directly, without another spline."""
        def __init__(self, lam_min=120.0, lam_max=2330.0, n_grid=2048):
            self.lam_min = float(lam_min)
            self.lam_max = float(lam_max)

        D = staticmethod(aligned.D_at)
        z = staticmethod(aligned.z_of_lambda)
        chi = staticmethod(aligned.chi_of_lambda)
        response_R = staticmethod(aligned.response_R)

        @staticmethod
        def a(lam):
            return 1.0 / (1.0 + np.asarray(aligned.z_of_lambda(lam)))

    # No subsequent call to _bootstrap.wire is allowed: it would reset K.
    ex._BG.Background = Background
    ex._DS._DEFAULT_BUILDER = None
    pa.TABLE_PATH = table.resolve(strict=True)
    pa._CACHE = None
    ex._DS._k3 = pa
    ff.LF = float(endpoint)
    ff.BG = Background()
    assert ex._DS._k3.TABLE_PATH == table.resolve()
    modules = [ex, fc, ff, ff_anti, ex._CORE, ex._DS, pa, ex._DS._corr_op]
    sources = [Path(module.__file__).resolve() for module in modules]
    sources += [response.resolve(), Path(__file__).resolve()]
    covariance = SACHS / "callables/C_propagator/corr_op/corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz"
    provenance = {
        "source_redshift": 5.0, "lambda_source_mpc": endpoint,
        "lower_support_mpc": 406.0, "ell_max_K": 15360,
        "table_path": str(table.resolve()), "table_sha256": sha256(table),
        "response_path": str(response.resolve()), "response_sha256": sha256(response),
        "covariance_path": str(covariance), "covariance_sha256": sha256(covariance),
        "source_sha256": {str(path): sha256(path) for path in sources},
        "Sigma2": "160-node cubic spline derivative of D^4 times equal-time corr_op, endpoint clipping, no anchor",
        "units": "physical Mpc, Born/FLRW numerical convention E0=1",
    }
    return SimpleNamespace(ex=ex, fc=fc, ff=ff, anti=ff_anti, core=ex._CORE,
                           aligned=aligned, provenance=provenance)


def check_response(runtime, lam, resp) -> float:
    expected = runtime.aligned.response_R(lam[1:], lam[:-1])
    np.testing.assert_allclose(resp[1:], expected, rtol=3e-15, atol=1e-15)
    return float(np.max(np.abs(resp[1:] - expected)))


def fold_values(path: Path, gammas, endpoint: float, table: Path, response: Path):
    """Read the new raw FK fold, with explicit interpolation provenance."""
    prepared_path = path.parent.parent / "prepared.json"
    prepared = json.loads(prepared_path.read_text())
    targets = [target for target in prepared["targets"]
               if Path(target["result"]).resolve() == path.resolve()]
    if len(targets) != 1 or Path(targets[0]["table"]).resolve() != table.resolve():
        raise ValueError("FK comparison preparation does not identify the same K input")
    if Path(prepared["source_mapping"]["response_module"]).resolve() != response.resolve():
        raise ValueError("FK comparison preparation uses a different response")
    pinned = {Path(item["path"]).resolve(): item["sha256"] for item in prepared["inputs"]}
    for input_path in (table, response):
        if pinned.get(input_path.resolve()) != sha256(input_path):
            raise ValueError(f"FK comparison preparation hash mismatch: {input_path}")
    with np.load(path, allow_pickle=True) as product:
        mask = (product["a"] == 0) & (product["b"] == 0) & (product["order"] == 2)
        final = np.asarray(product["t_final"], dtype=float)[mask]
        if not final.size or not np.allclose(final, endpoint, rtol=0, atol=1e-10):
            raise ValueError("FK comparison fold has a different source endpoint")
        x = np.asarray(list(product["x"][mask]), dtype=float)
        y = np.asarray(list(product["y"][mask]), dtype=float)
        cosines = np.sum(x * y, axis=1) / (np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1))
        angles = np.rad2deg(np.arccos(np.clip(cosines, -1, 1))) * 60
        order = np.argsort(angles)
        angles = angles[order]
        values = np.asarray(product["value"], dtype=float)[mask][order]
    if np.any(np.diff(angles) <= 0) or not np.all(np.isfinite(values)):
        raise ValueError("Comparison fold must contain one finite value per angle")
    requested = np.asarray(gammas, dtype=float)
    if requested.min() < angles[0] or requested.max() > angles[-1]:
        raise ValueError("Cannot extrapolate FK comparison beyond its angular grid")
    return np.interp(requested, angles, values), {
        "fold_path": str(path.resolve()), "fold_sha256": sha256(path),
        "fold_prepared_path": str(prepared_path.resolve()),
        "fold_prepared_sha256": sha256(prepared_path),
        "fold_interpolation": "linear in angular separation, matching historical plotting convention",
    }
