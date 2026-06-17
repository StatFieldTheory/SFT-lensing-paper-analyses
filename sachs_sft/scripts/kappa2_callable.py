"""sft-wick ``noise.kappa2`` source-cumulant module (full 3-component).

Component layout:
    a = 0  ->  Phi_00       (Ricci focusing, spin-0)
    a = 1  ->  Re Psi_0     (Weyl shear +, spin-2 in global tangent basis)
    a = 2  ->  Im Psi_0     (Weyl shear x, spin-2 in global tangent basis)

For each call ``kappa2(n1, t1, n2, t2) -> ndarray((3, 3))`` the wrapper
- maps the affine times ``t1, t2`` to comoving distances ``chi1, chi2``,
- finds the nearest (or interpolated) shell indices in the precomputed
  ``Cl_*(ell, chi_a, chi_b)`` tables,
- computes the great-circle frame angles ``(psi_1, psi_2)``,
- evaluates ``xi_plus, xi_minus, xi_TE`` from the C_l tables and Wigner-d,
- assembles the global-frame 3x3 block via ``spin_rotation``.

The .npz expected at the path specified by the environment variable
``SFT_WICK_KAPPA2_NPZ`` (default ``scripts/canoes_pipeline/inputs/kappa2_grid.npz``)
must contain:

    ells              (N_l,)   int     multipoles, sorted
    lambda_shells_Mpc (N_chi,) float   monotonically increasing
    chi_shells_Mpc    (N_chi,) float
    Cl_phi00_phi00    (N_l, N_chi, N_chi) float
    Cl_phi00_psi0     (N_l, N_chi, N_chi) float
    Cl_psi0_psi0      (N_l, N_chi, N_chi) float
"""

from __future__ import annotations

import os
import sys
import json
from pathlib import Path
from typing import Final

import numpy as np

# Allow this module to be loaded both as a package member and as a standalone file via
# sft-wick's ``importlib.util.spec_from_file_location``; in the latter case
# there is no parent package, so we resolve ``spin_rotation`` by adding the
# directory to ``sys.path``.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from spin_rotation import (  # noqa: E402
    GreatCircleAngles,
    great_circle_angles,
    scalar_spin2_pair,
    spin2_real_block,
    xi_plus_minus,
    xi_TE_scalar,
)

N_COMPONENTS: Final[int] = 3
DEFAULT_NPZ_RELATIVE: Final[str] = "scripts/canoes_pipeline/inputs/kappa2_grid.npz"
REWEIGHT_ENV: Final[str] = "SFT_WICK_KAPPA2_RADIAL_REWEIGHT"
REQUIRE_OPERATOR_ENV: Final[str] = "SFT_WICK_KAPPA2_REQUIRE_OPERATOR"
ROW_OPERATOR_SCHEMA: Final[str] = "canoes_kappa2_row_delta_bin_integral"
ROW_OPERATOR_COMPONENT_ORDER: Final[tuple[str, ...]] = (
    "phi00_phi00",
    "phi00_psi0",
    "psi0_phi00",
    "psi0_psi0",
)
UNSAFE_EXACT_SAMPLE_CONVENTIONS: Final[set[str]] = {
    "chi_sample_density",
    "lambda_sample_of_chi_density",
    "thin_shell_point",
}


def _resolve_npz_path() -> Path:
    env = os.environ.get("SFT_WICK_KAPPA2_NPZ")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    project_root = here.parent.parent.parent.parent
    return (project_root / DEFAULT_NPZ_RELATIVE).resolve()


_GRID: dict | None = None


def _decode_cosmo_meta(data: np.lib.npyio.NpzFile) -> dict:
    if "cosmo_meta" not in data.files:
        return {}
    raw = data["cosmo_meta"]
    try:
        if isinstance(raw, np.ndarray):
            if raw.dtype == np.uint8:
                blob = raw.tobytes()
            else:
                blob = raw.item() if raw.shape == () else bytes(raw.tolist())
        else:
            blob = raw
        if isinstance(blob, str):
            return json.loads(blob)
        if isinstance(blob, (bytes, bytearray)):
            return json.loads(bytes(blob).decode("utf-8"))
    except Exception:
        return {"cosmo_meta_decode_error": True}
    return {}


def _infer_radial_kernel_convention(meta: dict) -> str:
    explicit = meta.get("radial_kernel_convention")
    if explicit:
        return str(explicit)
    if meta.get("limber_delta_measure") == "lambda":
        return "lambda_delta_density"
    if meta.get("producer") == "canoes.sachs.compute_kappa2_grid":
        return "thin_shell_point"
    return "unknown"


def _decode_text_scalar(value) -> str:
    arr = np.asarray(value)
    if arr.shape == ():
        item = arr.item()
    else:
        item = arr.tolist()
    if isinstance(item, bytes):
        return item.decode("utf-8")
    return str(item)


def _detect_radial_operator_schema(
    data: np.lib.npyio.NpzFile,
    meta: dict,
) -> str:
    if "radial_operator_schema" in data.files:
        return _decode_text_scalar(data["radial_operator_schema"])
    schema = meta.get("radial_operator_schema", "")
    return "" if schema is None else str(schema)


def _env_bool(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def _dlambda_dchi_from_grid(data: np.lib.npyio.NpzFile, meta: dict) -> np.ndarray:
    """Estimate d lambda / d chi on the loaded shell centers."""
    if "d_lambda_d_chi" in meta:
        return np.asarray(meta["d_lambda_d_chi"], dtype=float)
    if "z_shells" in data.files:
        z = np.asarray(data["z_shells"], dtype=float)
        convention = str(meta.get("lambda_convention", "project"))
        if convention in {"project", "manuscript", "unknown"}:
            return (1.0 + z) ** -2
        if convention == "legacy":
            return 1.0 + z
    lam = np.asarray(data["lambda_shells_Mpc"], dtype=float)
    chi = np.asarray(data["chi_shells_Mpc"], dtype=float)
    return np.gradient(lam, chi)


def _apply_radial_reweight(grid: dict, data: np.lib.npyio.NpzFile) -> None:
    """Optional cheap Jacobian experiments for exact sample grids.

    These modes are diagnostics only; they avoid rebuilding canoes grids.
    """
    mode = os.environ.get(REWEIGHT_ENV, "").strip().lower()
    if mode in {"", "none", "0", "false"}:
        grid["radial_reweight"] = "none"
        return
    jac = _dlambda_dchi_from_grid(data, grid["cosmo_meta"])
    if np.any(jac <= 0.0):
        raise ValueError("d lambda / d chi must be positive for radial reweighting")

    if mode in {"sqrt_dlambda_dchi", "delta_lambda", "born_limber_delta"}:
        leg = np.sqrt(jac)
    elif mode in {"sqrt_dchi_dlambda", "inverse_delta_lambda"}:
        leg = 1.0 / np.sqrt(jac)
    elif mode in {"product_dlambda_dchi", "smooth_inverse"}:
        leg = jac
    elif mode in {"product_dchi_dlambda", "smooth_chi_to_lambda"}:
        leg = 1.0 / jac
    else:
        raise ValueError(
            f"unsupported {REWEIGHT_ENV}={mode!r}; expected one of "
            "sqrt_dlambda_dchi, sqrt_dchi_dlambda, "
            "product_dlambda_dchi, product_dchi_dlambda"
        )

    weights = np.outer(leg, leg)
    for key in ("Cl_phi00_phi00", "Cl_phi00_psi0", "Cl_psi0_psi0"):
        grid[key] *= weights[None, :, :]
    grid["radial_reweight"] = mode
    grid["radial_reweight_description"] = (
        f"Applied {mode} using d lambda / d chi shell factors from metadata/z_shells."
    )


def _load_grid() -> dict:
    global _GRID
    if _GRID is not None:
        return _GRID
    path = _resolve_npz_path()
    if not path.exists():
        raise FileNotFoundError(
            f"kappa2 grid not found at {path}. "
            "Run scripts/canoes_pipeline/build_cl.py first."
        )
    data = np.load(path, allow_pickle=False)
    meta = _decode_cosmo_meta(data)
    schema = _detect_radial_operator_schema(data, meta)
    if schema == ROW_OPERATOR_SCHEMA:
        _GRID = _load_row_operator_grid(path, data, meta)
        return _GRID

    radial_kernel_convention = _infer_radial_kernel_convention(meta)
    if (
        _env_bool(REQUIRE_OPERATOR_ENV)
        and radial_kernel_convention in UNSAFE_EXACT_SAMPLE_CONVENTIONS
    ):
        radial_discretization = meta.get("radial_discretization", "unknown")
        raise ValueError(
            f"Unsupported exact kappa2 radial convention in {path}:\n"
            f"  radial_kernel_convention={radial_kernel_convention!r}\n"
            f"  radial_discretization={radial_discretization!r}\n"
            "c_fast_gl production mode requires either\n"
            f"  radial_operator_schema={ROW_OPERATOR_SCHEMA!r}\n"
            "or a validated lambda bin/delta operator.\n"
            "Rebuild with scripts/canoes_pipeline/build_row_operator_cl.py."
        )

    _GRID = {
        "format": "dense_cl_v1",
        "path": str(path),
        "ells": np.asarray(data["ells"], dtype=int),
        "lambda_shells": np.asarray(data["lambda_shells_Mpc"], dtype=float),
        "chi_shells": np.asarray(data["chi_shells_Mpc"], dtype=float),
        "Cl_phi00_phi00": np.asarray(data["Cl_phi00_phi00"], dtype=float),
        "Cl_phi00_psi0": np.asarray(data["Cl_phi00_psi0"], dtype=float),
        "Cl_psi0_psi0": np.asarray(data["Cl_psi0_psi0"], dtype=float),
        "cosmo_meta": meta,
        "radial_kernel_convention": radial_kernel_convention,
        "scale_split_role": str(meta.get("scale_split_role", "")),
    }
    # Trivial-zero short-circuit for L_cut=0 stubs: when the writer advertises
    # ``ell_window='low_empty'`` (or simply stores zero ells) every ell-sum
    # is mathematically zero, so the (cos_gamma) caches downstream can be
    # collapsed to a single shared zero tensor.
    _GRID["is_trivial_zero"] = (
        _GRID["ells"].size == 0
        or str(meta.get("ell_window", "")) == "low_empty"
    )
    _apply_radial_reweight(_GRID, data)
    # Shape sanity.
    n_l = _GRID["ells"].size
    n_chi = _GRID["lambda_shells"].size
    for key in ("Cl_phi00_phi00", "Cl_phi00_psi0", "Cl_psi0_psi0"):
        if _GRID[key].shape != (n_l, n_chi, n_chi):
            raise ValueError(
                f"kappa2 grid {key} has shape {_GRID[key].shape}, expected {(n_l, n_chi, n_chi)}"
            )
    return _GRID


def _load_row_operator_grid(
    path: Path,
    data: np.lib.npyio.NpzFile,
    meta: dict,
) -> dict:
    version = int(np.asarray(data["radial_operator_schema_version"]).item())
    if version != 1:
        raise ValueError(
            f"unsupported row-operator schema version {version} in {path}; "
            "expected version 1"
        )
    component_order = tuple(
        _decode_text_scalar(v) for v in np.asarray(data["component_order"])
    )
    if component_order != ROW_OPERATOR_COMPONENT_ORDER:
        raise ValueError(
            f"unsupported row-operator component_order={component_order!r} "
            f"in {path}; expected {ROW_OPERATOR_COMPONENT_ORDER!r}"
        )
    op_bin_n_quad = (
        np.asarray(data["op_bin_n_quad"], dtype=np.int16)
        if "op_bin_n_quad" in data.files
        else np.zeros(data["op_bin_integral"].shape[0], dtype=np.int16)
    )

    grid = {
        "format": "row_operator_v1",
        "path": str(path),
        "radial_operator_schema": ROW_OPERATOR_SCHEMA,
        "radial_operator_schema_version": version,
        "component_order": component_order,
        "ells": np.asarray(data["ells"], dtype=int),
        "lambda_rows": np.asarray(data["lambda_rows_Mpc"], dtype=float),
        "chi_rows": np.asarray(data["chi_rows_Mpc"], dtype=float),
        "z_rows": np.asarray(data["z_rows"], dtype=float),
        "row_edges": np.asarray(data["row_edges_Mpc"], dtype=float),
        "row_weights": np.asarray(data["row_weights_Mpc"], dtype=float),
        "op_row_ptr": np.asarray(data["op_row_ptr"], dtype=np.int64),
        "op_bin_ell_index": np.asarray(data["op_bin_ell_index"], dtype=np.int32),
        "op_bin_row_index": np.asarray(data["op_bin_row_index"], dtype=np.int32),
        "op_bin_lambda_left": np.asarray(
            data["op_bin_lambda_left_Mpc"], dtype=float,
        ),
        "op_bin_lambda_right": np.asarray(
            data["op_bin_lambda_right_Mpc"], dtype=float,
        ),
        "op_bin_lambda_eff": np.asarray(
            data["op_bin_lambda_eff_Mpc"], dtype=float,
        ),
        "op_bin_n_quad": op_bin_n_quad,
        "op_bin_integral": np.asarray(data["op_bin_integral"], dtype=float),
        "op_bin_abs_integral_estimate": np.asarray(
            data["op_bin_abs_integral_estimate"], dtype=float,
        ),
        "cosmo_meta": meta,
        "radial_kernel_convention": str(
            meta.get("radial_kernel_convention", "lambda_row_bin_integral")
        ),
        "scale_split_role": str(meta.get("scale_split_role", "")),
    }
    n_ell = int(grid["ells"].size)
    n_row = int(grid["lambda_rows"].size)
    n_bin = int(grid["op_bin_integral"].shape[0])
    if grid["row_edges"].shape != (n_row + 1,):
        raise ValueError(
            f"row_edges_Mpc has shape {grid['row_edges'].shape}, expected {(n_row + 1,)}"
        )
    if grid["row_weights"].shape != (n_row,):
        raise ValueError(
            f"row_weights_Mpc has shape {grid['row_weights'].shape}, expected {(n_row,)}"
        )
    if grid["op_row_ptr"].shape != (n_ell * n_row + 1,):
        raise ValueError(
            "op_row_ptr has shape "
            f"{grid['op_row_ptr'].shape}, expected {(n_ell * n_row + 1,)}"
        )
    if grid["op_row_ptr"][0] != 0 or grid["op_row_ptr"][-1] != n_bin:
        raise ValueError("op_row_ptr endpoints do not match op_bin_integral length")
    if np.any(np.diff(grid["op_row_ptr"]) < 0):
        raise ValueError("op_row_ptr must be non-decreasing")
    for key in (
        "op_bin_ell_index",
        "op_bin_row_index",
        "op_bin_lambda_left",
        "op_bin_lambda_right",
        "op_bin_lambda_eff",
        "op_bin_n_quad",
    ):
        if grid[key].shape != (n_bin,):
            raise ValueError(f"{key} has shape {grid[key].shape}, expected {(n_bin,)}")
    if "op_bin_n_quad" in data.files and np.any(grid["op_bin_n_quad"] < 1):
        raise ValueError("op_bin_n_quad must be positive when present")
    if grid["op_bin_integral"].shape != (n_bin, len(ROW_OPERATOR_COMPONENT_ORDER)):
        raise ValueError(
            "op_bin_integral has shape "
            f"{grid['op_bin_integral'].shape}, expected "
            f"{(n_bin, len(ROW_OPERATOR_COMPONENT_ORDER))}"
        )
    if "op_bin_target_index" in data.files:
        target_index = np.asarray(data["op_bin_target_index"], dtype=np.int32)
        if target_index.shape != (n_bin,):
            raise ValueError(
                "op_bin_target_index has shape "
                f"{target_index.shape}, expected {(n_bin,)}"
            )
        if np.any((target_index < 0) | (target_index >= n_row)):
            raise ValueError("op_bin_target_index contains out-of-range row indices")
        grid["op_bin_target_index"] = target_index
    return grid


def _nearest_shell(lambda_query: float, lambda_shells: np.ndarray) -> int:
    return int(np.argmin(np.abs(lambda_shells - float(lambda_query))))


def _bracket_shell(lam: float, lambda_shells: np.ndarray) -> tuple[int, int, float]:
    """Linear-interpolation bracket for ``lam`` against the (sorted) shell grid.

    Returns ``(i_lo, i_hi, w_hi)`` where the interpolated value is
    ``(1 - w_hi) * f[i_lo] + w_hi * f[i_hi]``. Outside the grid the result is
    clamped (constant extrapolation, ``w_hi`` set to 0 or 1 accordingly).

    The smooth (linear) interpolation is essential for sft-wick's adaptive
    ``dblquad`` to converge: a nearest-neighbour lookup would make kappa2 a
    step function in (lambda_1, lambda_2), causing the quadrature to refine
    indefinitely at every shell boundary.
    """
    n = lambda_shells.size
    lam_f = float(lam)
    if lam_f <= lambda_shells[0]:
        return 0, 0, 0.0
    if lam_f >= lambda_shells[-1]:
        return n - 1, n - 1, 0.0
    i_hi = int(np.searchsorted(lambda_shells, lam_f, side="right"))
    i_hi = min(i_hi, n - 1)
    i_lo = i_hi - 1
    span = lambda_shells[i_hi] - lambda_shells[i_lo]
    w_hi = 0.0 if span <= 0 else float((lam_f - lambda_shells[i_lo]) / span)
    return i_lo, i_hi, w_hi


def _to_3d(n) -> np.ndarray:
    """Embed an arbitrary-dimensional unit vector in R^3 by padding with zeros.

    sft-wick's ``homogeneity: rotation`` builder passes 2-D vectors
    ``(1, 0)`` and ``(cos_gamma, sin_gamma)`` for the C-table sweep. The
    (3, 3) kappa2 matrix here is then expressed in the canonical equatorial
    frame (n in the xy-plane); rotation invariance makes the choice of frame
    irrelevant for ``kappa2[0, 0]`` (Phi_00 auto, frame-independent) and
    consistent for the spin-2 entries.
    """
    arr = np.asarray(n, dtype=float).reshape(-1)
    if arr.size >= 3:
        return arr[:3]
    out = np.zeros(3, dtype=float)
    out[: arr.size] = arr
    return out


_XI_CACHE: dict[tuple, tuple[float, float, float, float, float]] = {}
_PREFILLED_CG: set[float] = set()
_ANGLE_CACHE: dict[tuple, GreatCircleAngles] = {}
_CG_KEY_DIGITS: Final[int] = 12
_ZERO_XI: Final[tuple[float, float, float, float, float]] = (0.0, 0.0, 0.0, 0.0, 0.0)


def _angles_cached(n1, n2) -> GreatCircleAngles:
    """Memoized ``great_circle_angles`` keyed on the (n1, n2) component tuples.

    Profiling showed ~67% of a C-table cell's time was burned re-computing
    angles that don't change within a single ``dblquad``. Caching here is the
    single biggest win for sft-wick C-table builds.
    """
    n1a = _to_3d(n1)
    n2a = _to_3d(n2)
    key = (n1a[0], n1a[1], n1a[2], n2a[0], n2a[1], n2a[2])
    cached = _ANGLE_CACHE.get(key)
    if cached is not None:
        return cached
    angles = great_circle_angles(n1a, n2a)
    _ANGLE_CACHE[key] = angles
    return angles


def _prefill_all_shell_pairs(cg_key: float) -> None:
    """For a new cos_gamma, populate the cache for ALL (i1, i2) shell pairs.

    sft-wick's dblquad sweeps (lambda_1, lambda_2) across shell boundaries,
    so within ONE C-table cell at fixed cos_gamma the integrand visits many
    distinct (i1, i2) keys. Without prefill each new (i1, i2) hit triggers a
    ~3 ms cold compute; with prefill the entire cell is warm cache hits
    (~10 us each).
    """
    grid = _load_grid()
    if grid.get("is_trivial_zero", False):
        # Empty ell axis: every ell-sum is zero, so the per-(i1, i2) cache
        # has nothing to learn from this cos_gamma. _xi_at handles missing
        # cache entries by returning _ZERO_XI directly.
        _PREFILLED_CG.add(cg_key)
        return
    ells = grid["ells"]
    n_chi = grid["lambda_shells"].size
    weight = (2.0 * ells.astype(float) + 1.0) / (4.0 * np.pi)
    from scipy.special import eval_legendre

    cg_clip = max(-1.0, min(1.0, float(cg_key)))
    # Integer order is critical: scipy.special.eval_legendre with float order
    # returns NaN at x = -1 for orders >= ~1037 (generic Legendre evaluation
    # overflows). The integer-order path uses a stable recurrence.
    Pl = eval_legendre(ells, cg_clip)

    # Vectorise the Wigner-d evaluations across ells once per cos_gamma; reuse
    # for every (i1, i2) shell pair.
    from spin_rotation import _wigner_d_array
    d22 = _wigner_d_array(ells, 2, 2, cg_clip)
    d2m2 = _wigner_d_array(ells, 2, -2, cg_clip)
    d02 = _wigner_d_array(ells, 0, 2, cg_clip)
    w_d22 = weight * d22
    w_d2m2 = weight * d2m2
    w_d02 = weight * d02
    w_Pl = weight * Pl

    Cl_TT_grid = grid["Cl_phi00_phi00"]   # (N_l, N_chi, N_chi)
    Cl_TE_grid = grid["Cl_phi00_psi0"]    # (N_l, N_chi, N_chi)
    Cl_EE_grid = grid["Cl_psi0_psi0"]     # (N_l, N_chi, N_chi)

    for i1 in range(n_chi):
        for i2 in range(n_chi):
            key = (cg_key, i1, i2)
            if key in _XI_CACHE:
                continue
            Cl_TT = Cl_TT_grid[:, i1, i2]
            Cl_EE = Cl_EE_grid[:, i1, i2]
            Cl_TE = Cl_TE_grid[:, i1, i2]
            Cl_TE_swap = Cl_TE_grid[:, i2, i1]
            xi_TT = float(np.sum(w_Pl * Cl_TT))
            xi_p = float(np.sum(w_d22 * Cl_EE))
            xi_m = float(np.sum(w_d2m2 * Cl_EE))
            xi_TE_12 = float(np.sum(w_d02 * Cl_TE))
            xi_TE_21 = float(np.sum(w_d02 * Cl_TE_swap))
            _XI_CACHE[key] = (xi_TT, xi_p, xi_m, xi_TE_12, xi_TE_21)
    _PREFILLED_CG.add(cg_key)


def _xi_at(cg: float, i1: int, i2: int) -> tuple[float, float, float, float, float]:
    """Memoized 5-tuple ``(xi_TT, xi_+, xi_-, xi_TE_12, xi_TE_21)``.

    Triggers a one-shot prefill of ALL shell pairs the first time a new
    cos_gamma is seen, so subsequent (i1, i2) variations within the same
    dblquad cell are warm hits.
    """
    grid = _load_grid()
    if grid.get("is_trivial_zero", False):
        return _ZERO_XI
    cg_key = round(float(cg), _CG_KEY_DIGITS)
    if cg_key not in _PREFILLED_CG:
        _prefill_all_shell_pairs(cg_key)
    return _XI_CACHE[(cg_key, int(i1), int(i2))]


def kappa2(n1, t1: float, n2, t2: float) -> np.ndarray:
    """Bare source 2-pt cumulant, (3, 3) at one (n_1, t_1, n_2, t_2) point.

    Bilinearly interpolated in (lambda_1, lambda_2) across the 4 corner
    shell pairs at fixed cos_gamma. The smooth (vs nearest-neighbour)
    interpolation is essential for sft-wick's adaptive dblquad to converge
    in O(50) calls/cell instead of O(10^4-5) at every shell boundary.

    See module docstring for component layout.
    """
    grid = _load_grid()
    if grid.get("format") == "row_operator_v1":
        raise NotImplementedError(
            "kappa2_callable.kappa2 cannot evaluate point samples from a "
            "row-operator file. Use c_fast_gl.C_fn, which consumes row-bin "
            "integrals directly."
        )
    lambda_shells = grid["lambda_shells"]

    i1_lo, i1_hi, w1 = _bracket_shell(float(t1), lambda_shells)
    i2_lo, i2_hi, w2 = _bracket_shell(float(t2), lambda_shells)

    angles: GreatCircleAngles = _angles_cached(n1, n2)
    cg = angles.cos_gamma
    psi_1 = angles.psi_1
    psi_2 = angles.psi_2

    # Bilinear interpolation of the 5-tuple over the 4 corner shell pairs.
    xi_ll = _xi_at(cg, i1_lo, i2_lo)
    xi_lh = _xi_at(cg, i1_lo, i2_hi) if i2_hi != i2_lo else xi_ll
    xi_hl = _xi_at(cg, i1_hi, i2_lo) if i1_hi != i1_lo else xi_ll
    xi_hh = (_xi_at(cg, i1_hi, i2_hi)
             if (i1_hi != i1_lo or i2_hi != i2_lo) else xi_ll)

    a_ll = (1.0 - w1) * (1.0 - w2)
    a_lh = (1.0 - w1) * w2
    a_hl = w1 * (1.0 - w2)
    a_hh = w1 * w2

    xi_TT = a_ll * xi_ll[0] + a_lh * xi_lh[0] + a_hl * xi_hl[0] + a_hh * xi_hh[0]
    xi_p = a_ll * xi_ll[1] + a_lh * xi_lh[1] + a_hl * xi_hl[1] + a_hh * xi_hh[1]
    xi_m = a_ll * xi_ll[2] + a_lh * xi_lh[2] + a_hl * xi_hl[2] + a_hh * xi_hh[2]
    xi_TE_12 = a_ll * xi_ll[3] + a_lh * xi_lh[3] + a_hl * xi_hl[3] + a_hh * xi_hh[3]
    xi_TE_21 = a_ll * xi_ll[4] + a_lh * xi_lh[4] + a_hl * xi_hl[4] + a_hh * xi_hh[4]

    out = np.zeros((N_COMPONENTS, N_COMPONENTS), dtype=float)
    out[0, 0] = xi_TT
    qq_uu = spin2_real_block(xi_p, xi_m, psi_1, psi_2)
    out[1:3, 1:3] = qq_uu

    # Cross blocks: <T(1) Q(2)>, <T(1) U(2)>; only psi_2 enters since T(1)
    # is scalar. The transfer pipeline already produces Cl_phi00_psi0 in the
    # paper's new sign convention (Phi_00 = -(1/2) R k k), so no extra
    # sign factor is needed here.
    tq, tu = scalar_spin2_pair(xi_TE_12, psi_2)
    out[0, 1] = tq
    out[0, 2] = tu

    qt, ut = scalar_spin2_pair(xi_TE_21, psi_1)
    out[1, 0] = qt
    out[2, 0] = ut

    return out


def reset_xi_cache() -> None:
    """Drop the (cos_gamma, i1, i2) ell-sum cache (useful in tests)."""
    _XI_CACHE.clear()
    _PREFILLED_CG.clear()
    _ANGLE_CACHE.clear()


__all__ = ["N_COMPONENTS", "kappa2", "reset_xi_cache"]
