"""Fast vectorized C-propagator via 2D Gauss-Legendre quadrature.

Cache memory budget: ``_XI_TABLE_CACHE`` stores one ``(n_chi, n_chi, 5)``
tensor per distinct cos(gamma) value seen. Production sweeps visit ~40 cos
values so the cache stays under 100 MB at n_chi=200. Non-production callers
that walk cos(gamma) continuously (adaptive integrators, debug sweeps) are
protected by an LRU cap controlled via ``SFT_WICK_XI_TABLE_CACHE_MAX``
(default 128).


Plugs into sft-wick's L2 YAML at ``propagators.c_closed_form_module``:

    propagators:
      c_closed_form_module: c_fast_gl.py
      c_closed_form_attr: C_fn

Replaces sft-wick's internal ``dblquad``-driven C-table build (which calls
the integrand ~10000 times per cell via nested adaptive quadrature) with
a fixed-order 2D Gauss-Legendre rule on every cell — fully vectorised.
Per-cell cost drops from ~180 ms (post our cache+bilinear optimisations)
to ~200 us, a ~1000x reduction.

The integrand is

    f(lam1, lam2) = R(t1, lam1) * kappa2(n1, lam1, n2, lam2) * R(t2, lam2)

with the saddle-point response

    R(t, lam) = (D(lam) / D(t)) ** 2

derived from ``2 theta^(sa) = (log D^2)'`` (T5 in
``scripts/derive_jacobi_flat_flrw.wl``). ``D(lambda) = a(lambda) chi(lambda)``
is the closed-form Jacobi amplitude on flat FLRW (T1, see appendix
``append: D equals a chi`` of the paper); it is provided by
``D_callable.D_at`` and uses the fiducial-cosmology lightcone splines
in ``cosmo.py``. Working with the closed form avoids the 1/lambda
spline-overshoot pathology of the naive ``cumulative_trapezoid`` of
``2 theta^(sa)``; ``D(0)=0`` is exact, so ``R(t, 0)=0`` is preserved
without any near-vertex clamp.

Output is the (3, 3) C matrix.

Reuse:
- ``kappa2_callable._angles_cached``: per-(n1, n2) angle cache.
- ``kappa2_callable._xi_at`` and ``_prefill_all_shell_pairs``: per
  ``(cos_gamma, i1, i2)`` shell-pair cache of the ell-summed 5-tuple
  ``(xi_TT, xi_+, xi_-, xi_TE_12, xi_TE_21)``.
- ``kappa2_callable._load_grid``, ``_to_3d``: shell grid + 2D->3D embed.
- ``D_callable.D_at``: closed-form ``D(lambda) = chi(lambda)/(1+z(lambda))``,
  backed by the fiducial cosmology's lightcone splines in ``cosmo.py``. No
  prebuilt npz needed.

White-noise (sigma2) handling:
- Kappa2-only comparisons pass no sigma2 table; ``run_FF_single.py`` then sets
  ``SFT_WICK_GL_INCLUDE_SIGMA2=0``.
- Scale-split runs pass a canoes ``lambda_sigma2_delta`` table; the same
  runner sets ``SFT_WICK_GL_INCLUDE_SIGMA2=1`` and this module adds the
  collapsed one-dimensional high-ell contribution.

Tunables (env vars):
- ``SFT_WICK_GL_N`` (default 32): Gauss-Legendre order. Increase for
  smoother integrands; 32 gives ~10-12 digits for our smooth radial
  integrand. Keep <= 64 for memory locality.
"""

from __future__ import annotations

import os
import sys
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Final

import numpy as np

_HERE = Path(__file__).resolve().parent
_PIPELINE_DIR = _HERE.parent
for _p in (_HERE, _PIPELINE_DIR):
    if str(_p) not in sys.path:
        sys.path.append(str(_p))

# Reuse machinery from kappa2_callable.
from kappa2_callable import (  # noqa: E402
    _angles_cached, _bracket_shell, _CG_KEY_DIGITS, _load_grid,
    _prefill_all_shell_pairs, _XI_CACHE,
)

# Lazy sigma2 import: optional compatibility path from the original helper.
# The canoes pipeline normally runs with SFT_WICK_GL_INCLUDE_SIGMA2=0 and does
# not ship sigma2_callable.py, so this degrades gracefully to kappa2-only.
try:
    from sigma2_callable import sigma2 as _sigma2_fn  # noqa: E402
    try:
        from sigma2_callable import sigma2_batch as _sigma2_batch_fn  # noqa: E402
    except ImportError:  # pragma: no cover
        _sigma2_batch_fn = None
    _SIGMA2_AVAILABLE: Final[bool] = True
except (ImportError, FileNotFoundError):  # pragma: no cover
    _sigma2_fn = None
    _sigma2_batch_fn = None
    _SIGMA2_AVAILABLE = False

# User-facing toggle: ``SFT_WICK_GL_INCLUDE_SIGMA2=0`` disables the
# white-noise 1-D piece (e.g. for kappa2-only comparison runs).
_INCLUDE_SIGMA2: Final[bool] = bool(
    int(os.environ.get("SFT_WICK_GL_INCLUDE_SIGMA2", "0"))
)

N_COMPONENTS: Final[int] = 3
_GL_N: Final[int] = int(os.environ.get("SFT_WICK_GL_N", "32"))
_ROW_PARTIAL_POLICY_ENV: Final[str] = "SFT_WICK_ROW_OPERATOR_PARTIAL_BIN_POLICY"
_ROW_PARTIAL_TOL_ENV: Final[str] = "SFT_WICK_ROW_OPERATOR_PARTIAL_BIN_TOL"

# Closed-form D(lambda) = a(lambda) chi(lambda) is provided by D_callable.
# We delegate so that the Jacobi amplitude has a single source of truth
# across the package; the underlying lightcone splines are
# process-cached on the fiducial Cosmology (see cosmo.py).
from D_callable import D_at as _D_at  # noqa: E402

# --- Module-load: precompute GL nodes/weights. ---

_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(_GL_N)  # on [-1, 1]
_NODES_01 = 0.5 * (_NODES + 1.0)                            # on [0, 1]
_WEIGHTS_01 = 0.5 * _WEIGHTS                                # on [0, 1]


def _R_arr(lam_arr: np.ndarray, t: float) -> np.ndarray:
    """Saddle-point response R(t, lam) = (D(lam) / D(t))**2 (vectorised).

    Uses the closed form D(lambda) = a chi from D_callable.D_at; no
    Heaviside applied here (callers in this module integrate with lam in
    [0, t] only, so causality is automatic).
    """
    D_lam = _D_at(np.asarray(lam_arr, dtype=float))
    D_t = float(_D_at(float(t)))
    if D_t == 0.0:
        return np.zeros_like(D_lam)
    ratio = D_lam / D_t
    return ratio * ratio


# --- xi-table cache: assemble (n_chi, n_chi, 5) tensor per cos_gamma. ---

_XI_TABLE_CACHE_MAX: Final[int] = max(
    1, int(os.environ.get("SFT_WICK_XI_TABLE_CACHE_MAX", "128"))
)
_XI_TABLE_CACHE: "OrderedDict[float, np.ndarray]" = OrderedDict()
_RADIAL_KERNEL_WARNING_EMITTED = False
_ROW_DOMAIN_WARNING_EMITTED = False
# Lazily-allocated shared zero table for trivial-zero kappa2 stubs.
# A single (n_chi, n_chi, 5) tensor is returned for every cos_gamma so the
# cache footprint stays O(n_chi^2) instead of O(n_cg * n_chi^2).
_XI_TABLE_ZERO: np.ndarray | None = None


def _warn_if_radial_kernel_not_binned(grid: dict) -> None:
    """Warn once when the loaded kappa2 grid is not a quadrature-density kernel."""
    global _RADIAL_KERNEL_WARNING_EMITTED
    if _RADIAL_KERNEL_WARNING_EMITTED:
        return
    mode = str(grid.get("radial_kernel_convention", "unknown"))
    if mode in {
        "lambda_delta_density",
        "lambda_delta_bin_averaged",
        "lambda_cell_average_density",
        "lambda_bin_averaged",
        "lambda_binned_density",
    }:
        _RADIAL_KERNEL_WARNING_EMITTED = True
        return
    if grid.get("scale_split_role") == "exact_non_limber_low_ell_kappa2":
        _RADIAL_KERNEL_WARNING_EMITTED = True
        return
    warnings.warn(
        "Loaded kappa2 grid has radial_kernel_convention="
        f"{mode!r}. c_fast_gl integrates shell values with radial quadrature "
        "weights, which is only reliable for delta/bin-averaged kernels or "
        "for sample kernels whose radial peak is resolved by the native shell "
        "grid. Exact sample kernels should pass a row-area validation against "
        "the same trapezoid weights before production FF use.",
        RuntimeWarning,
        stacklevel=2,
    )
    _RADIAL_KERNEL_WARNING_EMITTED = True


def _xi_table_at(cg_key: float) -> np.ndarray:
    """Get the (n_chi, n_chi, 5) xi tensor at this cos_gamma (cached).

    Triggers ``_prefill_all_shell_pairs`` if needed, then assembles the
    cached scalars into a single ndarray for vectorised bilinear interp.

    For trivial-zero kappa2 stubs (L_cut=0 ``ell_window='low_empty'``),
    a single shared zero tensor is returned for every cos_gamma; the
    cache stays O(n_chi^2) total instead of O(n_cg * n_chi^2).
    """
    grid = _load_grid()
    if grid.get("is_trivial_zero", False):
        global _XI_TABLE_ZERO
        n_chi = grid["lambda_shells"].size
        if _XI_TABLE_ZERO is None or _XI_TABLE_ZERO.shape != (n_chi, n_chi, 5):
            _XI_TABLE_ZERO = np.zeros((n_chi, n_chi, 5), dtype=float)
            _XI_TABLE_ZERO.flags.writeable = False
        return _XI_TABLE_ZERO
    cached = _XI_TABLE_CACHE.get(cg_key)
    if cached is not None:
        # LRU bookkeeping: re-insert at the most-recently-used end.
        _XI_TABLE_CACHE.move_to_end(cg_key)
        return cached
    n_chi = grid["lambda_shells"].size
    _prefill_all_shell_pairs(cg_key)
    tab = np.empty((n_chi, n_chi, 5), dtype=float)
    for i1 in range(n_chi):
        for i2 in range(n_chi):
            tab[i1, i2] = _XI_CACHE[(cg_key, i1, i2)]
    _XI_TABLE_CACHE[cg_key] = tab
    # LRU eviction: drop the least-recently-used entry when over the cap.
    while len(_XI_TABLE_CACHE) > _XI_TABLE_CACHE_MAX:
        _XI_TABLE_CACHE.popitem(last=False)
    return tab


# --- Vectorised shell bracket. ---


def _bracket_arr(
    lam_arr: np.ndarray, lambda_shells: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised version of kappa2_callable._bracket_shell over a 1-D array.

    Returns ``(i_lo, i_hi, w_hi)`` arrays of shape ``lam_arr.shape``, with
    semantics: interpolated value = (1 - w_hi) * f[i_lo] + w_hi * f[i_hi].
    Outside the grid, clamped to the boundary (w_hi = 0).
    """
    n = lambda_shells.size
    lam = np.asarray(lam_arr, dtype=float)

    i_hi = np.searchsorted(lambda_shells, lam, side="right")
    i_hi = np.clip(i_hi, 0, n - 1)
    i_lo = np.maximum(i_hi - 1, 0)

    out_below = lam <= lambda_shells[0]
    out_above = lam >= lambda_shells[-1]
    in_range = ~(out_below | out_above)

    # Clamp i_hi/i_lo at the boundaries.
    i_hi = np.where(out_below, 0, i_hi)
    i_lo = np.where(out_above, n - 1, i_lo)
    i_hi = np.where(out_above, n - 1, i_hi)

    span = lambda_shells[i_hi] - lambda_shells[i_lo]
    safe_span = np.where(span > 0, span, 1.0)
    w_hi = np.where(
        in_range,
        (lam - lambda_shells[i_lo]) / safe_span,
        0.0,
    )
    return i_lo, i_hi, w_hi


# --- Vectorised spin-2 block constructors. ---


def _spin2_real_block_grid(
    xi_p: np.ndarray, xi_m: np.ndarray, psi_1: float, psi_2: float
) -> np.ndarray:
    """Vectorised ``spin_rotation.spin2_real_block``.

    Returns ``out[..., a, b]`` with ``a, b in {0=Q, 1=U}``, where
    ``...`` is the broadcast shape of ``xi_p``, ``xi_m``.
    """
    cp_diff = np.cos(2.0 * (psi_1 - psi_2))
    cp_sum = np.cos(2.0 * (psi_1 + psi_2))
    sp_diff = np.sin(2.0 * (psi_1 - psi_2))
    sp_sum = np.sin(2.0 * (psi_1 + psi_2))
    QQ = 0.5 * (cp_diff * xi_p + cp_sum * xi_m)
    UU = 0.5 * (cp_diff * xi_p - cp_sum * xi_m)
    QU = 0.5 * (sp_sum * xi_m - sp_diff * xi_p)
    UQ = 0.5 * (sp_sum * xi_m + sp_diff * xi_p)
    block = np.stack(
        [np.stack([QQ, QU], axis=-1), np.stack([UQ, UU], axis=-1)],
        axis=-2,
    )  # (..., 2, 2)
    return block


def _scalar_spin2_pair_grid(
    xi_TE: np.ndarray, psi_other: float
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised ``spin_rotation.scalar_spin2_pair``.

    Returns (Q-component, U-component), each broadcast-shaped like ``xi_TE``.
    """
    return (np.cos(2.0 * psi_other) * xi_TE, np.sin(2.0 * psi_other) * xi_TE)


# --- Main entry point. ---


def _trap_weights(lam: np.ndarray) -> np.ndarray:
    """1D trapezoidal weights on a (possibly non-uniform) ascending grid."""
    n = lam.size
    if n < 2:
        return np.zeros_like(lam)
    dl = np.diff(lam)
    w = np.empty(n)
    w[0] = 0.5 * dl[0]
    w[-1] = 0.5 * dl[-1]
    if n > 2:
        w[1:-1] = 0.5 * (dl[:-1] + dl[1:])
    return w


def _row_partial_bin_policy() -> str:
    policy = os.environ.get(_ROW_PARTIAL_POLICY_ENV, "fractional_width")
    policy = policy.strip().lower().replace("-", "_")
    if policy in {"fractional_width", "fractional", "width_fraction"}:
        return "fractional_width"
    if policy in {"error_if_large", "error"}:
        return "error_if_large"
    raise ValueError(
        f"unsupported {_ROW_PARTIAL_POLICY_ENV}={policy!r}; expected "
        "'fractional_width' or 'error_if_large'"
    )


def _row_partial_bin_tolerance() -> float:
    raw = os.environ.get(_ROW_PARTIAL_TOL_ENV, "1e-3")
    tol = float(raw)
    if tol < 0.0:
        raise ValueError(f"{_ROW_PARTIAL_TOL_ENV} must be non-negative")
    return tol


def _warn_if_row_operator_domain_truncated(grid: dict, t1: float, t2: float) -> None:
    """Warn once if source times extend beyond the stored row-operator domain."""
    global _ROW_DOMAIN_WARNING_EMITTED
    if _ROW_DOMAIN_WARNING_EMITTED:
        return
    meta = dict(grid.get("cosmo_meta", {}))
    if bool(meta.get("source_domain_intentionally_truncated", False)):
        _ROW_DOMAIN_WARNING_EMITTED = True
        return
    row_edges = np.asarray(grid.get("row_edges", []), dtype=float)
    if row_edges.size < 2:
        return
    lambda_max = float(row_edges[-1])
    row_width = float(np.nanmax(np.diff(row_edges)))
    source_max = max(float(t1), float(t2))
    if source_max > lambda_max + row_width:
        warnings.warn(
            "Row-operator source time exceeds the stored lambda domain: "
            f"source_max={source_max:.6g} Mpc, domain_max={lambda_max:.6g} Mpc. "
            "Rows outside the file are treated as zero; rebuild the row operator "
            "with a larger source domain for production use.",
            RuntimeWarning,
            stacklevel=3,
        )
        _ROW_DOMAIN_WARNING_EMITTED = True


def _dynamic_trapezoid_row_weights(
    lambda_rows: np.ndarray,
    t: float,
) -> np.ndarray | None:
    """Legacy dense-compatible moving-boundary trapezoid weights."""
    active = np.asarray(lambda_rows, dtype=float) <= float(t)
    if int(np.sum(active)) < 2:
        return None
    out = np.zeros(lambda_rows.shape, dtype=float)
    out[np.flatnonzero(active)] = _trap_weights(lambda_rows[active])
    return out


def _C_fn_row_operator(grid: dict, n1, t1: float, n2, t2: float) -> np.ndarray:
    """Consume a canoes row-bin operator without reconstructing dense K."""
    t1f = float(t1)
    t2f = float(t2)
    if t1f <= 0.0 or t2f <= 0.0:
        return np.zeros((N_COMPONENTS, N_COMPONENTS), dtype=float)
    _warn_if_row_operator_domain_truncated(grid, t1f, t2f)

    angles = _angles_cached(n1, n2)
    cg_clip = max(-1.0, min(1.0, float(angles.cos_gamma)))
    psi_1 = angles.psi_1
    psi_2 = angles.psi_2

    meta = dict(grid.get("cosmo_meta", {}))
    row_weight_kind = str(meta.get("row_weight_kind", "explicit_lambda_cell_widths"))
    lambda_rows_all = np.asarray(grid["lambda_rows"], dtype=float)
    row_idx = np.asarray(grid["op_bin_row_index"], dtype=np.intp)
    ell_idx = np.asarray(grid["op_bin_ell_index"], dtype=np.intp)
    row_lambda = lambda_rows_all[row_idx]
    stored_row_weights = np.asarray(grid["row_weights"], dtype=float)
    bin_left = np.asarray(grid["op_bin_lambda_left"], dtype=float)
    bin_right = np.asarray(grid["op_bin_lambda_right"], dtype=float)
    bin_eff = np.asarray(grid["op_bin_lambda_eff"], dtype=float)
    width = np.maximum(bin_right - bin_left, 0.0)

    values_all = np.asarray(grid["op_bin_integral"], dtype=float)
    dense_compat = (
        row_weight_kind == "trapezoid_on_lambda_rows"
        and "op_bin_target_index" in grid
    )
    if dense_compat:
        outer_weights = _dynamic_trapezoid_row_weights(lambda_rows_all, t1f)
        inner_weights = _dynamic_trapezoid_row_weights(lambda_rows_all, t2f)
        if outer_weights is None or inner_weights is None:
            return np.zeros((N_COMPONENTS, N_COMPONENTS), dtype=float)
        target_idx = np.asarray(grid["op_bin_target_index"], dtype=np.intp)
        base_inner = stored_row_weights[target_idx]
        inner_fraction = np.divide(
            inner_weights[target_idx],
            base_inner,
            out=np.zeros_like(base_inner, dtype=float),
            where=base_inner > 0.0,
        )
        row_weight_all = outer_weights[row_idx]
        bin_eff_all = lambda_rows_all[target_idx]
        active = (
            (row_weight_all > 0.0)
            & (inner_fraction > 0.0)
            & np.isfinite(inner_fraction)
        )
    else:
        row_edges = np.asarray(grid.get("row_edges", []), dtype=float)
        if (
            row_weight_kind == "explicit_lambda_cell_widths"
            and row_edges.shape == (lambda_rows_all.size + 1,)
        ):
            row_left = row_edges[:-1][row_idx]
            row_right = row_edges[1:][row_idx]
            row_width = np.maximum(row_right - row_left, 0.0)
            row_upper = np.minimum(row_right, t1f)
            outer_fraction = np.divide(
                np.maximum(row_upper - row_left, 0.0),
                row_width,
                out=np.zeros_like(row_width),
                where=row_width > 0.0,
            )
            row_weight_all = stored_row_weights[row_idx] * outer_fraction
        else:
            row_weight_all = stored_row_weights[row_idx] * (row_lambda <= t1f)

        inner_upper = np.minimum(bin_right, t2f)
        inner_width = np.maximum(inner_upper - bin_left, 0.0)
        inner_fraction = np.divide(
            inner_width,
            width,
            out=np.zeros_like(inner_width),
            where=width > 0.0,
        )
        bin_eff_all = np.where(
            inner_fraction >= 1.0,
            bin_eff,
            0.5 * (bin_left + inner_upper),
        )
        active = (
            (row_weight_all > 0.0)
            & (inner_fraction > 0.0)
            & np.isfinite(inner_fraction)
        )
    if not np.any(active):
        return np.zeros((N_COMPONENTS, N_COMPONENTS), dtype=float)

    ell_idx_a = ell_idx[active]
    row_lambda_a = row_lambda[active]
    row_weight_a = row_weight_all[active]
    fraction_a = inner_fraction[active]
    bin_eff_a = bin_eff_all[active]
    values = values_all[active]
    ells = np.asarray(grid["ells"], dtype=int)[ell_idx_a]

    from scipy.special import eval_legendre
    from spin_rotation import _wigner_d_array

    pref = (2.0 * ells.astype(float) + 1.0) / (4.0 * np.pi)
    w_Pl = pref * eval_legendre(ells, cg_clip)
    w_d22 = pref * _wigner_d_array(ells, 2, 2, cg_clip)
    w_d2m2 = pref * _wigner_d_array(ells, 2, -2, cg_clip)
    w_d02 = pref * _wigner_d_array(ells, 0, 2, cg_clip)

    R1 = _R_arr(row_lambda_a, t1f)
    R2 = _R_arr(bin_eff_a, t2f)
    radial = row_weight_a * R1 * fraction_a * R2
    crossing = fraction_a < 1.0
    if np.any(crossing) and not dense_compat:
        policy = _row_partial_bin_policy()
        if policy == "error_if_large":
            abs_values = np.asarray(
                grid.get("op_bin_abs_integral_estimate", np.abs(grid["op_bin_integral"])),
                dtype=float,
            )[active]
            angular_abs = np.column_stack([
                np.abs(w_Pl),
                np.abs(w_d02),
                np.abs(w_d02),
                np.maximum(np.abs(w_d22), np.abs(w_d2m2)),
            ])
            contribution_abs = (
                np.abs(row_weight_a * R1 * R2)[:, None]
                * angular_abs
                * abs_values
            )
            total_abs = float(np.sum(contribution_abs))
            crossing_abs = float(np.sum(contribution_abs[crossing]))
            tolerance = _row_partial_bin_tolerance()
            if crossing_abs > tolerance * max(total_abs, 1e-300):
                raise RuntimeError(
                    "row-operator inner source-time boundary crosses bins whose "
                    "estimated contribution is too large for fractional clipping: "
                    f"crossing_abs={crossing_abs:.6e}, total_abs={total_abs:.6e}, "
                    f"fraction={crossing_abs / max(total_abs, 1e-300):.6e}, "
                    f"tolerance={tolerance:.6e}. Rebuild with finer/CDF bins or "
                    f"set {_ROW_PARTIAL_POLICY_ENV}=fractional_width for a "
                    "diagnostic approximation."
                )

    xi_TT = float(np.sum(radial * w_Pl * values[:, 0]))
    xi_p = float(np.sum(radial * w_d22 * values[:, 3]))
    xi_m = float(np.sum(radial * w_d2m2 * values[:, 3]))
    xi_TE_12 = float(np.sum(radial * w_d02 * values[:, 1]))
    xi_TE_21 = float(np.sum(radial * w_d02 * values[:, 2]))

    C = np.zeros((N_COMPONENTS, N_COMPONENTS), dtype=float)
    C[0, 0] = xi_TT
    C[1:3, 1:3] = _spin2_real_block_grid(
        np.asarray(xi_p), np.asarray(xi_m), psi_1, psi_2
    )
    tq, tu = _scalar_spin2_pair_grid(np.asarray(xi_TE_12), psi_2)
    C[0, 1] = float(tq)
    C[0, 2] = float(tu)
    qt, ut = _scalar_spin2_pair_grid(np.asarray(xi_TE_21), psi_1)
    C[1, 0] = float(qt)
    C[2, 0] = float(ut)
    return C


def _add_sigma2_contribution(C: np.ndarray, n1, t1: float, n2, t2: float) -> np.ndarray:
    """Add the optional high-ell local sigma2 contribution.

    When ``sigma2_callable.sigma2_batch`` is available the 32-node Python
    loop over ``tau`` collapses to a single batched bilinear interp +
    spin-2 rotation, dropping the per-``C_fn`` cost from ~6 ms to <0.5 ms
    in L_cut=0 production sweeps.
    """
    if not (_INCLUDE_SIGMA2 and _SIGMA2_AVAILABLE and _sigma2_fn is not None):
        return C
    t_upper = min(float(t1), float(t2))
    if t_upper <= 0.0:
        return C
    tau = t_upper * _NODES_01
    w_tau = t_upper * _WEIGHTS_01
    R1_tau = _R_arr(tau, float(t1))
    R2_tau = _R_arr(tau, float(t2))
    if _sigma2_batch_fn is not None:
        sig_arr = _sigma2_batch_fn(n1, tau, n2)
    else:
        sig_arr = np.empty((tau.size, N_COMPONENTS, N_COMPONENTS), dtype=float)
        for i, tau_i in enumerate(tau):
            sig_arr[i] = _sigma2_fn(n1, float(tau_i), n2)
    wRR = w_tau * R1_tau * R2_tau
    C += np.einsum("i,iab->ab", wRR, sig_arr)
    return C


def C_fn(n1, t1: float, n2, t2: float) -> np.ndarray:
    """Return the (3, 3) C-propagator at (n1, t1, n2, t2).

    Uses 2D *trapezoidal* quadrature directly on the kappa2 shell grid
    (no bilinear interpolation). This is exact at the shell resolution
    and avoids two failure modes of the previous fixed-order GL path:

    1. **Diagonal smearing**: kappa2(lambda1, lambda2) is sharply peaked
       on the diagonal with width ~ chi/ell_max ~ 2 Mpc. With GL_N=32
       on a 1839 Mpc domain (node spacing ~60 Mpc), bilinear interp
       between dense (5 Mpc) shells could not resolve the peak.

    2. **Off-diagonal Bessel-Bessel noise aliasing**: the full
       Bessel-Bessel C_ell(chi_a, chi_b) at high ell has noisy
       off-diagonal residuals (~1e-3 of the diagonal) from rapidly
       oscillating Bessel products. Coarse GL nodes alias these into
       spurious O(10x) contributions; trapezoidal on the dense shell
       grid integrates them to ~zero (their true Limber value).

    Cost: O(n_shells**2) per cell. For n_shells=470 (dlambda=5 Mpc) and
    8 sweep cells, ~50ms / cell, well below the per-cell sft-wick
    overhead.
    """
    t1f = float(t1)
    t2f = float(t2)
    if t1f <= 0.0 or t2f <= 0.0:
        return np.zeros((N_COMPONENTS, N_COMPONENTS), dtype=float)

    grid = _load_grid()
    if grid.get("format") == "row_operator_v1":
        C = _C_fn_row_operator(grid, n1, t1f, n2, t2f)
        return _add_sigma2_contribution(C, n1, t1f, n2, t2f)

    _warn_if_radial_kernel_not_binned(grid)
    lambda_shells = grid["lambda_shells"]

    # Angles (cached on (n1, n2)).
    angles = _angles_cached(n1, n2)
    cg = angles.cos_gamma
    psi_1 = angles.psi_1
    psi_2 = angles.psi_2
    cg_key = round(float(cg), _CG_KEY_DIGITS)

    # Integration grid: ALL kappa2 shells in [shell_min, t]. We integrate
    # on the native shell positions to avoid bilinear interpolation; the
    # at-most-dlambda truncation between shell_max(<=t) and t is dropped
    # (negligible for dlambda ~ 5 Mpc).
    mask1 = lambda_shells <= t1f
    mask2 = lambda_shells <= t2f
    if mask1.sum() < 2 or mask2.sum() < 2:
        return np.zeros((N_COMPONENTS, N_COMPONENTS), dtype=float)
    lam1 = lambda_shells[mask1]                   # (n1,)
    lam2 = lambda_shells[mask2]                   # (n2,)
    idx1 = np.flatnonzero(mask1)                  # indices into lambda_shells
    idx2 = np.flatnonzero(mask2)
    w1_lam = _trap_weights(lam1)                  # (n1,)
    w2_lam = _trap_weights(lam2)                  # (n2,)

    # Response on the integration grid: R(t, lam) = (D(lam) / D(t))**2.
    R1 = _R_arr(lam1, t1f)  # (n1,)
    R2 = _R_arr(lam2, t2f)  # (n2,)

    # xi table at this cos_gamma: shape (n_chi, n_chi, 5). Slice to the
    # active shell sub-grid -- direct lookup, no interpolation.
    xi_tab = _xi_table_at(cg_key)
    xi_grid = xi_tab[np.ix_(idx1, idx2)]          # (n1, n2, 5)

    # Unpack the 5 channels.
    xi_TT = xi_grid[..., 0]      # (n1, n2)
    xi_p = xi_grid[..., 1]
    xi_m = xi_grid[..., 2]
    xi_TE_12 = xi_grid[..., 3]
    xi_TE_21 = xi_grid[..., 4]

    # Build the (N, N, 3, 3) kappa2 grid.
    kappa = np.zeros(xi_TT.shape + (N_COMPONENTS, N_COMPONENTS), dtype=float)
    kappa[..., 0, 0] = xi_TT
    kappa[..., 1:3, 1:3] = _spin2_real_block_grid(xi_p, xi_m, psi_1, psi_2)
    tq, tu = _scalar_spin2_pair_grid(xi_TE_12, psi_2)
    kappa[..., 0, 1] = tq
    kappa[..., 0, 2] = tu
    qt, ut = _scalar_spin2_pair_grid(xi_TE_21, psi_1)
    kappa[..., 1, 0] = qt
    kappa[..., 2, 0] = ut

    # Quadrature reduction: C[a, b] = sum_ij w1[i] R1[i] kappa[i,j,a,b] w2[j] R2[j].
    wR1 = w1_lam * R1  # (n1,)
    wR2 = w2_lam * R2  # (n2,)
    C = np.einsum("i,j,ijab->ab", wR1, wR2, kappa)

    # White-noise (δ-correlated) 1-D piece, added on top of the kappa2
    # 2-D integral. The δ collapses one of the time integrals so the
    # contribution is
    #     C_white[a, b] = ∫_0^{min(t1, t2)} R(t1, τ) σ²[a, b](n1, τ, n2) R(t2, τ) dτ
    # See ``evaluate.py:_C_value_direct`` (line ~1478) for the engine's
    # reference dblquad path; we duplicate it here on the same fixed-order
    # GL grid so c_closed_form_only=true users still get the high-ℓ
    # white-noise contribution that ell_cut truncates from kappa2.
    return _add_sigma2_contribution(C, n1, t1f, n2, t2f)


__all__ = ["N_COMPONENTS", "C_fn"]
