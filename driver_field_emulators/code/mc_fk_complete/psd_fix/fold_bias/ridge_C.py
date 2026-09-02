"""Ridge-aware (lam1, lam2) readings of a corr_op table.

corr_op serves its 21x21 lambda table through a BILINEAR
`RegularGridInterpolator`.  C(lam1, lam2) has a ridge along the equal-time
diagonal (the line-of-sight integral runs to min(lam1, lam2)), so at a cell
midpoint the bilinear value averages two on-diagonal corners with two
off-diagonal ones and undershoots by 4-27%.

This module rebuilds the SAME angular channels canoes builds (so the ell-sum,
Wigner-d table, taper and psi rotation are byte-identical to production) and
then applies a different (lam1, lam2) reading of `channels.M_table`:

    "bilinear"  : production; RegularGridInterpolator(method='linear').
    "cubic"     : RegularGridInterpolator(method='cubic').
    "tri"       : Freudenthal ("s >= t") triangulation of the same product
                  grid.  The main diagonal is a triangle EDGE, so on the
                  diagonal the interpolant reduces to the 1-D linear reading
                  of the exactly-stored diagonal values.  Continuous
                  everywhere, exact at every node, no new data.
    "tri_log"   : same triangulation, interpolating log|C| with a common sign
                  (falls back to "tri" wherever the three vertices are not
                  strictly same-signed).
    "ridge"     : the recommended cure.  Factor out the exactly-stored
                  equal-time diagonal of the Phi00 channel,

                      C_ab(l1, l2) = D(m) * rho_ab(l1, l2),
                      m = (l1 + l2) / 2,  D = log-PCHIP through C_00(l_i, l_i),

                  and interpolate the residual rho -- which is O(1), smooth and
                  slowly varying, because the steep growth of C toward the
                  source plane has been divided out -- on the SAME product grid
                  with the diagonal-aligned triangulation.  Exact at every node
                  (rho_ij = M_ij / D(m_ij) by construction), exactly equal to
                  the 1-D log-PCHIP reading of the diagonal on l1 = l2, and
                  continuous everywhere.  No new table data.
    "diag_only" : diagnostic; every query is answered by the 1-D reading of
                  the diagonal at (lam1+lam2)/2.  Not a candidate cure.

All schemes keep corr_op's own out-of-table contract: zero when either source
lies at or below the table's smallest node.
"""
from __future__ import annotations

import sys
from typing import Literal

import numpy as np
from numpy.typing import NDArray

sys.path.insert(0, "/Users/zzhang/projects/angular_statistics/canoes/src")

from canoes.sachs.sft_input.corr_op.table import (  # noqa: E402
    _get_or_build_channels,
    _great_circle_angles,
    load_table_2d_from_npz,
)
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator  # noqa: E402

Scheme = Literal["bilinear", "cubic", "tri", "tri_log", "ridge", "diag_only"]


def _locate(grid: NDArray[np.float64], v: NDArray[np.float64]):
    """Cell index and local coordinate, clamped to the grid."""
    i = np.clip(np.searchsorted(grid, v) - 1, 0, grid.size - 2)
    h = grid[i + 1] - grid[i]
    return i, (v - grid[i]) / h


def _tri_weights(s, t):
    """Freudenthal split of the unit square along s = t.

    Returns the four corner weights (w00, w10, w01, w11).  On s == t the two
    triangles agree and the weights reduce to (1-s, 0, 0, s) -- a pure 1-D
    linear reading along the diagonal.
    """
    lower = s >= t                       # triangle (0,0) (1,0) (1,1)
    w00 = np.where(lower, 1.0 - s, 1.0 - t)
    w10 = np.where(lower, s - t, 0.0)
    w01 = np.where(lower, 0.0, t - s)
    w11 = np.where(lower, t, s)
    return w00, w10, w01, w11


def _tri_eval(M, grid, l1, l2, log_space: bool = False):
    """Triangular interpolation of M (..., N, N) at (l1, l2)."""
    i, s = _locate(grid, l1)
    j, t = _locate(grid, l2)
    w00, w10, w01, w11 = _tri_weights(s, t)
    f00 = M[..., i, j]
    f10 = M[..., i + 1, j]
    f01 = M[..., i, j + 1]
    f11 = M[..., i + 1, j + 1]
    lin = w00 * f00 + w10 * f10 + w01 * f01 + w11 * f11
    if not log_space:
        return lin
    # Geometric variant: only where every contributing corner is same-signed
    # and nonzero (weights can be exactly 0, in which case that corner is
    # irrelevant -- mask it in by copying the sign of the others).
    ws = np.stack([w00, w10, w01, w11], axis=0)
    fs = np.stack([f00, f10, f01, f11], axis=0)
    # values carry leading component axes; give the weights matching trailing-
    # axis broadcasting (weights vary only along the sample axis).
    ws = ws.reshape(ws.shape[0], *((1,) * (fs.ndim - ws.ndim)), *ws.shape[1:])
    active = ws > 0
    sgn = np.sign(fs)
    ok = np.all(np.where(active, sgn == sgn[0], True), axis=0) & (sgn[0] != 0)
    ok &= np.all(np.where(active, fs != 0.0, True), axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.sum(np.where(active, ws * np.log(np.abs(np.where(fs == 0, 1.0, fs))), 0.0), axis=0)
        geo = sgn[0] * np.exp(lg)
    return np.where(ok, geo, lin)


class RidgeC:
    """A corr_op C_fn with a selectable (lam1, lam2) reading."""

    def __init__(self, table_path: str, scheme: Scheme = "bilinear") -> None:
        self.table = load_table_2d_from_npz(str(table_path))
        self.scheme: Scheme = scheme
        self.grid = np.asarray(self.table.cl_table.lambda_grid, dtype=float)
        self.lam_min = float(self.grid.min())
        self._rgi: dict[tuple, RegularGridInterpolator] = {}
        self._rho: dict[tuple, tuple] = {}
        if str(self.table.config.source_coord) != "lambda_project_mpc":
            raise ValueError("expected a lambda_project_mpc table")

    # -- channel table (identical to production up to the lam reading) ------
    def _channels(self, n1, n2):
        cos_g, psi_1, psi_2 = _great_circle_angles(n1, n2)
        return _get_or_build_channels(self.table, cos_g, psi_1, psi_2,
                                      interp_method="linear")

    def _stacked_rgi(self, ch, method: str):
        key = (round(ch.cos_g, 10), round(ch.psi_1, 10), round(ch.psi_2, 10), method)
        r = self._rgi.get(key)
        if r is None:
            stacked = np.stack([ch.M_table[i, j] for i in range(3) for j in range(3)],
                               axis=-1)
            # scipy 1.17's cubic RegularGridInterpolator silently returns
            # EXACTLY ZERO for values of order 1e-10 (verified on a 6x6 toy
            # grid: the same data scaled by 1e-10 gives 0.0 instead of 7.25).
            # Rescale to O(1) before the spline solve and undo afterwards.
            scale = float(np.max(np.abs(stacked))) or 1.0
            r = RegularGridInterpolator((self.grid, self.grid), stacked / scale,
                                        method=method, bounds_error=False,
                                        fill_value=None)
            r._corr_op_scale = scale  # type: ignore[attr-defined]
            self._rgi[key] = r
        return r

    def _ridge_parts(self, ch):
        """(D_pchip, rho_table) for the ridge-factorised scheme."""
        key = (round(ch.cos_g, 10), round(ch.psi_1, 10), round(ch.psi_2, 10))
        got = self._rho.get(key)
        if got is not None:
            return got
        diag = np.diag(ch.M_table[0, 0]).astype(float)   # C_00(l_i, l_i) > 0
        if np.any(diag <= 0):
            raise ValueError("Phi00 equal-time diagonal is not strictly positive; "
                             "the ridge normaliser is undefined for this table.")
        logD = PchipInterpolator(self.grid, np.log(diag), extrapolate=True)

        def D(x):
            return np.exp(logD(np.clip(x, self.grid[0], self.grid[-1])))

        m_ij = 0.5 * (self.grid[:, None] + self.grid[None, :])
        rho = ch.M_table / D(m_ij)                        # (3, 3, N, N)
        self._rho[key] = (D, rho)
        return self._rho[key]

    # -- public ------------------------------------------------------------
    def C_fn(self, n1, t1, n2, t2) -> NDArray[np.float64]:
        out = self.C_fn_batch(n1, np.array([float(t1)]), n2, np.array([float(t2)]))
        return out[0]

    def C_fn_batch(self, n1_arr, t1_arr, n2_arr, t2_arr) -> NDArray[np.float64]:
        t1 = np.atleast_1d(np.asarray(t1_arr, dtype=float))
        t2 = np.atleast_1d(np.asarray(t2_arr, dtype=float))
        n = max(t1.size, t2.size)
        t1 = np.broadcast_to(t1, (n,))
        t2 = np.broadcast_to(t2, (n,))
        a1 = np.broadcast_to(np.atleast_2d(np.asarray(n1_arr, float)), (n, 3))
        a2 = np.broadcast_to(np.atleast_2d(np.asarray(n2_arr, float)), (n, 3))
        out = np.zeros((n, 3, 3), dtype=float)
        live = (t1 > 0) & (t2 > 0) & (t1 >= self.lam_min) & (t2 >= self.lam_min)
        if not np.any(live):
            return out
        keys: dict[tuple, list[int]] = {}
        for k in np.flatnonzero(live):
            cg, p1, p2 = _great_circle_angles(a1[k], a2[k])
            keys.setdefault((round(cg, 8), round(p1, 8), round(p2, 8)), []).append(int(k))
        for (cg, p1, p2), idxs in keys.items():
            ch = _get_or_build_channels(self.table, cg, p1, p2, interp_method="linear")
            idx = np.asarray(idxs, dtype=int)
            l1, l2 = t1[idx], t2[idx]
            if self.scheme in ("bilinear", "cubic"):
                m = "linear" if self.scheme == "bilinear" else "cubic"
                rgi = self._stacked_rgi(ch, m)
                vals = rgi(np.stack([l1, l2], axis=-1)) * rgi._corr_op_scale
                out[idx] = vals.reshape(-1, 3, 3)
            elif self.scheme in ("tri", "tri_log"):
                v = _tri_eval(ch.M_table, self.grid, l1, l2,
                              log_space=(self.scheme == "tri_log"))
                out[idx] = np.moveaxis(v, -1, 0)
            elif self.scheme == "ridge":
                D, rho = self._ridge_parts(ch)
                v = _tri_eval(rho, self.grid, l1, l2) * D(0.5 * (l1 + l2))
                out[idx] = np.moveaxis(v, -1, 0)
            elif self.scheme == "diag_only":
                m = 0.5 * (l1 + l2)
                v = _tri_eval(ch.M_table, self.grid, m, m)
                out[idx] = np.moveaxis(v, -1, 0)
            else:
                raise ValueError(f"unknown scheme {self.scheme!r}")
        return out
