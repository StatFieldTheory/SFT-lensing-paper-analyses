"""Candidate repairs for the non-PSD Sigma2 density, as drop-in subclasses.

ROOT CAUSE (established in d1-d4).  ``corr_op.C_fn`` interpolates its table
between 21 lambda nodes and is only C^1 there: its second derivative is
piecewise constant and jumps at every node (factor 1.55 at lambda = 1300).
``Sigma2Builder`` fits a globally C^2 interpolating cubic spline through 160
samples of that function and differentiates it.  A C^2 spline cannot represent a
curvature jump, so it RINGS in the neighbouring cells; the ringing is measured at
+13, -22, +47, +7, -77, -132, -56, -2, +10 percent across the cells bracketing
lambda = 1300.  The -132% overshoot drives the smallest eigenvalue negative.
Eleven such windows exist, one 3-10 Mpc after each interior table node.

The true density is continuous with a SLOPE discontinuity at each table node
(``C`` is C^1, so ``dC/dlam`` is continuous).  Every repair below respects that.

R0  reference : per-cell high-order differentiation, knots at the 21 table nodes.
                Slow (many corr_op calls) but as close to exact as we can get;
                used as the yardstick, not as production.
R1  knots     : same interpolate-then-differentiate recipe, but with the spline
                broken at the 21 table nodes so no cell straddles a kink.
R2  density   : sample the DENSITY by central differences inside each cell, then
                interpolate the density itself (never differentiate an
                interpolant).
R3  smoothed  : keep the current recipe but use a penalised smoothing spline, so
                the ringing is damped rather than avoided.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline, PchipInterpolator, make_smoothing_spline

import driver_stats as ds
import corr_op_C_callable as _corr_op

N_COMP = ds.N_COMP
TABLE_LAM: NDArray[np.float64] = np.load(
    _corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)


def _cells(lam_lo: float, lam_hi: float) -> NDArray[np.float64]:
    """Interior corr_op table nodes inside the sampling span, as break points."""
    inner = TABLE_LAM[(TABLE_LAM > lam_lo + 1e-9) & (TABLE_LAM < lam_hi - 1e-9)]
    return np.concatenate([[lam_lo], inner, [lam_hi]])


class _Base(ds.Sigma2Builder):
    """Shares the constructor and the cache; only ``matrix`` changes."""

    def _sample_f(self, cos_gamma: float, lam: NDArray[np.float64]):
        """``D^4 C_ab(t,t)`` at the requested lambda."""
        n2 = ds._n2_of_cos(cos_gamma)
        C = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t)) for t in lam])
        return (np.asarray(self.bg.D(lam), float) ** 4)[:, None, None] * C

    def _finish(self, out: NDArray[np.float64]) -> NDArray[np.float64]:
        out = 0.5 * (out + out.T)
        if self.apply_c0:
            out = out * ds.ANCHOR_C0
        return out


class R1Knots(_Base):
    """Interpolate-then-differentiate, but never across a corr_op table node."""

    def _density_splines(self, cos_gamma: float):
        key = round(float(cos_gamma), 12)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        breaks = _cells(self.lam_lo, self.lam_hi)
        per_cell = []
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            n = max(6, int(np.ceil((hi - lo) / 12.0)) + 1)
            xs = np.linspace(lo, hi, n)
            f = self._sample_f(cos_gamma, xs)
            per_cell.append((lo, hi, [[CubicSpline(xs, f[:, a, b])
                                       for b in range(N_COMP)]
                                      for a in range(N_COMP)]))
        self._cache[key] = per_cell
        return per_cell

    def matrix(self, cos_gamma: float, lam: float) -> NDArray[np.float64]:
        cells = self._density_splines(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        for lo, hi, sp in cells:
            if lo - 1e-9 <= x <= hi + 1e-9:
                break
        d4 = float(self.bg.D(x)) ** 4
        out = np.array([[sp[a][b](x, 1) for b in range(N_COMP)]
                        for a in range(N_COMP)]) / d4
        return self._finish(out)


class R2Density(_Base):
    """Sample the density by central differences, then interpolate the density."""

    def __init__(self, *a, fd_h: float = 2.0, n_per_cell: int = 13, **kw):
        super().__init__(*a, **kw)
        self.fd_h = float(fd_h)
        self.n_per_cell = int(n_per_cell)

    def _density_splines(self, cos_gamma: float):
        key = round(float(cos_gamma), 12)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        breaks = _cells(self.lam_lo, self.lam_hi)
        xs_all, dens_all = [], []
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            # keep the FD stencil strictly inside the cell so it never straddles a
            # kink; the corr_op grid has cells as narrow as 4 Mpc near the source,
            # so the step has to adapt to the cell.
            h = min(self.fd_h, (hi - lo) / 8.0)
            n = max(4, min(self.n_per_cell, int((hi - lo) / max(h, 1e-9))))
            xs = np.linspace(lo + h, hi - h, n)
            fp = self._sample_f(cos_gamma, xs + h)
            fm = self._sample_f(cos_gamma, xs - h)
            d4 = np.asarray(self.bg.D(xs), float) ** 4
            dens = (fp - fm) / (2.0 * h) / d4[:, None, None]
            # extend to the cell edges by linear extrapolation in lambda
            xs_all.append(np.concatenate([[lo], xs, [hi]]))
            e0 = dens[0] + (dens[1] - dens[0]) * (lo - xs[0]) / (xs[1] - xs[0])
            e1 = dens[-1] + (dens[-1] - dens[-2]) * (hi - xs[-1]) / (xs[-1] - xs[-2])
            dens_all.append(np.concatenate([e0[None], dens, e1[None]]))
        x = np.concatenate(xs_all)
        y = np.concatenate(dens_all)
        o = np.argsort(x, kind="stable")
        x, y = x[o], y[o]
        keep = np.concatenate([[True], np.diff(x) > 1e-9])
        x, y = x[keep], y[keep]
        # PCHIP: C^1, shape-preserving, no overshoot -- the density has kinks
        flat = [PchipInterpolator(x, y[:, a, b])
                for a in range(N_COMP) for b in range(N_COMP)]
        self._cache[key] = flat
        return flat

    def matrix(self, cos_gamma: float, lam: float) -> NDArray[np.float64]:
        flat = self._density_splines(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        out = np.array([[flat[a * N_COMP + b](x) for b in range(N_COMP)]
                        for a in range(N_COMP)])
        return self._finish(out)


class R3Smoothed(_Base):
    """Current recipe, but with a penalised smoothing spline before differentiating."""

    def _density_splines(self, cos_gamma: float):
        key = round(float(cos_gamma), 12)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        f = self._sample_f(cos_gamma, self._nodes)
        flat = []
        for a in range(N_COMP):
            for b in range(N_COMP):
                y = f[:, a, b]
                if np.allclose(y, 0.0):
                    flat.append(CubicSpline(self._nodes, y))
                else:
                    flat.append(make_smoothing_spline(self._nodes, y))
        self._cache[key] = flat
        return flat

    def matrix(self, cos_gamma: float, lam: float) -> NDArray[np.float64]:
        flat = self._density_splines(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        d4 = float(self.bg.D(x)) ** 4
        out = np.empty((N_COMP, N_COMP))
        for a in range(N_COMP):
            for b in range(N_COMP):
                s = flat[a * N_COMP + b]
                out[a, b] = (s(x, 1) if isinstance(s, CubicSpline)
                             else s.derivative()(x)) / d4
        return self._finish(out)


class R0Reference(_Base):
    """Yardstick: dense central differences inside the containing cell, no
    interpolation of the result at all.  Expensive; for validation only."""

    def __init__(self, *a, fd_h: float = 0.5, **kw):
        super().__init__(*a, **kw)
        self.fd_h = float(fd_h)

    def matrix(self, cos_gamma: float, lam: float) -> NDArray[np.float64]:
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        breaks = _cells(self.lam_lo, self.lam_hi)
        j = int(np.searchsorted(breaks, x) - 1)
        lo, hi = breaks[max(j, 0)], breaks[min(j + 1, len(breaks) - 1)]
        h = min(self.fd_h, max((hi - lo) / 8.0, 1e-3))
        xq = float(np.clip(x, lo + h, hi - h))
        fp = self._sample_f(cos_gamma, np.array([xq + h]))[0]
        fm = self._sample_f(cos_gamma, np.array([xq - h]))[0]
        d4 = float(self.bg.D(xq)) ** 4
        return self._finish((fp - fm) / (2.0 * h) / d4)
