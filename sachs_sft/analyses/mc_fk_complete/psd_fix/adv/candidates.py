"""Adversarial candidates: an EXACT closed-form density, a properly fixed R2,
and the two schemes the attack names (Akima; C^1 Hermite through the nodes).

Key structural fact (verified in a1_structure.py to 1.2e-14):
    corr_op interpolates its 21x21 table BILINEARLY, so restricted to the
    equal-time diagonal C(t,t) is EXACTLY a quadratic Bezier on each cell,
        C(u) = (1-u)^2 A + u(1-u) S + u^2 B,   u = (lam - l_i)/h,
        A = M[i,i],  B = M[i+1,i+1],  S = M[i,i+1] + M[i+1,i].
    Hence  dC/dlam  is available in CLOSED FORM and the density
        Sigma2 = dC/dlam + 4 (D'/D) C
    needs no differentiation of any interpolant at all.
"""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import PchipInterpolator, Akima1DInterpolator, CubicHermiteSpline

import driver_stats as ds
import corr_op_C_callable as _corr_op

N_COMP = ds.N_COMP
TABLE_LAM: NDArray[np.float64] = np.load(
    _corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)

# ---- instrumented corr_op counter -----------------------------------------
class _Counter:
    n = 0
COUNT = _Counter()

def _C(n1, t1, n2, t2):
    COUNT.n += 1
    return np.asarray(_corr_op.C_fn(n1, float(t1), n2, float(t2)), float)


def _cells(lam_lo: float, lam_hi: float) -> NDArray[np.float64]:
    inner = TABLE_LAM[(TABLE_LAM > lam_lo + 1e-9) & (TABLE_LAM < lam_hi - 1e-9)]
    return np.concatenate([[lam_lo], inner, [lam_hi]])


class _Base(ds.Sigma2Builder):
    def _sample_f(self, cos_gamma, lam):
        n2 = ds._n2_of_cos(cos_gamma)
        C = np.array([_C(ds._N1, float(t), n2, float(t)) for t in lam])
        return (np.asarray(self.bg.D(lam), float) ** 4)[:, None, None] * C

    def _finish(self, out):
        out = 0.5 * (out + out.T)
        if self.apply_c0:
            out = out * ds.ANCHOR_C0
        return out


# ---------------------------------------------------------------------------
class R5Exact(_Base):
    """Closed-form per-cell analytic density.  No interpolation of anything.

    Construction cost: 21 diagonal + 20 off-diagonal corr_op evaluations per
    cos_gamma (41), versus 160 for the stock builder.
    """

    def __init__(self, *a, use_transpose_symmetry: bool = True, **kw):
        super().__init__(*a, **kw)
        self.use_transpose_symmetry = bool(use_transpose_symmetry)
        # exact node values of D'/D from the background's own D spline
        self._lam_nodes = TABLE_LAM

    def _theta(self, x):
        # d ln D / dlam taken from the SAME D spline the production builder
        # differentiates through D^4 (NOT bg.theta_sa, which splines ln D
        # separately and differs from it by up to 1e-2 near the grid edge).
        sp = self.bg._D_spline
        return float(sp(x, 1) / sp(x))

    def verify_structure(self, cos_gamma: float, atol: float = 1e-10) -> float:
        """FAIL-LOUD guard on the one assumption R5Exact rests on.

        R5Exact is exact only while corr_op interpolates its table BILINEARLY
        (``interp_method='linear'``, the canoes default).  If the table is ever
        rebuilt or served with ``interp_method='cubic'`` the diagonal stops being
        a quadratic Bezier and this class silently becomes wrong.  The guard
        re-derives C at three interior abscissae per cell from the closed form
        and compares against the served C_fn.
        """
        L, diag, S = self._coeffs(cos_gamma)
        n2 = ds._n2_of_cos(cos_gamma)
        worst = 0.0
        for i in range(L.size - 1):
            h = L[i + 1] - L[i]
            for u in (0.25, 0.5, 0.75):
                pred = ((1 - u) ** 2 * diag[i] + u * (1 - u) * S[i]
                        + u ** 2 * diag[i + 1])
                act = _C(ds._N1, float(L[i] + u * h), n2, float(L[i] + u * h))
                sc = max(float(np.max(np.abs(act))), 1e-300)
                worst = max(worst, float(np.max(np.abs(pred - act))) / sc)
        if worst > atol:
            raise RuntimeError(
                f"R5Exact structural guard FAILED: max rel |bezier - C_fn| = "
                f"{worst:.3e} > {atol:.1e}.  corr_op is no longer bilinear on "
                f"its lambda table (interp_method changed?), so the closed-form "
                f"density is invalid.  Do not use R5Exact with this table.")
        return worst

    def _coeffs(self, cos_gamma):
        key = ("R5", round(float(cos_gamma), 12))
        c = self._cache.get(key)
        if c is not None:
            return c
        n2 = ds._n2_of_cos(cos_gamma)
        L = TABLE_LAM
        diag = np.array([_C(ds._N1, float(t), n2, float(t)) for t in L])
        S = np.empty((L.size - 1, N_COMP, N_COMP))
        for i in range(L.size - 1):
            up = _C(ds._N1, float(L[i]), n2, float(L[i + 1]))
            if self.use_transpose_symmetry:
                # C_ab(t1,t2) == C_ba(t2,t1); verified to 3.5e-16 in a2.
                S[i] = up + up.T
            else:
                S[i] = up + _C(ds._N1, float(L[i + 1]), n2, float(L[i]))
        out = (L, diag, S)
        self._cache[key] = out
        return out

    def matrix(self, cos_gamma, lam):
        L, diag, S = self._coeffs(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        i = int(np.clip(np.searchsorted(L, x, side="right") - 1, 0, L.size - 2))
        h = L[i + 1] - L[i]
        u = (x - L[i]) / h
        A, B, Si = diag[i], diag[i + 1], S[i]
        C = (1 - u) ** 2 * A + u * (1 - u) * Si + u ** 2 * B
        dC = (-2 * (1 - u) * A + (1 - 2 * u) * Si + 2 * u * B) / h
        return self._finish(dC + 4.0 * self._theta(x) * C)

    def matrix_side(self, cos_gamma, lam, side="right"):
        """Density at an exact table node, choosing which cell to use."""
        L, diag, S = self._coeffs(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        if side == "left":
            i = int(np.clip(np.searchsorted(L, x, side="left") - 1, 0, L.size - 2))
            u = 1.0 if abs(L[i + 1] - x) < 1e-9 else (x - L[i]) / (L[i + 1] - L[i])
        else:
            i = int(np.clip(np.searchsorted(L, x, side="right") - 1, 0, L.size - 2))
            u = (x - L[i]) / (L[i + 1] - L[i])
        h = L[i + 1] - L[i]
        A, B, Si = diag[i], diag[i + 1], S[i]
        C = (1 - u) ** 2 * A + u * (1 - u) * Si + u ** 2 * B
        dC = (-2 * (1 - u) * A + (1 - 2 * u) * Si + 2 * u * B) / h
        return self._finish(dC + 4.0 * self._theta(x) * C)


# ---------------------------------------------------------------------------
class R2Fixed(_Base):
    """R2Density with the narrow-cell edge handling fixed properly.

    Two changes vs the shipped R2Density:
      (a) the PCHIP is built PER CELL, so the interior node values are never
          de-duplicated away.  The shipped version concatenates all cells and
          then drops duplicate abscissae, which forces a single value at each
          table node and smears the true jump.
      (b) no linear extrapolation to the cell edges; the per-cell PCHIP is
          evaluated by its own (cubic) extrapolation over the <= h margin.
    """

    def __init__(self, *a, fd_h: float = 2.0, n_per_cell: int = 13, **kw):
        super().__init__(*a, **kw)
        self.fd_h = float(fd_h)
        self.n_per_cell = int(n_per_cell)

    def _density_splines(self, cos_gamma):
        key = ("R2F", round(float(cos_gamma), 12))
        c = self._cache.get(key)
        if c is not None:
            return c
        breaks = _cells(self.lam_lo, self.lam_hi)
        per_cell = []
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            h = min(self.fd_h, (hi - lo) / 8.0)
            n = max(4, min(self.n_per_cell, int((hi - lo) / max(h, 1e-9))))
            xs = np.linspace(lo + h, hi - h, n)
            fp = self._sample_f(cos_gamma, xs + h)
            fm = self._sample_f(cos_gamma, xs - h)
            d4 = np.asarray(self.bg.D(xs), float) ** 4
            dens = (fp - fm) / (2.0 * h) / d4[:, None, None]
            per_cell.append((lo, hi, [PchipInterpolator(xs, dens[:, a, b],
                                                        extrapolate=True)
                                      for a in range(N_COMP)
                                      for b in range(N_COMP)]))
        self._cache[key] = per_cell
        return per_cell

    def matrix(self, cos_gamma, lam):
        cells = self._density_splines(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        for lo, hi, sp in cells:
            if lo - 1e-9 <= x <= hi + 1e-9:
                break
        out = np.array([[sp[a * N_COMP + b](x) for b in range(N_COMP)]
                        for a in range(N_COMP)])
        return self._finish(out)


# ---------------------------------------------------------------------------
class R7Akima(_Base):
    """'A different scheme entirely': Akima on D^4 C over the 160 stock nodes,
    then differentiate.  Akima is the standard anti-ringing interpolant."""

    def _density_splines(self, cos_gamma):
        key = ("R7", round(float(cos_gamma), 12))
        c = self._cache.get(key)
        if c is not None:
            return c
        f = self._sample_f(cos_gamma, self._nodes)
        flat = [Akima1DInterpolator(self._nodes, f[:, a, b])
                for a in range(N_COMP) for b in range(N_COMP)]
        self._cache[key] = flat
        return flat

    def matrix(self, cos_gamma, lam):
        flat = self._density_splines(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        d4 = float(self.bg.D(x)) ** 4
        out = np.array([[flat[a * N_COMP + b].derivative()(x)
                         for b in range(N_COMP)] for a in range(N_COMP)]) / d4
        return self._finish(out)


class R8HermiteNodes(_Base):
    """'A C^1 Hermite through the table nodes': sample D^4 C at the 21 table
    nodes only, give each node a single PCHIP-style slope, differentiate.

    This is the scheme that assumes the density is CONTINUOUS at the nodes.
    """

    def _density_splines(self, cos_gamma):
        key = ("R8", round(float(cos_gamma), 12))
        c = self._cache.get(key)
        if c is not None:
            return c
        L = TABLE_LAM
        f = self._sample_f(cos_gamma, L)
        flat = []
        for a in range(N_COMP):
            for b in range(N_COMP):
                y = f[:, a, b]
                dydx = PchipInterpolator(L, y).derivative()(L)
                flat.append(CubicHermiteSpline(L, y, dydx))
        self._cache[key] = flat
        return flat

    def matrix(self, cos_gamma, lam):
        flat = self._density_splines(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        d4 = float(self.bg.D(x)) ** 4
        out = np.array([[flat[a * N_COMP + b](x, 1) for b in range(N_COMP)]
                        for a in range(N_COMP)]) / d4
        return self._finish(out)
