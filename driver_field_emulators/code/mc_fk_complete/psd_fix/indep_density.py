"""Independent re-derivation of the Sigma2 equal-time driving density.

Clean-room: written without reading ``psd_fix/sigma2_repaired.py``.

WHAT THE OBJECT IS
------------------
``driver_stats`` defines the bare equal-time driving covariance density by
inverting the equal-time value of the corr_op propagator (module docstring,
eq. 4):

    Sigma2_ab(cos; lam) = (1 / D(lam)^4) d/dlam [ D(lam)^4 C_ab(cos; lam, lam) ]

with ``C_ab(cos; t1, t2) = corr_op.C_fn(n1, t1, n2, t2)``.

WHAT C ACTUALLY IS (measured, not assumed -- see indep_probe2/3.py)
-------------------------------------------------------------------
``corr_op`` interpolates a 21 x 21 lambda table with
``scipy.interpolate.RegularGridInterpolator(..., method="linear")``
(canoes ``table.py:958-965``; ``make_C_fn_lambda_project`` takes the default
``interp_method="linear"``).  A **bilinear** interpolant restricted to the
diagonal ``t1 = t2 = t`` is, inside the cell ``[x_i, x_{i+1}]`` with
``u = (t - x_i)/h_i``,

    C(t,t) = C_ii (1-u)^2 + 2 C_{i,i+1} u(1-u) + C_{i+1,i+1} u^2 ,

i.e. **exactly a quadratic in t inside every cell**, continuous in value at the
nodes but with a JUMPING FIRST DERIVATIVE:

    dC/dt|_{x_i^-} = 2 (C_ii - C_{i-1,i}) / h_{i-1}
    dC/dt|_{x_i^+} = 2 (C_{i,i+1} - C_ii) / h_i .

Both statements are verified numerically in ``__main__`` below (deg-2 fit
residual ~1e-15 relative in all 20 cells x 9 entries; value jump ~1e-16;
derivative jump of order 100%).

THE EVALUATOR
-------------
Because C is exactly quadratic inside a cell, the derivative needs no finite
difference at all.  Writing the product rule for the D^4 factor,

    Sigma2 = (1/D^4) d/dlam [D^4 C] = dC/dlam + 4 (D'/D) C ,           (*)

the ONLY numerical ingredients are (a) the per-cell quadratic, obtained by an
exact local fit that never crosses a node, and (b) the background log-derivative
D'/D.  Route (*) is the primary evaluator (``matrix``).

An entirely separate route -- a 5-point central finite difference of
``g(t) = D(t)^4 C(t,t)`` with the whole stencil confined to one cell -- is
implemented in ``matrix_fd`` and used to *demonstrate* convergence rather than
assume it.  The two routes agree to ~1e-11 relative (see ``__main__``).

Because the density is genuinely DISCONTINUOUS at the 21 table nodes, a query
exactly at a node has two answers.  ``side="left"/"right"`` selects one;
``side="auto"`` (default) takes the cell containing lam, resolving an exact node
hit to the right cell -- matching what a cell-local scheme does.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Final

import numpy as np
from numpy.typing import NDArray

_MC_FK = Path(__file__).resolve().parent.parent
if str(_MC_FK) not in sys.path:
    sys.path.insert(0, str(_MC_FK))
import _bootstrap  # noqa: E402

_ds, _bg, _core = _bootstrap.wire()
import corr_op_C_callable as _corr_op  # noqa: E402

N_COMP: Final[int] = 3
ANCHOR_C0: Final[float] = _ds.ANCHOR_C0
_N1: Final[NDArray[np.float64]] = np.array([0.0, 0.0, 1.0])

# The corr_op table's own lambda nodes -- the only places where the diagonal
# of the bilinear interpolant is non-smooth.
NODES: Final[NDArray[np.float64]] = np.asarray(
    _corr_op._TABLE.cl_table.lambda_grid, dtype=float)


def _n2_of_cos(cos_gamma: float) -> NDArray[np.float64]:
    c = float(np.clip(cos_gamma, -1.0, 1.0))
    s = float(np.sqrt(max(0.0, 1.0 - c * c)))
    return np.array([s, 0.0, c])


# WARNING (measured, see indep_check_batch3.py): corr_op's batch and scalar
# paths are NOT interchangeable.  C_fn_batch groups samples by a ROUNDED
# (cos gamma, psi1, psi2) key and caches the resulting angular-channel table;
# whichever path touches a given cos FIRST determines what BOTH return
# afterwards.  Batch-first vs scalar-only differ by up to ~0.17% at gamma = 1'.
# driver_stats.Sigma2Builder uses the SCALAR C_fn, so the default here is
# scalar too -- anything else would compare two different angles.
USE_BATCH_DEFAULT = False


def _diag(t: NDArray[np.float64], cos_gamma: float,
          use_batch: bool = USE_BATCH_DEFAULT) -> NDArray[np.float64]:
    """C_ab(cos; t, t) for an array of t -> (n, 3, 3)."""
    ta = np.asarray(t, dtype=float)
    n2 = _n2_of_cos(cos_gamma)
    if use_batch:
        return np.asarray(_corr_op.C_fn_batch(_N1, ta, n2, ta), dtype=float)
    return np.asarray([_corr_op.C_fn(_N1, float(x), n2, float(x)) for x in ta],
                      dtype=float)


_diag_batch = _diag


class IndepDensity:
    """Drop-in replacement for ``driver_stats.Sigma2Builder``.

    Same public surface (``bg``, ``lam_lo``, ``lam_hi``, ``apply_c0``,
    ``matrix(cos_gamma, lam)``) so it can be handed to ``driver_stats.order0_mc``
    and ``sachs_mc_core._assemble_6x6`` unchanged.  ``n_nodes`` is accepted and
    ignored (there is no global sampling grid here).
    """

    def __init__(self, background=None, n_nodes: int = 0,
                 lam_lo: float = 406.0, lam_hi: float = 2328.0,
                 apply_c0: bool = True, n_fit: int = 9, deg: int = 2,
                 inset: float = 0.02,
                 use_batch: bool = USE_BATCH_DEFAULT) -> None:
        self.bg = background if background is not None else _bg.Background()
        self.lam_lo = float(lam_lo)
        self.lam_hi = float(lam_hi)
        self.apply_c0 = bool(apply_c0)
        self.n_fit = int(n_fit)
        self.deg = int(deg)
        self.inset = float(inset)
        self.use_batch = bool(use_batch)
        self._dD = self.bg._D_spline.derivative(1)   # D'(lam) from the same spline
        self._cache: dict[float, tuple] = {}
        self._resid: dict[float, float] = {}

    # -- per-cell exact local polynomial of C_ab(t,t) -------------------------
    def _fit(self, cos_gamma: float):
        key = round(float(cos_gamma), 12)
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        ncell = NODES.size - 1
        # coefficients in the LOCAL variable u = (t - mid) / half
        coef = np.empty((ncell, N_COMP, N_COMP, self.deg + 1), dtype=float)
        worst = 0.0
        ts = []
        for i in range(ncell):
            lo, hi = NODES[i], NODES[i + 1]
            eps = self.inset * (hi - lo)
            ts.append(np.linspace(lo + eps, hi - eps, self.n_fit))
        Call = _diag(np.concatenate(ts), cos_gamma, self.use_batch)
        off = 0
        for i in range(ncell):
            lo, hi = NODES[i], NODES[i + 1]
            mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
            t = ts[i]
            u = (t - mid) / half
            blk = Call[off:off + self.n_fit]
            off += self.n_fit
            for a in range(N_COMP):
                for b in range(N_COMP):
                    y = blk[:, a, b]
                    c = np.polyfit(u, y, self.deg)
                    coef[i, a, b] = c
                    sc = np.abs(y).max()
                    if sc > 0.0:
                        worst = max(worst, np.abs(y - np.polyval(c, u)).max() / sc)
        out = (coef, worst)
        self._cache[key] = out
        self._resid[key] = worst
        return out

    def _cell(self, lam: float, side: str) -> int:
        lam = float(np.clip(lam, self.lam_lo, self.lam_hi))
        i = int(np.searchsorted(NODES, lam, side="right")) - 1
        i = min(max(i, 0), NODES.size - 2)
        if side == "left" and np.isclose(lam, NODES[i], rtol=0, atol=1e-9) and i > 0:
            i -= 1
        if side == "right" and np.isclose(lam, NODES[i + 1], rtol=0,
                                          atol=1e-9) and i < NODES.size - 2:
            i += 1
        return i

    def C_and_dC(self, cos_gamma: float, lam: float, side: str = "auto"):
        """Return ``(C_ab, dC_ab/dlam)`` at ``lam`` from the cell-local polynomial."""
        coef, _ = self._fit(cos_gamma)
        i = self._cell(lam, side)
        lo, hi = NODES[i], NODES[i + 1]
        mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
        u = (float(np.clip(lam, self.lam_lo, self.lam_hi)) - mid) / half
        C = np.empty((N_COMP, N_COMP)); dC = np.empty((N_COMP, N_COMP))
        for a in range(N_COMP):
            for b in range(N_COMP):
                c = coef[i, a, b]
                C[a, b] = np.polyval(c, u)
                dC[a, b] = np.polyval(np.polyder(c), u) / half
        return C, dC

    def matrix(self, cos_gamma: float, lam: float,
               side: str = "auto") -> NDArray[np.float64]:
        """Bare equal-time (3,3) density  dC/dlam + 4 (D'/D) C  (eq. *)."""
        lam_c = float(np.clip(lam, self.lam_lo, self.lam_hi))
        C, dC = self.C_and_dC(cos_gamma, lam_c, side)
        Dv = float(self.bg.D(lam_c))
        dlnD = float(self._dD(lam_c)) / Dv
        out = dC + 4.0 * dlnD * C
        out = 0.5 * (out + out.T)
        if self.apply_c0:
            out = out * ANCHOR_C0
        return out

    # -- independent cross-check route: FD of g = D^4 C, confined to one cell --
    def matrix_fd(self, cos_gamma: float, lam: float, h: float | None = None,
                  side: str = "auto") -> NDArray[np.float64]:
        """5-point central FD of ``d/dlam[D^4 C]`` / D^4, stencil inside one cell."""
        lam_c = float(np.clip(lam, self.lam_lo, self.lam_hi))
        i = self._cell(lam_c, side)
        lo, hi = NODES[i], NODES[i + 1]
        room = min(lam_c - lo, hi - lam_c)
        if room <= 0.0:                       # exactly on a node: shift inward
            lam_c = lam_c + 1e-6 * (hi - lo) if side != "left" else lam_c - 1e-6 * (hi - lo)
            room = min(lam_c - lo, hi - lam_c)
        hh = float(h) if h is not None else min(1.0, 0.4 * room)
        hh = min(hh, 0.49 * room)
        t = lam_c + hh * np.array([-2.0, -1.0, 1.0, 2.0])
        C = _diag(t, cos_gamma, self.use_batch)
        D4 = np.asarray(self.bg.D(t), dtype=float) ** 4
        g = D4[:, None, None] * C
        gp = (g[0] - 8.0 * g[1] + 8.0 * g[2] - g[3]) / (12.0 * hh)
        out = gp / float(self.bg.D(lam_c)) ** 4
        out = 0.5 * (out + out.T)
        if self.apply_c0:
            out = out * ANCHOR_C0
        return out


def assemble6(builder, cos_gamma: float, lam: float, side: str = "auto"):
    """(6,6) two-ray covariance, same layout as sachs_mc_core._assemble_6x6."""
    try:
        w = builder.matrix(1.0, lam, side=side)
        x = builder.matrix(float(cos_gamma), lam, side=side)
    except TypeError:
        w = builder.matrix(1.0, lam)
        x = builder.matrix(float(cos_gamma), lam)
    out = np.zeros((6, 6))
    out[0:3, 0:3] = w; out[3:6, 3:6] = w
    out[0:3, 3:6] = x; out[3:6, 0:3] = x.T
    return out


if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=False)
    bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
    ev = IndepDensity(background=bgd, apply_c0=False)

    print("corr_op table nodes (%d):" % NODES.size)
    print("  ", np.array2string(NODES, precision=0, max_line_width=200))

    print("\n=== A. cell-local polynomial is EXACT (no truncation to converge) ===")
    for gam in (0.5, 1.0, 5.0, 17.3, 60.0):
        cg = float(np.cos(np.deg2rad(gam / 60.0)))
        r2 = IndepDensity(background=bgd, apply_c0=False, deg=2)._fit(cg)[1]
        r3 = IndepDensity(background=bgd, apply_c0=False, deg=3)._fit(cg)[1]
        r4 = IndepDensity(background=bgd, apply_c0=False, deg=4)._fit(cg)[1]
        print(f"  gamma={gam:6.2f}'  worst rel fit residual over 20 cells x 9 entries: "
              f"deg2 {r2:.2e}  deg3 {r3:.2e}  deg4 {r4:.2e}")

    print("\n=== B. degree convergence of the DERIVATIVE (deg 2 vs 3 vs 4) ===")
    cg = float(np.cos(np.deg2rad(17.3 / 60.0)))
    e2 = IndepDensity(background=bgd, apply_c0=False, deg=2)
    e3 = IndepDensity(background=bgd, apply_c0=False, deg=3, n_fit=11)
    e4 = IndepDensity(background=bgd, apply_c0=False, deg=4, n_fit=13)
    lams = np.array([420.0, 700.5, 1111.0, 1299.0, 1301.0, 1800.0, 2295.0, 2305.0])
    print(f"  {'lam':>8} {'deg2 [0,0]':>14} {'|d3-d2|/d2':>12} {'|d4-d2|/d2':>12}")
    for l in lams:
        v2 = e2.matrix(cg, l)[0, 0]; v3 = e3.matrix(cg, l)[0, 0]; v4 = e4.matrix(cg, l)[0, 0]
        print(f"  {l:>8.1f} {v2:>14.7e} {abs(v3-v2)/abs(v2):>12.2e} "
              f"{abs(v4-v2)/abs(v2):>12.2e}")

    print("\n=== C. FD step convergence of the independent route (lam=1111, 17.3') ===")
    lam = 1111.0
    ref = ev.matrix(cg, lam)[0, 0]
    print(f"  analytic route  dC/dlam + 4(D'/D)C = {ref:.10e}")
    print(f"  {'h [Mpc]':>10} {'FD 5-point':>18} {'rel diff':>12}")
    for h in (30.0, 10.0, 3.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 1e-3, 1e-4):
        v = ev.matrix_fd(cg, lam, h=h)[0, 0]
        print(f"  {h:>10.4g} {v:>18.10e} {abs(v-ref)/abs(ref):>12.2e}")

    print("\n=== D. analytic vs FD over many lam and both blocks ===")
    rel = []
    for cgx in (1.0, cg):
        for l in np.linspace(410.0, 2310.0, 120):
            a = ev.matrix(cgx, l); b = ev.matrix_fd(cgx, l)
            sc = np.abs(a).max()
            rel.append(np.abs(a - b).max() / sc)
    rel = np.array(rel)
    print(f"  240 (block, lam) pairs: max rel diff {rel.max():.3e}, "
          f"median {np.median(rel):.3e}")

    print("\n=== E. the density is DISCONTINUOUS at the table nodes ===")
    print(f"  {'node':>8} {'Sigma00 left':>16} {'Sigma00 right':>16} {'jump/left':>12}")
    for x in NODES[1:-1]:
        L = ev.matrix(cg, float(x), side="left")[0, 0]
        R = ev.matrix(cg, float(x), side="right")[0, 0]
        print(f"  {x:>8.0f} {L:>16.7e} {R:>16.7e} {(R-L)/abs(L):>12.4f}")
