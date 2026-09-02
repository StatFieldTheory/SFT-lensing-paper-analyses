"""Bare equal-time driving-field cumulants for the Monte-Carlo Sachs integrator.

This module exposes, for the two-ray bundle ``f = (Phi_00, Re Psi_0, Im Psi_0)``
sampled along two sky directions, the *bare* white-noise driving cumulants that
the MC ODE Green's function then windows exactly once:

* ``sigma2_matrix(cos_gamma, lam)``  -- equal-time (3,3) angular covariance
  density (Part A; the anchor target).
* ``sigma2_6x6(cos_gamma, lam)``     -- the two-ray (6,6) block covariance.
* ``zeta_tensor(n1, n2, n3, lam)``   -- equal-time (3,3,3) third cumulant
  (Part B; wraps the project kappa3 vertex).
* ``zeta6(cos_gamma, lam)``          -- the two-ray (6,6,6) third cumulant.
* ``solve_Q(M2, Zeta6)``             -- the demo2 skewness deformation tensor
  (Part C).

Component convention everywhere: 0 = Phi_00 (kappa, spin-0), 1 = Re Psi_0
(gamma_+, spin-2), 2 = Im Psi_0 (gamma_x, spin-2).  Physical Mpc; lambda affine,
lambda = 0 at the observer.


THE 2-CUMULANT SOURCE AND THE ANCHOR  (read this before changing Sigma2)
=======================================================================

The MC's exact linear (Order-0) prediction is the line-of-sight double integral

    O0(gamma) = int_0^{lam_f} int_0^{lam_f} C_00(n1,t1; n2,t2) dt1 dt2          (1)

with C the Sachs correlation propagator ``<Phi(n1,t1) Phi(n2,t2)>``.  This is the
``integrate_over: "all"`` weak-lensing observable: kappa is the LOS-accumulated
Sachs scalar ``kappa = int_0^{lam_f} s(t) dt``, so the white-noise driver feeds
the ODE Green's function once and is then integrated along the ray.  Numerically
(1) reproduces analysis-3 ``C_corr_op_O0`` to 6 significant figures (verified by
running the sft-wick O0 workflow on a 3-gamma grid: 8.443189e-04 at 0.5').

For an equal-time (delta-collapsed-in-lambda) driver the per-leg response is the
project Green's function ``R(t, la) = [D(la)/D(t)]^2`` and (1) collapses to the
single integral

    O0(gamma) = int_0^{lam_f} Sigma2_00(cos; la) G(la)^2 dla,                   (2)
    G(la)     = int_la^{lam_f} [D(la)/D(t)]^2 dt,                               (3)

where Sigma2_ab(cos; la) is the BARE equal-time angular covariance density we
build here.  ``Sigma2`` is the object the MC samples; ``G`` is the LOS window the
ODE supplies for free.


WHY Sigma2 IS BUILT FROM corr_op, NOT FROM THE LIMBER C TABLE
------------------------------------------------------------
The DESIGN's first idea was to read Sigma2 from the diagonal of the bare limber
C table (``callables/C_propagator/limber/...``, ``cl_pp[:, j, j]`` etc.) and sum
it via the great-circle spin-2 Legendre/Wigner-d basis.  That FAILS the anchor,
and the failure is not a constant:

* The limber diagonal under-predicts O0 by a factor that *varies* with both ell
  (1.3e3 at ell=2 -> 5.7e2 at ell=1000) and lambda (Sigma_eff/Sigma_limber runs
  13 -> 100 across la), and even flips the sign of O0 at gamma ~ 90'.
* Root cause: the limber table stores the BARE delta-collapsed angular power
  density, while analysis-3's O0 (and the corr_op propagator) carry the full
  ASZ-2017 ``theta(lambda - chi)`` source-window convolution, which is
  ell-dependent and mixes lambda shells (the corr_op cl_pp is NOT diagonal in
  the two lambda axes, unlike the limber cl_pp).  No single ``[D/D]^n (1+z)^m``
  window reconciles the two.

The corr_op propagator IS that windowed object, and ``C_fn(n1, lam, n2, lam)``
returns the full (3,3) ``<Phi(lam) Phi(lam)>`` already including the window.  The
correct, anchor-passing bare equal-time density is therefore extracted from the
corr_op equal-time diagonal.  Since R = [D/D]^2 is scalar (iso_R), the same D^4
factor de-windows every component:

    Sigma2_ab(cos; lam) = (1 / D(lam)^4) d/dlam [ D(lam)^4 C_ab(cos; lam, lam) ]. (4)

Identity check: C_ab(cos; t, t) = int_0^t [D(la)/D(t)]^4 Sigma2_ab(cos; la) dla,
so D(t)^4 C(t,t) = int_0^t D(la)^4 Sigma2(la) dla and (4) is its exact inverse.
Sigma2_00 extracted this way is positive and smooth.


THE ANCHOR CONSTANT C0
----------------------
Feeding the corr_op-extracted Sigma2_00 through (2)-(3) reproduces analysis-3 O0
with a CLEAN, gamma-INDEPENDENT ratio over the physically-relevant small-angle
regime (0.5' .. ~150', where O0 > 0 and dominant):

    O0_eqtime(gamma) / O0_analysis3(gamma) = 0.881 +/- 0.014  (median 0.881).

The ~12% deficit is the off-diagonal (cross-shell) part of the corr_op kernel
that the strict equal-time / white-noise collapse drops; it is the same finite
lambda-correlation tail the DESIGN parameterises by ``sigma_lambda``.  Because it
is a single O(1) constant we absorb it into the driver as ``ANCHOR_C0 = 1/0.881``
so that the white-noise MC matches analysis-3 O0 in absolute normalisation.  (At
gamma > 200' both curves are ~1e-8 and oscillate about zero; the ratio there is
numerical noise on a vanishing signal and is not used to fix C0.)
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Final

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline

import background as _bg

# --- locate the project callables (corr_op propagator + kappa3 vertex) -------
_SACHS_SFT = Path(__file__).resolve().parents[2]  # .../sachs_sft
_CORR_OP_DIR = _SACHS_SFT / "callables" / "C_propagator" / "corr_op"
_KAPPA3_DIR = _SACHS_SFT / "callables" / "kappa3_vertex" / "equal_time_limber_cut15360_permaware"
for _d in (_CORR_OP_DIR, _KAPPA3_DIR):
    if str(_d) not in sys.path:
        sys.path.insert(0, str(_d))

import corr_op_C_callable as _corr_op  # noqa: E402  C_fn(n1, t1, n2, t2) -> (3,3)
# Permutation-aware vertex at ell_max = 15360 since 2026-09-02 (was the cut1000
# equal_time_limber_kappa3_callable); only a fresh Monte-Carlo run is affected.
import perm_aware_kappa3_callable as _k3  # noqa: E402  coupling_fn -> (3,3,3)

N_COMP: Final[int] = 3
N_RAY: Final[int] = 2
N6: Final[int] = N_COMP * N_RAY

# Anchor constant: corrects the strict-equal-time (white-noise) collapse to the
# analysis-3 O0 normalisation.  See module docstring "THE ANCHOR CONSTANT C0".
ANCHOR_RATIO: Final[float] = 0.881
ANCHOR_C0: Final[float] = 1.0 / ANCHOR_RATIO

# Reference rays for the two-ray reduced problem.
_N1: Final[NDArray[np.float64]] = np.array([0.0, 0.0, 1.0])


def _n2_of_cos(cos_gamma: float) -> NDArray[np.float64]:
    """Second ray at angular separation ``arccos(cos_gamma)`` in the x-z plane."""
    c = float(np.clip(cos_gamma, -1.0, 1.0))
    s = float(np.sqrt(max(0.0, 1.0 - c * c)))
    return np.array([s, 0.0, c])


# ---------------------------------------------------------------------------
# Part A -- equal-time (3,3) angular covariance density Sigma2_ab(cos; lam)
# ---------------------------------------------------------------------------
class Sigma2Builder:
    """Equal-time bare driving covariance density built from the corr_op propagator.

    For a fixed ``cos_gamma`` the (3,3) corr_op diagonal ``C_ab(t, t)`` is sampled
    on a dense lambda grid, multiplied by ``D(t)^4``, splined, and differentiated
    to give the bare equal-time density (eq. 4 in the module docstring).  The
    anchor constant ``ANCHOR_C0`` is applied so the MC matches analysis-3 O0.

    Parameters
    ----------
    background : background.Background, optional
        FLRW background providing ``D(lam)``.  A default is constructed if None.
    n_nodes : int
        Number of lambda nodes for the corr_op-diagonal sampling spline.
    lam_lo, lam_hi : float
        Sampling span (physical Mpc).  Defaults track the corr_op table support
        (406 .. 2328 Mpc); queries outside are clamped to the spline ends.
    apply_c0 : bool
        Multiply Sigma2 by ``ANCHOR_C0`` (default True).  Set False to inspect
        the raw equal-time density.
    """

    def __init__(self, background: _bg.Background | None = None, n_nodes: int = 160,
                 lam_lo: float = 406.0, lam_hi: float = 2328.0,
                 apply_c0: bool = True) -> None:
        self.bg = background if background is not None else _bg.Background()
        self.lam_lo = float(lam_lo)
        self.lam_hi = float(lam_hi)
        self.apply_c0 = bool(apply_c0)
        self._nodes = np.linspace(self.lam_lo, self.lam_hi, int(n_nodes))
        self._D4_nodes = np.asarray(self.bg.D(self._nodes), dtype=float) ** 4
        self._cache: dict[float, list[CubicSpline]] = {}

    # -- internal: per-cos spline of d/dlam[D^4 C_ab(t,t)] / D^4 ---------------
    def _density_splines(self, cos_gamma: float) -> list[CubicSpline]:
        key = round(float(cos_gamma), 12)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        n2 = _n2_of_cos(cos_gamma)
        # Sample the corr_op equal-time (3,3) diagonal on the lambda nodes.
        c_diag = np.empty((self._nodes.size, N_COMP, N_COMP), dtype=float)
        for i, t in enumerate(self._nodes):
            c_diag[i] = _corr_op.C_fn(_N1, float(t), n2, float(t))
        f = self._D4_nodes[:, None, None] * c_diag  # D^4 * C_ab(t,t)
        splines: list[CubicSpline] = []
        for a in range(N_COMP):
            row: list[CubicSpline] = []
            for b in range(N_COMP):
                row.append(CubicSpline(self._nodes, f[:, a, b]))
            splines.append(row)  # type: ignore[arg-type]
        # Flatten to a 9-list for cache compactness.
        flat = [splines[a][b] for a in range(N_COMP) for b in range(N_COMP)]
        self._cache[key] = flat
        return flat

    def matrix(self, cos_gamma: float, lam: float) -> NDArray[np.float64]:
        """Bare equal-time (3,3) covariance density at separation cos and shell lam."""
        flat = self._density_splines(cos_gamma)
        lam_c = float(np.clip(lam, self.lam_lo, self.lam_hi))
        d4 = float(self.bg.D(lam_c)) ** 4
        out = np.empty((N_COMP, N_COMP), dtype=float)
        for a in range(N_COMP):
            for b in range(N_COMP):
                out[a, b] = flat[a * N_COMP + b](lam_c, 1) / d4
        out = 0.5 * (out + out.T)  # enforce exact symmetry
        if self.apply_c0:
            out *= ANCHOR_C0
        return out


# Module-level default builder (lazy) so the free functions need no object.
_DEFAULT_BUILDER: Sigma2Builder | None = None


def _builder() -> Sigma2Builder:
    global _DEFAULT_BUILDER
    if _DEFAULT_BUILDER is None:
        _DEFAULT_BUILDER = Sigma2Builder()
    return _DEFAULT_BUILDER


def sigma2_matrix(cos_gamma: float, lam: float) -> NDArray[np.float64]:
    """Bare equal-time (3,3) angular covariance density (the Part-A object).

    ``cos_gamma = 1`` returns the within-ray block; ``cos_gamma = cos(gamma)``
    returns the cross-ray block.  See module docstring for the construction and
    the anchor.
    """
    return _builder().matrix(cos_gamma, lam)


def sigma2_blocks(cos_gamma: float, lam: float
                  ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Convenience: return ``(within_ray_block, cross_ray_block)`` at one lam.

    The within-ray block is ``sigma2_matrix(1, lam)``; the cross-ray block is
    ``sigma2_matrix(cos_gamma, lam)``.
    """
    return sigma2_matrix(1.0, lam), sigma2_matrix(cos_gamma, lam)


def sigma2_6x6(cos_gamma: float, lam: float) -> NDArray[np.float64]:
    """Two-ray (6,6) equal-time covariance density.

    Layout: indices 0..2 are ray 1's (Phi_00, Re Psi_0, Im Psi_0); indices 3..5
    are ray 2's.  Diagonal 3x3 blocks use the within-ray (cos=1) covariance; the
    off-diagonal blocks use the cross-ray (cos=gamma) covariance.
    """
    within, cross = sigma2_blocks(cos_gamma, lam)
    out = np.zeros((N6, N6), dtype=float)
    out[0:N_COMP, 0:N_COMP] = within
    out[N_COMP:N6, N_COMP:N6] = within
    out[0:N_COMP, N_COMP:N6] = cross
    out[N_COMP:N6, 0:N_COMP] = cross.T
    return out


# ---------------------------------------------------------------------------
# Part A anchor -- Order-0 linear prediction O0_MC(gamma)
# ---------------------------------------------------------------------------
def order0_window(background: _bg.Background, lam: NDArray[np.float64],
                  lam_f: float, n_gauss: int = 48) -> NDArray[np.float64]:
    """LOS window ``G(la) = int_la^{lam_f} [D(la)/D(t)]^2 dt`` (eq. 3), vectorised."""
    nodes, weights = leggauss(int(n_gauss))
    la = np.asarray(lam, dtype=float)
    D_la = np.asarray(background.D(la), dtype=float)
    out = np.empty_like(la)
    for i, l in enumerate(la):
        if l >= lam_f:
            out[i] = 0.0
            continue
        t = 0.5 * (lam_f - l) * nodes + 0.5 * (lam_f + l)
        jw = 0.5 * (lam_f - l) * weights
        out[i] = float(np.sum(jw * (D_la[i] / np.asarray(background.D(t))) ** 2))
    return out


def order0_mc(cos_gamma: float, lam_f: float | None = None,
              builder: Sigma2Builder | None = None,
              n_gauss: int = 48) -> float:
    """Analytic Order-0 prediction (eq. 2): ``int Sigma2_00 G^2 dla``.

    This is the MC's exact linear limit and the anchor target.  ``lam_f`` defaults
    to the analysis-3 baseline source (2313.029 Mpc).
    """
    b = builder if builder is not None else _builder()
    lf = float(lam_f) if lam_f is not None else _bg.LAM_SOURCE_BASELINE
    nodes, weights = leggauss(int(n_gauss))
    la = 0.5 * (lf - b.lam_lo) * nodes + 0.5 * (lf + b.lam_lo)
    jw = 0.5 * (lf - b.lam_lo) * weights
    sig00 = np.array([b.matrix(cos_gamma, float(l))[0, 0] for l in la])
    G = order0_window(b.bg, la, lf, n_gauss=n_gauss)
    return float(np.sum(jw * sig00 * G ** 2))


# ---------------------------------------------------------------------------
# Part B -- equal-time (3,3,3) third cumulant zeta_abc(cos-triple; lambda)
# ---------------------------------------------------------------------------
def zeta_tensor(n1: NDArray[np.float64], n2: NDArray[np.float64],
                n3: NDArray[np.float64], lam: float) -> NDArray[np.float64]:
    """Equal-time bare (3,3,3) third cumulant at the three directions and shell lam.

    Thin wrapper over the project's ``equal_time_limber`` kappa3 vertex
    (``ALREADY_R_CONTRACTED=False``, ``EQUAL_TIME=True``, lambda-density, per-leg
    (1+z)^4 Sachs prefactor already baked in).  Returns the LOCAL real-basis
    cumulant ``K_abc`` in component order (Phi_00, Re Psi_0, Im Psi_0).
    """
    return np.asarray(
        _k3.coupling_fn([n1, n2, n3], [float(lam)] * 3), dtype=float,
    )


def _ray_of(idx6: int) -> int:
    """Which ray (0 or 1) does a 6-index leg belong to."""
    return idx6 // N_COMP


def _comp_of(idx6: int) -> int:
    """Component (0..2) within a ray for a 6-index leg."""
    return idx6 % N_COMP


def zeta6(cos_gamma: float, lam: float) -> NDArray[np.float64]:
    """Two-ray (6,6,6) third cumulant density assembled from the 4 ray triples.

    Each 6-index leg (0..5) carries a ray (0 or 1) and a component (0..2).  For a
    given (a6, b6, c6) the three rays select one of the four direction triples
    (n1n1n1, n1n1n2, n1n2n2, n2n2n2); the (3,3,3) zeta for that triple supplies
    the value at the legs' component indices.  Symmetric by construction (the
    kappa3 vertex is symmetric and the placement respects leg identity).
    """
    n1 = _N1
    n2 = _n2_of_cos(cos_gamma)
    rays = (n1, n2)
    # Precompute the 4 distinct (ray-pattern) zeta tensors keyed by sorted rays
    # is not enough -- leg order matters for the placement, but the vertex is
    # symmetric so we evaluate per (r_a, r_b, r_c) ordered triple (8 calls; the
    # vertex caches identical-cos triples internally).
    out = np.zeros((N6, N6, N6), dtype=float)
    zcache: dict[tuple[int, int, int], NDArray[np.float64]] = {}
    for a6 in range(N6):
        ra, ca = _ray_of(a6), _comp_of(a6)
        for b6 in range(N6):
            rb, cb = _ray_of(b6), _comp_of(b6)
            for c6 in range(N6):
                rc, cc = _ray_of(c6), _comp_of(c6)
                key = (ra, rb, rc)
                z = zcache.get(key)
                if z is None:
                    z = zeta_tensor(rays[ra], rays[rb], rays[rc], lam)
                    zcache[key] = z
                out[a6, b6, c6] = z[ca, cb, cc]
    return out


# ---------------------------------------------------------------------------
# Part C -- skewness deformation tensor Q
# ---------------------------------------------------------------------------
def solve_Q(M2: NDArray[np.float64], Zeta6: NDArray[np.float64],
            rcond: float = 1e-8) -> NDArray[np.float64]:
    """Solve the demo2 skewness deformation relation for Q.

    Given the per-step Gaussian covariance ``M2`` (symmetric, PD up to a tiny
    ridge) and the target third cumulant ``Zeta6``, find ``Q`` (symmetric in its
    last two indices) such that

        Zeta6_abc = Q_amn M2_mb M2_nc + Q_bmn M2_ma M2_nc + Q_cmn M2_ma M2_nb.   (5)

    Method (M2-eigenbasis diagonalisation).  In the basis that diagonalises M2
    (``M2 = U diag(d) U^T``), define ``Q'_pqr = U_ap U_bq U_cr Q_abc`` and
    ``Z'_pqr = U_ap U_bq U_cr Zeta6_abc``.  Equation (5) becomes diagonal:

        Z'_pqr = Q'_pqr d_q d_r + Q'_qpr d_p d_r + Q'_rqp d_p d_q,

    and with ``Q'`` symmetric in its last two indices (q,r) the three terms have
    coefficients ``(d_q d_r, d_p d_r, d_p d_q)``.  Symmetrising Z' over all three
    indices (the cumulant is fully symmetric) and solving componentwise gives

        Q'_pqr = Zsym'_pqr / (d_p d_q + d_q d_r + d_r d_p),

    then ``Q = U Q' U^T`` (rotate back on every index).  Reduces to the scalar
    demo2 relation ``z = 3 Q m^2 => Q = z/(3 m^2)`` (one eigenvalue d=m).
    """
    M2 = np.asarray(M2, dtype=float)
    Zeta6 = np.asarray(Zeta6, dtype=float)
    n = M2.shape[0]
    if M2.shape != (n, n):
        raise ValueError(f"M2 must be square; got {M2.shape}")
    if Zeta6.shape != (n, n, n):
        raise ValueError(f"Zeta6 must be ({n},{n},{n}); got {Zeta6.shape}")

    # Tiny ridge for numerical PD-ness, scaled to the matrix magnitude.
    scale = float(np.max(np.abs(np.diag(M2)))) or 1.0
    M2r = M2 + np.eye(n) * (1e-12 * scale)
    d, U = np.linalg.eigh(M2r)  # M2r = U diag(d) U^T, U orthonormal

    # Rotate Zeta6 into the eigenbasis: Z'_pqr = U_ap U_bq U_cr Zeta6_abc.
    Zp = np.einsum("ap,bq,cr,abc->pqr", U, U, U, Zeta6, optimize=True)
    # Fully symmetrise (the cumulant is symmetric; guards against input asymmetry).
    Zp = (Zp + Zp.transpose(0, 2, 1) + Zp.transpose(1, 0, 2)
          + Zp.transpose(1, 2, 0) + Zp.transpose(2, 0, 1)
          + Zp.transpose(2, 1, 0)) / 6.0
    denom = (d[:, None, None] * d[None, :, None]
             + d[None, :, None] * d[None, None, :]
             + d[None, None, :] * d[:, None, None])
    Qp = Zp / denom  # denom > 0 since d > 0 after the ridge
    # Regularise rank-deficient M2 (e.g. two near-aligned rays -> antisymmetric
    # combos have ~0 variance and ~0 cumulant): zero Q' in any eigendirection
    # below a relative floor so 1/(small d)^2 does not blow up the deformation.
    good = d > (rcond * float(d.max()))
    keep = good[:, None, None] & good[None, :, None] & good[None, None, :]
    Qp = np.where(keep, Qp, 0.0)
    # Rotate back: Q_abc = U_ap U_bq U_cr Q'_pqr.
    Q = np.einsum("ap,bq,cr,pqr->abc", U, U, U, Qp, optimize=True)
    return 0.5 * (Q + Q.transpose(0, 2, 1))  # symmetric in last two indices


def cum3_from_Q(Q: NDArray[np.float64], M2: NDArray[np.float64]) -> NDArray[np.float64]:
    """Forward map (eq. 5): reconstruct the third cumulant produced by ``Q``.

    Used to verify ``solve_Q`` round-trips.
    """
    t1 = np.einsum("amn,mb,nc->abc", Q, M2, M2, optimize=True)
    t2 = np.einsum("bmn,ma,nc->abc", Q, M2, M2, optimize=True)
    t3 = np.einsum("cmn,ma,nb->abc", Q, M2, M2, optimize=True)
    return t1 + t2 + t3


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
def _self_test() -> None:
    import math

    print("=" * 72)
    print("driver_stats self-test")
    print("=" * 72)

    # ---- Part A: O0 anchor against analysis-3 ----------------------------
    anc_path = (_SACHS_SFT / "sftwick_outputs" / "2PCF" / "C_corr_op_O0"
                / "xi_C_corr_op_O0.npz")
    # allow_pickle: trusted local analysis-3 artifact; x/y are object arrays of
    # direction tuples written by this project's own sft-wick sweep.
    ad = np.load(anc_path, allow_pickle=True)
    m = (ad["a"] == 0) & (ad["b"] == 0) & (ad["order"] == 0)
    xs, ys, vs = ad["x"][m], ad["y"][m], ad["value"][m]
    gammas, coss = [], []
    for xi, yi in zip(xs, ys):
        xi = np.asarray(xi, float); yi = np.asarray(yi, float)
        c = float(np.clip(np.dot(xi / np.linalg.norm(xi), yi / np.linalg.norm(yi)),
                          -1.0, 1.0))
        gammas.append(math.degrees(math.acos(c)) * 60.0)
        coss.append(c)
    order = np.argsort(gammas)
    gammas = np.asarray(gammas)[order]
    coss = np.asarray(coss)[order]
    vs = np.asarray(vs, float)[order]

    builder = _builder()
    o0 = np.array([order0_mc(c, lam_f=float(ad["t_final"][0]), builder=builder)
                   for c in coss])
    ratio = o0 / vs

    print("\n[Part A] Order-0 anchor  O0_MC vs analysis-3 C_corr_op_O0 (kk):")
    print(f"{'gamma[arcmin]':>14} {'O0_MC':>13} {'analysis3':>13} {'ratio':>9}")
    for label_g in (0.5, 1.0, 10.0, 100.0, 1000.0):
        i = int(np.argmin(np.abs(gammas - label_g)))
        print(f"{gammas[i]:14.3f} {o0[i]:13.5e} {vs[i]:13.5e} {ratio[i]:9.4f}")
    pos = (vs > 0) & (gammas < 200.0)
    print(f"\n  small-gamma ratio (gamma<200', O0>0): "
          f"median={np.median(ratio[pos]):.4f} "
          f"std={np.std(ratio[pos]):.4f} "
          f"[{ratio[pos].min():.4f}, {ratio[pos].max():.4f}]")
    ok_anchor = abs(np.median(ratio[pos]) - 1.0) < 0.10 and np.std(ratio[pos]) < 0.05
    print(f"  ANCHOR {'PASS' if ok_anchor else 'FAIL'} "
          f"(median within 10% of 1, spread < 5%)")

    # ---- Part B: zeta magnitudes -----------------------------------------
    print("\n[Part B] equal-time zeta_abc magnitudes (two-ray triples):")
    n1 = _N1
    lam_mid = 1361.6
    for g_arcmin in (1.0, 100.0):
        g = g_arcmin / 60.0 * math.pi / 180.0
        n2 = _n2_of_cos(math.cos(g))
        print(f"  gamma = {g_arcmin:6.1f}'  (lam={lam_mid:.1f}):")
        for trip, label in (((n1, n1, n1), "n1n1n1"), ((n1, n1, n2), "n1n1n2"),
                            ((n1, n2, n2), "n1n2n2"), ((n2, n2, n2), "n2n2n2")):
            z = zeta_tensor(*trip, lam_mid)
            print(f"    {label}: TTT={z[0,0,0]:.3e}  max|zeta|={np.max(np.abs(z)):.3e}")

    # ---- Part C: solve_Q round-trip + scalar reduction -------------------
    print("\n[Part C] solve_Q round-trip:")
    rng = np.random.default_rng(0)
    # Build a realistic (6,6) M2 and (6,6,6) Zeta6 from the drivers.
    cosg = math.cos(10.0 / 60.0 * math.pi / 180.0)
    dlam = 5.0
    M2 = sigma2_6x6(cosg, lam_mid) * dlam
    Zeta6 = zeta6(cosg, lam_mid) * dlam
    Q = solve_Q(M2, Zeta6)
    recon = cum3_from_Q(Q, M2)
    err = np.max(np.abs(recon - Zeta6)) / max(np.max(np.abs(Zeta6)), 1e-300)
    print(f"  physical (6,6,6): max rel round-trip err = {err:.2e}")

    # Random well-conditioned case (the driver zeta is near-degenerate).
    A = rng.standard_normal((N6, N6))
    M2b = A @ A.T + np.eye(N6)
    Zb = rng.standard_normal((N6, N6, N6))
    Zb = (Zb + Zb.transpose(0, 2, 1) + Zb.transpose(1, 0, 2)
          + Zb.transpose(1, 2, 0) + Zb.transpose(2, 0, 1)
          + Zb.transpose(2, 1, 0)) / 6.0
    Qb = solve_Q(M2b, Zb)
    errb = np.max(np.abs(cum3_from_Q(Qb, M2b) - Zb)) / np.max(np.abs(Zb))
    print(f"  random  (6,6,6): max rel round-trip err = {errb:.2e}")

    # Scalar reduction: z = 3 Q m^2  =>  Q = z / (3 m^2).
    m_s, z_s = 0.7, 0.013
    Qs = solve_Q(np.array([[m_s]]), np.array([[[z_s]]]))
    Q_expected = z_s / (3.0 * m_s ** 2)
    print(f"  scalar: solve_Q -> {Qs[0,0,0]:.6e}, "
          f"z/(3 m^2) -> {Q_expected:.6e}, "
          f"match={'OK' if abs(Qs[0,0,0]-Q_expected) < 1e-12 else 'FAIL'}")

    ok_q = errb < 1e-8 and abs(Qs[0, 0, 0] - Q_expected) < 1e-12
    print(f"  solve_Q {'PASS' if ok_q else 'FAIL'}")

    # ---- Optional anchor plot --------------------------------------------
    _try_plot(gammas, o0, vs, ratio)
    print("\n" + ("ALL PASS" if (ok_anchor and ok_q) else "CHECK FAILURES ABOVE"))


def _try_plot(gammas: NDArray[np.float64], o0: NDArray[np.float64],
              vs: NDArray[np.float64], ratio: NDArray[np.float64]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - plotting is optional
        print(f"  [plot skipped: {exc}]")
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.loglog(gammas, np.abs(o0), "o-", ms=3, label="O0_MC (driver_stats)")
    ax1.loglog(gammas, np.abs(vs), "s--", ms=3, label="analysis-3 C_corr_op_O0")
    ax1.set_xlabel("gamma [arcmin]")
    ax1.set_ylabel("|O0_kk|")
    ax1.set_title("Order-0 anchor (kk)")
    ax1.legend()
    ax1.grid(True, which="both", alpha=0.3)
    ax2.semilogx(gammas, ratio, "o-", ms=3)
    ax2.axhline(1.0, color="k", lw=0.8)
    ax2.set_xlabel("gamma [arcmin]")
    ax2.set_ylabel("O0_MC / analysis-3")
    ax2.set_ylim(0.0, 2.0)
    ax2.set_title("anchor ratio")
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    out_dir = Path(__file__).resolve().parent / "figures"
    out_dir.mkdir(exist_ok=True)
    out = out_dir / "anchor_O0_check.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"  [anchor plot -> {out}]")


__all__ = [
    "Sigma2Builder",
    "sigma2_matrix",
    "sigma2_blocks",
    "sigma2_6x6",
    "order0_mc",
    "order0_window",
    "zeta_tensor",
    "zeta6",
    "solve_Q",
    "cum3_from_Q",
    "ANCHOR_C0",
]


if __name__ == "__main__":
    _self_test()
