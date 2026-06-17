"""Independent analytic equal-time FK anchor for the MC Sachs 2-point check.

WHAT THIS IS
============
A from-scratch quadrature of the Order-2 FK channel of ``<kappa_a kappa_b>(gamma)``
that reproduces analysis-3's ``C_corr_op_K_limber_FK`` WITHOUT going through the
sft-wick diagram engine.  It uses the SAME bare equal-time zeta callable
(``equal_time_limber``), the SAME scalar response ``R(t, la) = [D(la)/D(t)]^2``,
and the SAME F-tensor.  So it cross-checks the sft-wick FK *diagram assembly*
(the equal-time time-collapse, the causal measure, the F-K index contraction and
the cosine-triple geometry) -- it does NOT independently re-derive the zeta
table's internal units (both this quadrature and analysis-3 read the same table).


THE FK DIAGRAM (the expression my first handoff got wrong)
==========================================================
In MSR the driver 3-cumulant enters as a 3-response (psi-psi-psi) non-local
vertex K; the F-vertex is local (1 psi, 2 phi).  The Order-2 FK term of
``<phi_a(n1) phi_b(n2)>`` has FOUR response (R) propagators and ZERO C
propagators (4 psi legs -> 4 phi legs; psi-psi = 0 in MSR).  The connected
topology (verbatim structure of demo2's validated ``_fk_spatial_integral``):

    obs a  --R-->  F-vertex  ==R,R==>  two K-legs ;   third K-leg  --R-->  obs b

so for a full-spacetime K it is the 4-D time integral

    xi_FK_ab = sum_{diagrams x<->y} sum_{m,n}
        F_{a m n}  int dtau0 dtau1 dtau2 dtau3
            R(t_f,tau0) R(t_f,tau1) R(tau0,tau2) R(tau0,tau3)
            K_{b m n}(tau1, tau2, tau3 ; directions) .

Two of the three K legs (tau2, tau3) sit on the F-vertex (direction n_a); the
third (tau1) sits on the OTHER observable (direction n_b).  Hence the external
index a rides the F-vertex (F_{a m n}) and only b rides a K leg (K_{b m n}); the
F's two phi indices m,n contract the two co-located K legs.  The cosine triple of
the three K legs is therefore (n_a.n_a, n_a.n_b, n_a.n_b) = (1, cos gamma,
cos gamma) -- two legs co-located, one separated -- NOT a generic triple.

EQUAL-TIME COLLAPSE.  Our zeta is equal-time (tau1 = tau2 = tau3 = lam_v), so the
4-D integral collapses to 2-D.  For the convergence observable
(``integrate_over='all'``, kappa = int_0^{lam_f} s dt) each external R becomes the
LOS window W:

    W(la) = int_la^{lam_f} [D(la)/D(t)]^2 dt        (driver_stats.order0_window)

and the two F<->K responses give [R(tau,lam_v)]^2 = [D(lam_v)/D(tau)]^4.  Result:

    xi_FK_ab(gamma) = int_0^{lam_f} dtau int_0^{tau} dlam_v
            W(tau) W(lam_v) [D(lam_v)/D(tau)]^4  G_ab(lam_v; gamma) ,      (FK)

    G_ab(lam_v; gamma) = sum_{m,n} [ F_{a m n} Z_{b m n} + F_{b m n} Z_{a m n} ]
                       = (M + M^T)_{ab},   M_{ab} = sum_{m,n} F_{a m n} Z_{b m n},

with Z = Z(lam_v; cos gamma) the bare equal-time (3,3,3) cumulant at the triple
(1, cos gamma, cos gamma).  The (M + M^T) is the x<->y symmetrization (the two
diagrams share the same cosine triple, so they share Z).  The overall MSR/
combinatorial prefactor is 1 in this convention (calibrated on demo2's validated
FK, where K is diagonal and (FK) reduces to (F[a,b,b]+F[b,a,a]) * integral).

For kk (a=b=0):  G_00 = 2 * sum_{m,n} F_{0 m n} Z_{0 m n}
                      = -2 ( Z_000 + Z_011 + Z_022 )  (since F_{0mn}=-delta_{mn})
                      = -2 ( zeta_TTT + zeta_Bmod ) .
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Final

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.typing import NDArray
from scipy.interpolate import CubicSpline

import background as _bg
import driver_stats as _ds

# Paper Table-1 F-vertex (= scripts/inputs/F_tensor.npy; verified equal upstream).
_F_PATH: Final[Path] = (
    _ds._SACHS_SFT / "scripts" / "inputs" / "F_tensor.npy"
)
_F: Final[NDArray[np.float64]] = np.load(_F_PATH).astype(float)

_N1: Final[NDArray[np.float64]] = np.array([0.0, 0.0, 1.0])

# Zeta table radial support (Mpc); the integrand is identically zero outside it.
_LAM_LO: Final[float] = 396.63
_LAM_HI: Final[float] = 2326.61


def _n2_of_cos(cos_gamma: float) -> NDArray[np.float64]:
    c = float(np.clip(cos_gamma, -1.0, 1.0))
    s = float(np.sqrt(max(0.0, 1.0 - c * c)))
    return np.array([s, 0.0, c])


def _G_tensor(cos_gamma: float, lam_v: float) -> NDArray[np.float64]:
    """G_ab(lam_v; gamma) = (M + M^T), M_ab = sum_{mn} F_{a m n} Z_{b m n}.

    Z is the bare equal-time (3,3,3) cumulant at the cosine triple
    (1, cos gamma, cos gamma): two legs co-located on ray n1, one on ray n2.
    """
    n2 = _n2_of_cos(cos_gamma)
    Z = _ds.zeta_tensor(_N1, _N1, n2, float(lam_v))  # triple (1, cosg, cosg)
    M = np.einsum("amn,bmn->ab", _F, Z, optimize=True)
    return M + M.T


def fk_analytic(cos_gamma: float, lam_f: float,
                n_lam: int = 96, n_win: int = 64) -> NDArray[np.float64]:
    """Analytic equal-time FK (3,3) matrix at separation ``cos_gamma``.

    Implements eq. (FK).  The double integral is reduced to two 1-D quadratures
    by swapping the order:

        xi_FK_ab = int dlam_v  W(lam_v) u(lam_v)^4 G_ab(lam_v) * H(lam_v),
        H(lam_v) = int_{lam_v}^{lam_f} W(tau) u(tau)^{-4} dtau,    u = D / D(lam_f)

    (the D(lam_f) normalisation keeps u ~ O(1) so u^{+/-4} never overflows).

    Parameters
    ----------
    cos_gamma : float
        Cosine of the angular separation.
    lam_f : float
        Source affine parameter (analysis-3 baseline 2313.029 Mpc).
    n_lam : int
        Gauss-Legendre nodes for the outer (lam_v) integral.
    n_win : int
        Gauss-Legendre nodes for each LOS window W.
    """
    bg = _ds._builder().bg  # share the cached Background
    lf = float(lam_f)
    lo = _LAM_LO
    Df = float(bg.D(lf))

    # Dense grid for H(lam_v) via a splined antiderivative of W(tau) u(tau)^-4.
    dense = np.linspace(lo, lf, 400)
    W_dense = _ds.order0_window(bg, dense, lf, n_gauss=n_win)
    u_dense = np.asarray(bg.D(dense), float) / Df
    h_integrand = W_dense * u_dense ** (-4)
    h_spline = CubicSpline(dense, h_integrand)
    H_anti = h_spline.antiderivative()
    H_at_lf = float(H_anti(lf))

    # Outer Gauss-Legendre over lam_v in [lo, lf].
    x, w = leggauss(int(n_lam))
    lam_v = 0.5 * (lf - lo) * x + 0.5 * (lf + lo)
    jac = 0.5 * (lf - lo) * w

    W_v = _ds.order0_window(bg, lam_v, lf, n_gauss=n_win)
    u_v = np.asarray(bg.D(lam_v), float) / Df
    H_v = H_at_lf - np.asarray([float(H_anti(l)) for l in lam_v], float)

    out = np.zeros((3, 3), dtype=float)
    for i, lv in enumerate(lam_v):
        G = _G_tensor(cos_gamma, float(lv))
        out += jac[i] * W_v[i] * u_v[i] ** 4 * H_v[i] * G
    return out


# ---------------------------------------------------------------------------
# Anchor against analysis-3 C_corr_op_K_limber_FK (kk)
# ---------------------------------------------------------------------------
def _load_analysis3_fk():
    path = (_ds._SACHS_SFT / "sftwick_outputs" / "2PCF"
            / "C_corr_op_K_limber_FK" / "xi_C_corr_op_K_limber_FK.npz")
    d = np.load(path, allow_pickle=True)  # trusted project artifact
    m = (d["a"] == 0) & (d["b"] == 0)
    xs, ys = d["x"][m], d["y"][m]
    vs = np.asarray(d["value"], float)[m]
    lam_f = float(np.unique(d["t_final"])[0])
    gammas, coss = [], []
    for xi, yi in zip(xs, ys):
        xi = np.asarray(xi, float); yi = np.asarray(yi, float)
        c = float(np.clip(np.dot(xi / np.linalg.norm(xi),
                                 yi / np.linalg.norm(yi)), -1.0, 1.0))
        gammas.append(math.degrees(math.acos(c)) * 60.0)
        coss.append(c)
    order = np.argsort(gammas)
    return (np.asarray(gammas)[order], np.asarray(coss)[order],
            vs[order], lam_f)


def _self_test() -> None:
    print("=" * 78)
    print("fk_analytic self-test  (independent equal-time FK vs analysis-3)")
    print("=" * 78)
    gammas, coss, vs, lam_f = _load_analysis3_fk()
    fk = np.array([fk_analytic(c, lam_f)[0, 0] for c in coss])
    ratio = fk / np.where(vs != 0, vs, np.nan)

    print(f"\n{'gamma[arcmin]':>14} {'fk_analytic':>14} {'analysis3 FK':>14} "
          f"{'ratio':>9}")
    for label_g in (0.5, 1.0, 2.6, 10.0, 44.0, 100.0):
        i = int(np.argmin(np.abs(gammas - label_g)))
        print(f"{gammas[i]:14.3f} {fk[i]:14.5e} {vs[i]:14.5e} {ratio[i]:9.4f}")

    pos = (vs > 0) & (gammas < 200.0)
    med = float(np.median(ratio[pos]))
    spread = float(np.std(ratio[pos]))
    print(f"\n  small-gamma (gamma<200', FK>0): median ratio = {med:.4f} "
          f"std = {spread:.4f}  [{ratio[pos].min():.4f}, {ratio[pos].max():.4f}]")
    print(f"  fk_analytic kk peak = {np.nanmax(fk):.5e}  "
          f"(analysis-3 peak {np.nanmax(vs):.5e})")
    ok = abs(med - 1.0) < 0.20 and spread < 0.10
    print(f"\n  {'PASS' if ok else 'CHECK'} -- the FK expression reproduces "
          f"analysis-3 {'to O(1) with a flat ratio' if ok else 'NOT flat'}.")

    _plot(gammas, fk, vs, ratio)


def _plot(gammas, fk, vs, ratio) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print(f"  [plot skipped: {exc}]")
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.loglog(gammas, np.abs(fk), "o-", ms=3, label="fk_analytic (this module)")
    ax1.loglog(gammas, np.abs(vs), "s--", ms=3, label="analysis-3 FK")
    ax1.set_xlabel("gamma [arcmin]"); ax1.set_ylabel("|FK_kk|")
    ax1.set_title("Equal-time FK anchor (kk)"); ax1.legend()
    ax1.grid(True, which="both", alpha=0.3)
    ax2.semilogx(gammas, ratio, "o-", ms=3)
    ax2.axhline(1.0, color="k", lw=0.8)
    ax2.set_xlabel("gamma [arcmin]"); ax2.set_ylabel("fk_analytic / analysis-3")
    ax2.set_ylim(0.0, 2.0); ax2.set_title("anchor ratio")
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    out_dir = Path(__file__).resolve().parent / "figures"
    out_dir.mkdir(exist_ok=True)
    out = out_dir / "anchor_FK_check.pdf"
    fig.savefig(out); plt.close(fig)
    print(f"  [anchor FK plot -> {out}]")


__all__ = ["fk_analytic"]


if __name__ == "__main__":
    _self_test()
