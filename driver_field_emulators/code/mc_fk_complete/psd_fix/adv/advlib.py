"""Adversarial-review library: a parameter-free EXACT density, tunable variants
of R1/R0/stock, and three quadrature rules for the Order-0 anchor.

The equal-time diagonal C_ab(t,t) served by corr_op is EXACTLY piecewise
quadratic in lambda (verified to 1.7e-15 in a0_structure.py), because the table
is interpolated BILINEARLY.  Therefore
    Sigma2_ab(lam) = d/dlam[D^4 C_ab]/D^4 = dC_ab/dlam + 4 (D'/D) C_ab
is analytic with ZERO free parameters: no finite-difference step, no spline, no
knot choice.  ``Exact`` below implements it.  It is the yardstick that the
attack claims cannot exist.
"""
from __future__ import annotations
import sys, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R
import corr_op_C_callable as CO
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import leggauss

N_COMP = ds.N_COMP
TABLE_LAM = np.load(CO.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)
LAM_F = bgm.LAM_SOURCE_BASELINE
BG = bgm.Background(lam_min=120.0, lam_max=LAM_F + 5.0)


class Exact(R._Base):
    """Parameter-free: exact per-cell quadratic of C on the diagonal."""

    def _density_splines(self, cos_gamma):
        key = round(float(cos_gamma), 12)
        c = self._cache.get(key)
        if c is not None:
            return c
        breaks = R._cells(self.lam_lo, self.lam_hi)
        cells = []
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            # 3 interior points -> exact quadratic (residual checked in a0)
            xs = lo + (hi - lo) * np.array([0.25, 0.5, 0.75])
            F = self._sample_f_C(cos_gamma, xs)          # raw C, no D^4
            u = (xs - lo) / (hi - lo)
            V = np.vander(u, 3)                           # [u^2, u, 1]
            coef = np.linalg.solve(V, F.reshape(3, -1)).reshape(3, N_COMP, N_COMP)
            cells.append((lo, hi, coef))
        self._cache[key] = cells
        return cells

    def _sample_f_C(self, cos_gamma, lam):
        n2 = ds._n2_of_cos(cos_gamma)
        return np.array([CO.C_fn(ds._N1, float(t), n2, float(t)) for t in lam])

    def _cell(self, cells, x, side="right"):
        for lo, hi, coef in cells:
            if side == "right":
                if lo - 1e-12 <= x < hi - 1e-12 or (x >= cells[-1][1] - 1e-12 and hi == cells[-1][1]):
                    return lo, hi, coef
            else:
                if lo + 1e-12 < x <= hi + 1e-12:
                    return lo, hi, coef
        return cells[-1]

    def C_and_dC(self, cos_gamma, lam, side="right"):
        cells = self._density_splines(cos_gamma)
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        lo, hi, coef = self._cell(cells, x, side)
        u = (x - lo) / (hi - lo)
        C = coef[0] * u * u + coef[1] * u + coef[2]
        dC = (2.0 * coef[0] * u + coef[1]) / (hi - lo)
        return C, dC

    def matrix(self, cos_gamma, lam, side="right"):
        C, dC = self.C_and_dC(cos_gamma, lam, side)
        th = float(BG.theta_sa(float(np.clip(lam, self.lam_lo, self.lam_hi))))
        return self._finish(dC + 4.0 * th * C)


class R1Var(R.R1Knots):
    """R1 with a tunable per-cell sampling spacing (production is dx = 12 Mpc)."""

    def __init__(self, *a, dx: float = 12.0, nmin: int = 6, **kw):
        super().__init__(*a, **kw)
        self.dx = float(dx); self.nmin = int(nmin)

    def _density_splines(self, cos_gamma):
        key = round(float(cos_gamma), 12)
        c = self._cache.get(key)
        if c is not None:
            return c
        breaks = R._cells(self.lam_lo, self.lam_hi)
        per_cell = []
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            n = max(self.nmin, int(np.ceil((hi - lo) / self.dx)) + 1)
            xs = np.linspace(lo, hi, n)
            f = self._sample_f(cos_gamma, xs)
            per_cell.append((lo, hi, [[CubicSpline(xs, f[:, a, b])
                                       for b in range(N_COMP)]
                                      for a in range(N_COMP)]))
        self._cache[key] = per_cell
        return per_cell


# ---------------------------------------------------------------------------
# Order-0 anchor with a choice of outer quadrature rule.  G is always
# n_gauss=256 (smooth integrand, converged) so only the outer rule varies.
# ---------------------------------------------------------------------------
def _integ(builder, cos_gamma, la, w):
    sig = np.array([builder.matrix(cos_gamma, float(l))[0, 0] for l in la])
    G = ds.order0_window(builder.bg, la, LAM_F, n_gauss=256)
    return float(np.sum(w * sig * G ** 2))


def rule_gl(lam_lo, lam_f, n):
    x, w = leggauss(int(n))
    return 0.5 * (lam_f - lam_lo) * x + 0.5 * (lam_f + lam_lo), 0.5 * (lam_f - lam_lo) * w


def rule_panel(lam_lo, lam_f, n_per=48):
    """Jump-aware: break at every corr_op table node inside the range."""
    inner = TABLE_LAM[(TABLE_LAM > lam_lo + 1e-9) & (TABLE_LAM < lam_f - 1e-9)]
    edges = np.concatenate([[lam_lo], inner, [lam_f]])
    x, w = leggauss(int(n_per))
    LA, W = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        LA.append(0.5 * (hi - lo) * x + 0.5 * (hi + lo)); W.append(0.5 * (hi - lo) * w)
    return np.concatenate(LA), np.concatenate(W)


def rule_trap(lam_lo, lam_f, n):
    """Uniform trapezoid -- what the MC's own uniform lambda grid realises.
    Endpoints nudged inward by 1e-6 to dodge the R0Reference endpoint bug."""
    la = np.linspace(lam_lo, lam_f, int(n))
    la[0] += 1e-6; la[-1] -= 1e-6
    w = np.gradient(np.linspace(lam_lo, lam_f, int(n)))
    w[0] *= 0.5 * 2 / 2; # gradient already gives trapezoid-like weights
    dl = (lam_f - lam_lo) / (n - 1)
    w = np.full(int(n), dl); w[0] *= 0.5; w[-1] *= 0.5
    return la, w


GAM = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3, 293.9, 596.9])
_d = np.load("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz",
             allow_pickle=True)
G_A3, O0_A3 = np.asarray(_d["g"], float), np.asarray(_d["o0"], float)


def anchor(builder, rule, gam=GAM):
    la, w = rule(builder.lam_lo, LAM_F)
    r = []
    for g in gam:
        cosg = float(np.cos(np.deg2rad(g / 60.0)))
        a3 = float(np.interp(g, G_A3, O0_A3))
        r.append(_integ(builder, cosg, la, w) / a3)
    return np.array(r)
