#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""stage1_bkappa_phase_reference.py
=====================================
Interpreter (ABSOLUTE): /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
(only needs numpy + scipy + astropy; the sft-wick env has all three. PyCCL env
also works.)

THE DECISIVE REFERENCE
----------------------
An INDEPENDENT, transparent cosmic-shear 3PCF natural-component reference built
STRAIGHT FROM the convergence bispectrum B_kappa and the STANDARD flat-sky
spin-2 projection

        gamma~(l) = e^{2 i phi_l} kappa~(l)          [phi_l = polar angle of l]

which IS the definition of the shear 3PCF (Schneider & Lombardi 2003, A&A 408,
829; Schneider, Kilbinger & Lombardi 2005, A&A 431, 9; Sugiyama+2024,
arXiv:2407.01798).  This implementation is independent of BOTH

  * canoes' driving-field spin-2 cumulant zeta_D route  (object 2 under test), AND
  * fastnc's 2DFFTLog multipole machinery                (object 3),

sharing ONLY the physics inputs that must be common: the SAME tree-level SPT
B_delta, the SAME single-plane convergence Limber kernel, the SAME P(k)
(PCAMBz0.txt), the SAME cosmology, z_s = 5, and the CENTROID projection.

DEFINITION (real-space SL03, centroid projection -- no fastnc x-frame needed)
----------------------------------------------------------------------------
Three galaxy positions X1, X2, X3 form the SAS isoceles triangle.  The shear at
each point, projected onto the vertex->centroid direction, is

        g|_j = gamma(X_j) e^{-2 i alpha_j},  alpha_j = polar(centroid - X_j).

The natural components are (s_j^mu = +1 unconjugated leg, -1 conjugated leg):

        Gamma^0 = < g|_1   g|_2   g|_3  >    (s = +,+,+)
        Gamma^1 = < g|_1^* g|_2   g|_3  >    (s = -,+,+)
        Gamma^2 = < g|_1   g|_2^* g|_3  >    (s = +,-,+)
        Gamma^3 = < g|_1   g|_2   g|_3^*>    (s = +,+,-)

In Fourier space, with  <k~(l1) k~(l2) k~(l3)> = (2pi)^2 delta(l1+l2+l3) B_kappa,

  Gamma^mu(X1,X2,X3)
    = INT d^2 l2/(2pi)^2 d^2 l3/(2pi)^2  B_kappa(|l1|,|l2|,|l3|)
        * exp[ i l2.(X2-X1) + i l3.(X3-X1) ]
        * exp[ 2 i ( s1 (phi_{l1}-alpha1) + s2 (phi_{l2}-alpha2)
                     + s3 (phi_{l3}-alpha3) ) ],          l1 = -(l2+l3).

This is the CENTROID projection BY CONSTRUCTION (each leg's spin reference is
its vertex->centroid direction); it does NOT use fastnc's x2cent bridge, so it
is an independent check of the projection too.

B_kappa: IDENTICAL to fastnc's Limber chain (stage1_fastnc_tree.py), in h-units
---------------------------------------------------------------------------------
  chi(z)  = astropy wCDM comoving_distance(z) * h        [Mpc/h]   (bispectrum.py:293)
  g(chi)  = (3/2)(100/c)^2 Om0 (1 - chi/chi_s)           [h-units] (bispectrum.py:369-374)
  weight  = g(chi)^3 * (1/chi) * (1+z)^3                            (bispectrum.py:533)
  B_delta = D(z)^4 [ 2 F2(k1,k2)P(k1)P(k2) + 2 cyc ]    z=0 P(k)
  F2      = 5/7 + (1/2) cosθ (ka/kb+kb/ka) + (2/7) cosθ^2   (TRUE SPT, not F2_eff)
  k_j     = l_j / chi                                    [h/Mpc]   (bispectrum.py:611-612)
  B_kappa(l1,l2,l3) = INT dchi  weight * B_delta(l1/chi,l2/chi,l3/chi; z)

Growth D(z)/D(0): the SAME analytic flat-wCDM growth fastnc uses (stage1
fastnc_tree.build_growth).  z_s = 5 single delta source plane.

UNITS: everything in h-units internally (chi [Mpc/h], k [h/Mpc], P [(Mpc/h)^3]).
B_kappa is dimensionless (convergence bispectrum), so Gamma^mu is dimensionless
-- directly comparable to fastnc's Gamma (which is built the same way).

METHOD / CONVERGENCE
--------------------
The 2D Fourier integral is evaluated in polar coordinates (l2, phi2, l3, phi3):
4 nested quadratures (two radial log-spaced Gauss-like trapezoid, two angular
uniform).  B_kappa is precomputed on an (l1,l2,l3) triangle table is NOT used;
instead B_kappa is evaluated on-the-fly via a fast Limber routine with a cached
chi-grid.  Internal convergence is reported by refining (n_l, n_ang, n_chi).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from astropy.cosmology import wCDM

HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "outputs"
OUTDIR.mkdir(exist_ok=True)

with open(HERE / "stage1_grid.json") as f:
    GRID = json.load(f)

Z_S = GRID["z_s"]
COSMO = GRID["cosmology"]
OMEGA_M = COSMO["Omega_m"]
H = COSMO["h"]
N_S = COSMO["n_s"]
W0 = COSMO["w0"]
OMEGA_DE = 1.0 - OMEGA_M
PK_FILE = GRID["pk_file"]
C_KM_S = 299792.458
H0_HUNIT = 100.0 / C_KM_S  # h-unit H0 (chi in Mpc/h)

ARCMIN2RAD = np.pi / 10800.0


# ---------------------------------------------------------------------------
# 1. Cosmology helpers (match fastnc exactly)
# ---------------------------------------------------------------------------
def _E(z):
    a = 1.0 / (1.0 + z)
    return np.sqrt(OMEGA_M / a**3 + OMEGA_DE * a ** (-3.0 * (1.0 + W0)))


def _growth_unnorm(z):
    val, _ = quad(lambda zp: (1.0 + zp) / _E(zp) ** 3, z, 1000.0)
    return _E(z) * val


def build_growth_spline():
    z_arr = np.linspace(0.0, Z_S, 400)
    D = np.array([_growth_unnorm(z) for z in z_arr])
    D /= _growth_unnorm(0.0)
    return ius(z_arr, D, ext=1)


def build_chi_splines():
    """chi(z) [Mpc/h] from astropy wCDM, identical to fastnc bispectrum.py:293."""
    cosmo = wCDM(H0=100.0 * H, Om0=OMEGA_M, Ode0=OMEGA_DE, w0=W0)
    z = np.linspace(0.0, Z_S, 400)
    chi = cosmo.comoving_distance(z).value * H  # Mpc/h
    z2chi = ius(z, chi, ext=1)
    chi2z = ius(chi, z, ext=1)
    return z2chi, chi2z, float(z2chi(Z_S))


def load_pklin_spline():
    arr = np.loadtxt(PK_FILE)
    k_h, pk_h = arr[:, 0], arr[:, 1]
    kmin, kmax = k_h.min(), k_h.max()
    sp = ius(np.log(k_h), np.log(pk_h), ext=0)

    def pk(k):
        k = np.asarray(k, dtype=float)
        inside = (k >= kmin) & (k <= kmax)
        logk = np.log(np.where(inside, k, kmin))   # avoid extrapolation overflow
        out = np.exp(sp(logk))
        return np.where(inside, out, 0.0)

    return pk, kmin, kmax


def _F2(k1, k2, k3):
    """Standard SPT F2(k1,k2); k3 = |k1+k2|.  cosθ = (k3^2-k1^2-k2^2)/(2 k1 k2).

    Guarded against degenerate k (when an ell magnitude -> 0): F2 is multiplied by
    P(k_a)P(k_b) in B_delta and P(0)=0 (outside the table) kills those cells, so we
    just return a finite placeholder where k1 or k2 vanish.
    """
    k1s = np.where(k1 > 0, k1, 1.0)
    k2s = np.where(k2 > 0, k2, 1.0)
    cos_t = (k3**2 - k1s**2 - k2s**2) / (2.0 * k1s * k2s)
    f2 = 5.0 / 7.0 + 0.5 * cos_t * (k1s / k2s + k2s / k1s) + (2.0 / 7.0) * cos_t**2
    return np.where((k1 > 0) & (k2 > 0), f2, 0.0)


# ---------------------------------------------------------------------------
# 2. B_kappa Limber routine (vectorised, identical kernel to fastnc)
# ---------------------------------------------------------------------------
class BkappaLimber:
    def __init__(self, n_chi=160):
        self.pk, self.kmin, self.kmax = load_pklin_spline()
        self.Dz = build_growth_spline()
        self.z2chi, self.chi2z, self.chi_s = build_chi_splines()
        # LOS chi grid 0..chi_s (avoid chi=0 singularity in 1/chi)
        self.n_chi = n_chi
        zl = np.linspace(1e-4, Z_S, n_chi)
        self.zl = zl
        self.chil = self.z2chi(zl)
        # convergence kernel g(chi) (h-units), single source plane z_s
        pref = 1.5 * (H0_HUNIT**2) * OMEGA_M
        self.g = pref * (1.0 - self.chil / self.chi_s)  # comoving lensing efficiency
        self.weight = (self.g**3) * (1.0 / self.chil) * (1.0 + zl) ** 3
        self.D4 = self.Dz(zl) ** 4

    def __call__(self, l1, l2, l3):
        """B_kappa(l1,l2,l3) for arbitrary-shape ell magnitude arrays (broadcast).

        Returns array of shape broadcast(l1,l2,l3).  Integrates over the cached
        chi grid (trapezoid)."""
        l1 = np.asarray(l1, dtype=float)
        l2 = np.asarray(l2, dtype=float)
        l3 = np.asarray(l3, dtype=float)
        shape = np.broadcast(l1, l2, l3).shape
        l1 = np.broadcast_to(l1, shape).ravel()
        l2 = np.broadcast_to(l2, shape).ravel()
        l3 = np.broadcast_to(l3, shape).ravel()
        n_ell = l1.size
        chi = self.chil[None, :]                      # (1, n_chi)
        k1 = l1[:, None] / chi
        k2 = l2[:, None] / chi
        k3 = l3[:, None] / chi
        p1 = self.pk(k1)
        p2 = self.pk(k2)
        p3 = self.pk(k3)
        bdelta = (2.0 * _F2(k1, k2, k3) * p1 * p2
                  + 2.0 * _F2(k2, k3, k1) * p2 * p3
                  + 2.0 * _F2(k3, k1, k2) * p3 * p1)
        bdelta *= self.D4[None, :]
        integrand = self.weight[None, :] * bdelta     # (n_ell, n_chi)
        bk = np.trapezoid(integrand, self.chil, axis=1)
        return bk.reshape(shape)


class BkappaInterp:
    """Fast B_kappa(l1,l2,l3) via a 3D spline in (ln l2, ln l3, cos23).

    B_kappa depends only on the three ell magnitudes; with l1^2 = l2^2+l3^2+
    2 l2 l3 cos23, the natural reduced coordinates are (l2, l3, cos23).  We tabulate
    the exact Limber B_kappa on a log-log-linear grid and trilinearly interpolate
    inside the inner Fourier loop -- O(1e3) exact evaluations replace O(1e8).
    Accuracy is validated against the exact routine in the convergence study.
    """

    def __init__(self, bk_exact, l_min=5e-2, l_max=3.0e4, n_l=96, n_c=64):
        from scipy.interpolate import RegularGridInterpolator as RGI
        self.bk_exact = bk_exact
        self.lnl = np.linspace(np.log(l_min), np.log(l_max), n_l)
        self.cosg = np.linspace(-1.0, 1.0, n_c)
        lr = np.exp(self.lnl)
        L2, L3, C = np.meshgrid(lr, lr, self.cosg, indexing="ij")
        L1 = np.sqrt(np.maximum(L2**2 + L3**2 + 2.0 * L2 * L3 * C, 0.0))
        vals = bk_exact(L1, L2, L3)
        # store signed-log for smooth interpolation (B_kappa>0 for tree SPT here,
        # but guard with sign anyway)
        self.sign = np.sign(vals)
        self.logabs = np.log(np.abs(vals) + 1e-300)
        self._rgi_log = RGI((self.lnl, self.lnl, self.cosg), self.logabs,
                            bounds_error=False, fill_value=None)
        self._rgi_sgn = RGI((self.lnl, self.lnl, self.cosg), self.sign,
                            bounds_error=False, fill_value=None)
        self.l_min, self.l_max = l_min, l_max
        self.chi_s = bk_exact.chi_s
        self.Dz = bk_exact.Dz

    def __call__(self, l1, l2, l3):
        l1 = np.asarray(l1, dtype=float)
        l2 = np.asarray(l2, dtype=float)
        l3 = np.asarray(l3, dtype=float)
        shape = np.broadcast(l1, l2, l3).shape
        l1b = np.broadcast_to(l1, shape)
        l2b = np.broadcast_to(l2, shape)
        l3b = np.broadcast_to(l3, shape)
        # reduced cos23 from the three magnitudes
        denom = 2.0 * l2b * l3b
        cos23 = np.where(denom > 0,
                         (np.asarray(l1b) ** 2 - l2b**2 - l3b**2) / np.where(denom > 0, denom, 1.0),
                         0.0)
        cos23 = np.clip(cos23, -1.0, 1.0)
        lnl2 = np.log(np.clip(l2b, self.l_min, self.l_max))
        lnl3 = np.log(np.clip(l3b, self.l_min, self.l_max))
        pts = np.stack([lnl2.ravel(), lnl3.ravel(), cos23.ravel()], axis=-1)
        logabs = self._rgi_log(pts)
        sgn = np.sign(self._rgi_sgn(pts))
        out = sgn * np.exp(logabs)
        # zero out cells outside the tabulated ell range (no power) for any leg
        inside = ((l2b >= self.l_min) & (l2b <= self.l_max)
                  & (l3b >= self.l_min) & (l3b <= self.l_max)
                  & (l1b >= self.l_min) & (l1b <= self.l_max)).ravel()
        out = np.where(inside, out, 0.0)
        return out.reshape(shape)


# ---------------------------------------------------------------------------
# 3. Geometry: SAS isoceles vertices + centroid reference angles
# ---------------------------------------------------------------------------
def triangle_geometry(gamma_rad, phi_rad):
    """Return X1,X2,X3 (apex at origin) and d2,d3,alpha1,alpha2,alpha3.

    SAS isoceles: equal sides t1=t2=gamma meet at APEX (vertex 3) with opening
    angle phi.  Base vertices 1,2; t3 = 2 gamma sin(phi/2).  (Matches the
    symbolic derivation in spin2_centroid_direct.wl.)
    """
    g = gamma_rad
    X3 = np.array([0.0, 0.0])
    X1 = g * np.array([np.cos(phi_rad / 2.0),  np.sin(phi_rad / 2.0)])
    X2 = g * np.array([np.cos(phi_rad / 2.0), -np.sin(phi_rad / 2.0)])
    cen = (X1 + X2 + X3) / 3.0
    alpha1 = np.arctan2((cen - X1)[1], (cen - X1)[0])
    alpha2 = np.arctan2((cen - X2)[1], (cen - X2)[0])
    alpha3 = np.arctan2((cen - X3)[1], (cen - X3)[0])
    d2 = X2 - X1
    d3 = X3 - X1
    return X1, X2, X3, d2, d3, alpha1, alpha2, alpha3


# ---------------------------------------------------------------------------
# 4. The 2D Fourier-plane natural-component integral
# ---------------------------------------------------------------------------
# sign patterns (s1,s2,s3) per natural component mu
_SIGNS = {0: (1, 1, 1), 1: (-1, 1, 1), 2: (1, -1, 1), 3: (1, 1, -1)}


def natural_components(bk, gamma_rad, phi_rad,
                       n_l=220, n_ang=128, l_min=1e-1, l_max=2.0e4):
    """Compute Gamma^0..3 at one (gamma, phi) by the 2D Fourier integral.

    Polar quadrature: l2,l3 log-spaced (radial), phi2,phi3 uniform (angular).
    l1 = -(l2+l3), phi_{l1} = polar(l1).  Returns dict mu->complex.

    Memory-safe: the (l2, phi2) inner plane is fully vectorised; we loop over the
    outer (l3, phi3) pairs and accumulate, so peak memory is O(n_l * n_ang), NOT
    O(n_l^2 * n_ang^2).
    """
    X1, X2, X3, d2, d3, a1, a2, a3 = triangle_geometry(gamma_rad, phi_rad)

    # radial nodes (log) + log measure  int f d^2 l = int f l dl dphi = int f l^2 dphi dln l
    lr = np.geomspace(l_min, l_max, n_l)
    dlnl = np.gradient(np.log(lr))
    wl = lr * lr * dlnl                       # l^2 dln l = l dl
    pa = np.linspace(0.0, 2.0 * np.pi, n_ang, endpoint=False)
    dpa = 2.0 * np.pi / n_ang

    # inner (l2, phi2) plane  (n_l, n_ang)
    L2 = lr[:, None]
    P2 = pa[None, :]
    l2x = L2 * np.cos(P2)
    l2y = L2 * np.sin(P2)
    w2 = (wl[:, None] * dpa)                   # inner measure (n_l, n_ang)
    pos2 = l2x * d2[0] + l2y * d2[1]           # l2.d2
    # accumulators per mu
    acc = {mu: 0.0 + 0.0j for mu in _SIGNS}
    norm = 1.0 / (2.0 * np.pi) ** 4

    for i3 in range(n_l):
        l3 = lr[i3]
        w3r = wl[i3] * dpa
        for j3 in range(n_ang):
            p3 = pa[j3]
            l3x = l3 * np.cos(p3)
            l3y = l3 * np.sin(p3)
            l1x = -(l2x + l3x)
            l1y = -(l2y + l3y)
            L1 = np.sqrt(l1x**2 + l1y**2)
            phi_l1 = np.arctan2(l1y, l1x)
            Bk = bk(L1, L2, l3)               # (n_l, n_ang)
            pos = pos2 + (l3x * d3[0] + l3y * d3[1])
            base = Bk * np.exp(1j * pos) * w2 * w3r * norm
            for mu, (s1, s2, s3) in _SIGNS.items():
                spin = 2.0 * (s1 * (phi_l1 - a1)
                              + s2 * (P2 - a2)
                              + s3 * (p3 - a3))
                acc[mu] += np.sum(base * np.exp(1j * spin))
    return {mu: complex(acc[mu]) for mu in _SIGNS}


# ---------------------------------------------------------------------------
# 5. Driver
# ---------------------------------------------------------------------------
def fmt_c(z):
    return f"({z.real:+.4e}{z.imag:+.4e}j)"


def main():
    print("=" * 78)
    print("DECISIVE B_kappa-phase REFERENCE: shear 3PCF natural components")
    print(f"z_s={Z_S}  Om={OMEGA_M:.4f}  h={H}  projection=cent (by construction)")
    print("=" * 78)

    # configs: primary first, then 2-3 more for structure
    configs = [
        (10.0, 60.0),   # PRIMARY
        (10.0, 30.0),
        (10.0, 120.0),
        (50.0, 60.0),
    ]

    # B_kappa: exact Limber chi-integral, tabulated + trilinearly interpolated for
    # the inner Fourier loop (interp validated to ~1e-3 median rel err vs exact).
    print("\n--- building B_kappa (exact Limber + interpolator) ---")
    bk_exact = BkappaLimber(n_chi=160)
    bk = BkappaInterp(bk_exact, n_l=96, n_c=64)
    print(f"chi_s = {bk.chi_s:.2f} Mpc/h,  D(z_s)={bk.Dz(Z_S):.4f}")
    print(f"B_kappa(100,100,100)={bk_exact(100.,100.,100.):.4e}  "
          f"B_kappa(1000,1000,1732)={bk_exact(1000.,1000.,1732.):.4e}")
    # interp-vs-exact spot check
    rng = np.random.default_rng(0)
    errs = []
    for _ in range(300):
        l2 = 10 ** rng.uniform(0.5, 4.0); l3 = 10 ** rng.uniform(0.5, 4.0)
        c = rng.uniform(-1, 1); l1 = np.sqrt(max(l2**2 + l3**2 + 2 * l2 * l3 * c, 0.0))
        if l1 < 1:
            continue
        ex = bk_exact(l1, l2, l3)
        if abs(ex) > 1e-25:
            errs.append(abs(float(bk(l1, l2, l3)) - ex) / abs(ex))
    errs = np.array(errs)
    print(f"interp vs exact B_kappa rel err: median={np.median(errs):.2e} "
          f"p90={np.percentile(errs,90):.2e} max={errs.max():.2e}")

    # ----- internal convergence study at the PRIMARY point -----------------
    print("\n--- internal convergence (PRIMARY gamma=10', phi=60deg) ---")
    g0_p = 10.0 * ARCMIN2RAD
    p0_p = np.radians(60.0)
    conv_settings = [
        dict(n_l=64, n_ang=64, l_max=2.0e4),
        dict(n_l=96, n_ang=96, l_max=2.0e4),
        dict(n_l=128, n_ang=128, l_max=3.0e4),
    ]
    last = None
    for s in conv_settings:
        res = natural_components(bk, g0_p, p0_p,
                                 n_l=s["n_l"], n_ang=s["n_ang"], l_max=s["l_max"])
        fi = np.sqrt(sum(abs(res[m]) ** 2 for m in range(4)))
        rd = ""
        if last is not None:
            rd = f"  d(frame-inv)={abs(fi-last)/max(last,1e-300):.2e}"
        print(f"  n_l={s['n_l']:3d} n_ang={s['n_ang']:3d} "
              f"lmax={s['l_max']:.0e} | frame-inv={fi:.4e}{rd}")
        print(f"      G0={fmt_c(res[0])} G1={fmt_c(res[1])} "
              f"G2={fmt_c(res[2])} G3={fmt_c(res[3])}")
        last = fi

    # production setting (converged: 128/128 within ~3% of 96/96)
    PROD = dict(n_l=128, n_ang=128, l_max=3.0e4)
    print(f"\n--- production setting {PROD} ---")

    results = {}
    for (g_am, p_deg) in configs:
        g_rad = g_am * ARCMIN2RAD
        p_rad = np.radians(p_deg)
        res = natural_components(bk, g_rad, p_rad,
                                 n_l=PROD["n_l"], n_ang=PROD["n_ang"],
                                 l_max=PROD["l_max"])
        results[(g_am, p_deg)] = res
        fi = np.sqrt(sum(abs(res[m]) ** 2 for m in range(4)))
        print(f"\n[ref] gamma={g_am:.0f}' phi={p_deg:.0f}deg  frame-inv={fi:.4e}")
        for m in range(4):
            print(f"    Gamma^{m} = {fmt_c(res[m])}  |.|={abs(res[m]):.4e}")

    # ----- save npz --------------------------------------------------------
    gam_arr = np.array([c[0] for c in configs])
    phi_arr = np.array([c[1] for c in configs])
    G = {m: np.array([results[c][m] for c in configs], dtype=np.complex128)
         for m in range(4)}
    outpath = OUTDIR / "stage1_bkappa_phase_ref.npz"
    np.savez(
        outpath,
        gamma_arcmin=gam_arr, phi_deg=phi_arr,
        Gamma0=G[0], Gamma1=G[1], Gamma2=G[2], Gamma3=G[3],
        projection="cent",
        prod_setting=json.dumps({**PROD, "n_chi": 160, "interp_nl": 96, "interp_nc": 64}),
        z_s=Z_S, Omega_m=OMEGA_M, h=H, n_s=N_S,
        convention_note=(
            "B_kappa-phase DECISIVE reference: shear 3PCF natural components from "
            "the 2D Fourier integral of B_kappa with centroid-projected spin-2 "
            "phases e^{2i s_j (phi_{l_j}-alpha_j)}, alpha_j=polar(centroid-X_j). "
            "B_kappa = fastnc-identical single-plane Limber: g=(3/2)(100/c)^2 Om "
            "(1-chi/chi_s), weight=g^3/chi (1+z)^3, B_delta=D^4[2 F2 PP+2cyc] true "
            "SPT F2, PCAMBz0.txt z=0 P(k), h-units. Independent of canoes zeta_D AND "
            "fastnc 2DFFTLog. Configs: (10,60)=PRIMARY,(10,30),(10,120),(50,60)."),
    )
    print(f"\n[ref] saved -> {outpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
