#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""probe_J2_radial_swap.py -- DIAGNOSTIC A: decouple radial from angular.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

GOAL
----
Feed the PROVEN-correct J_2 orientation kernel (the centroid-projected spin-2
2D-Fourier assembly `natural_components` from stage1_bkappa_phase_reference.py)
TWO radial inputs in turn, holding the angular kernel BYTE-IDENTICAL:

  (a) refB  = B_kappa(l1,l2,l3)   -- the trusted B_kappa phase reference radial shape
  (b) canB  = W_can(l1,l2,l3)     -- the canoes HIGH B_delta convergence-3pt scalar
                                     weight (the deployed zeta_D radial integrand,
                                     PPP response), as a general triangle function.

Both are the SCALAR coefficient that multiplies the SAME orientation kernel, so
passing each through the identical `natural_components` isolates the radial
ell-shape difference. The canoes B_delta differs from B_kappa by the +0.44
log-log ell-slope transfer R(l1,l2,l3)=W_can/B_kappa documented in
outputs/Wcan_vs_Bkappa_decomp.npz.

OUTPUT
------
At gamma=10', phi in {30,60,120} deg, compute the frame-invariant per-config
value for each radial input, form the ratio canB-through-J2 / refB-through-J2,
and check whether it reproduces the target J_2 residual phi-trend
[0.68 @30, 0.43 @60, 2.03 @120]. Report the residual AFTER dividing it out, and
quantify any leftover angular/channel piece.

Both radial shapes share the SAME (u,v,phi) angular measure, (2pi)^-4 norm, and
SAS triangle geometry; only the radial function B(l1,l2,l3) is swapped. So any
phi-dependence in the ratio is 100% attributable to the radial ell-shape.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"
sys.path.insert(0, str(ANALYSIS))

# the PROVEN-correct J_2 assembly + the reference B_kappa radial builder
from stage1_bkappa_phase_reference import (  # noqa: E402
    BkappaLimber, BkappaInterp, natural_components, ARCMIN2RAD, H0_HUNIT, OMEGA_M,
)

# canoes radial integrand
sys.path.insert(0, str(Path("/Users/zzhang/projects/canoes")))
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.sachs.kappa3 import (  # noqa: E402
    chi_of_lambda, build_lightcone, _f2_kernel_numpy,
    _kappa3_limber_response_product_h, _lambda_jacobian,
    _kappa3_radial_density_conversion, _H0_H_PER_MPC,
)


# ---------------------------------------------------------------------------
# 1. canoes scalar weight W_can(l1,l2,l3) -- general triangle, vectorised over
#    the LOS lambda grid (mirrors probe_HIGH_B_ratio_ellscaling.W_can exactly).
# ---------------------------------------------------------------------------
class WcanLimber:
    """canoes HIGH convergence-3pt scalar weight as a function of (l1,l2,l3).

    Identical chain to the deployed zeta_D radial fold (PPP response):
      W_can = INT dlambda K_g^3 * growth^4 * B_delta(k_i) * resp_PPP / chi^4 * rf * h^6
    with k_i = l_i / chi.  This is the canoes B_delta radial shape; it shares the
    (u,v,phi) measure + (2pi)^-4 + orientation kernel with B_kappa.
    """

    def __init__(self, n_z=40):
        self.cosmo = FiducialCosmology()
        self.pk = load_pk_delta(n_s=self.cosmo.n_s)
        zf = np.linspace(0, 5, 6000)
        chif = np.array([chi_of_z(z, Omega_m=self.cosmo.Omega_m, h=self.cosmo.h)
                         for z in zf])
        af = 1.0 / (1.0 + zf)
        lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
        z = np.linspace(5.0 / n_z, 5.0, n_z)
        chi_phys = np.interp(z, zf, chif)
        lam_phys = np.interp(z, zf, lam_f)
        lam_h = lam_phys * self.cosmo.h
        chi_s = float(chif[-1])
        a = 1.0 / (1.0 + z)
        self.lam_phys = lam_phys
        self.Kg = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s
        self.chi_arr = np.maximum(
            chi_of_lambda(self.cosmo, lam_h, convention="project"), 1e-6)
        lc, opz = build_lightcone(self.cosmo, self.chi_arr)
        self.opz = opz
        poisson_amp = -1.5 * self.cosmo.Omega_m * (_H0_H_PER_MPC ** 2) * opz
        djdc = _lambda_jacobian(opz, convention="project")
        rf, _, _ = _kappa3_radial_density_conversion(djdc, "lambda")
        self.rf = rf
        self.growth = np.asarray(lc.M[0] / opz, dtype=np.float64)
        self.resp = np.array([
            _kappa3_limber_response_product_h("PPP", float(poisson_amp[i]),
                                              float(opz[i]))
            for i in range(len(opz))
        ])
        self.h6 = self.cosmo.h ** 6
        # precompute per-shell radial prefactor common to all triangles:
        # K_g^3 * growth^4 * resp / chi^4 * rf * h^6   (b_delta multiplied in per-call)
        self.pref = (self.Kg ** 3) * (self.growth ** 4) * self.resp \
            / (self.chi_arr ** 4) * self.rf * self.h6

    def _scalar(self, l1, l2, l3):
        """W_can for a single triangle (scalar ell magnitudes)."""
        l1 = float(l1); l2 = float(l2); l3 = float(l3)
        if min(l1, l2, l3) <= 0:
            return 0.0
        c12 = (l3 ** 2 - l1 ** 2 - l2 ** 2) / (2 * l1 * l2)
        c13 = (l2 ** 2 - l1 ** 2 - l3 ** 2) / (2 * l1 * l3)
        c23 = (l1 ** 2 - l2 ** 2 - l3 ** 2) / (2 * l2 * l3)
        chi = self.chi_arr
        k1 = l1 / chi; k2 = l2 / chi; k3 = l3 / chi
        p1 = self.pk(k1); p2 = self.pk(k2); p3 = self.pk(k3)
        f12 = _f2_kernel_numpy(k1, k2, np.full_like(chi, c12))
        f13 = _f2_kernel_numpy(k1, k3, np.full_like(chi, c13))
        f23 = _f2_kernel_numpy(k2, k3, np.full_like(chi, c23))
        bdelta = 2.0 * (f12 * p1 * p2 + f13 * p1 * p3 + f23 * p2 * p3)
        integ = self.pref * bdelta
        return float(np.trapezoid(integ, self.lam_phys))

    def grid(self, l1, l2, l3):
        """W_can over an array of triangles, vectorised over LOS.

        l1,l2,l3: 1D arrays of ell magnitudes (one per triangle, same length).
        Returns 1D array of W_can.  Batched in chunks to bound memory at
        O(chunk * n_chi)."""
        l1 = np.asarray(l1, dtype=float).ravel()
        l2 = np.asarray(l2, dtype=float).ravel()
        l3 = np.asarray(l3, dtype=float).ravel()
        n = l1.size
        out = np.zeros(n, dtype=float)
        chi = self.chi_arr[None, :]                      # (1, n_chi)
        lam = self.lam_phys
        pref = self.pref[None, :]                         # (1, n_chi)
        chunk = 8192
        for s in range(0, n, chunk):
            e = min(s + chunk, n)
            a1 = l1[s:e][:, None]; a2 = l2[s:e][:, None]; a3 = l3[s:e][:, None]
            valid = (l1[s:e] > 0) & (l2[s:e] > 0) & (l3[s:e] > 0)
            denom12 = 2 * np.where(a1 * a2 > 0, a1 * a2, 1.0)
            denom13 = 2 * np.where(a1 * a3 > 0, a1 * a3, 1.0)
            denom23 = 2 * np.where(a2 * a3 > 0, a2 * a3, 1.0)
            c12 = (a3 ** 2 - a1 ** 2 - a2 ** 2) / denom12   # (m, 1)
            c13 = (a2 ** 2 - a1 ** 2 - a3 ** 2) / denom13
            c23 = (a1 ** 2 - a2 ** 2 - a3 ** 2) / denom23
            k1 = a1 / chi; k2 = a2 / chi; k3 = a3 / chi     # (m, n_chi)
            p1 = self.pk(k1); p2 = self.pk(k2); p3 = self.pk(k3)
            c12b = np.broadcast_to(c12, k1.shape)
            c13b = np.broadcast_to(c13, k1.shape)
            c23b = np.broadcast_to(c23, k1.shape)
            f12 = _f2_kernel_numpy(k1, k2, c12b)
            f13 = _f2_kernel_numpy(k1, k3, c13b)
            f23 = _f2_kernel_numpy(k2, k3, c23b)
            bdelta = 2.0 * (f12 * p1 * p2 + f13 * p1 * p3 + f23 * p2 * p3)
            integ = pref * bdelta                           # (m, n_chi)
            wv = np.trapezoid(integ, lam, axis=1)
            out[s:e] = np.where(valid, wv, 0.0)
        return out


# ---------------------------------------------------------------------------
# 2. Interpolator wrapper so the inner Fourier loop is fast (same interface as
#    BkappaInterp; tabulate W_can on (ln l2, ln l3, cos23) and trilinear).
# ---------------------------------------------------------------------------
class WcanInterp:
    def __init__(self, wcan_grid, chi_s, Dz, l_min=5e-2, l_max=3.0e4,
                 n_l=96, n_c=64):
        from scipy.interpolate import RegularGridInterpolator as RGI
        self.lnl = np.linspace(np.log(l_min), np.log(l_max), n_l)
        self.cosg = np.linspace(-1.0, 1.0, n_c)
        lr = np.exp(self.lnl)
        L2, L3, C = np.meshgrid(lr, lr, self.cosg, indexing="ij")
        L1 = np.sqrt(np.maximum(L2 ** 2 + L3 ** 2 + 2.0 * L2 * L3 * C, 0.0))
        vals = wcan_grid(L1.ravel(), L2.ravel(), L3.ravel()).reshape(L1.shape)
        self.sign = np.sign(vals)
        self.logabs = np.log(np.abs(vals) + 1e-300)
        self._rgi_log = RGI((self.lnl, self.lnl, self.cosg), self.logabs,
                            bounds_error=False, fill_value=None)
        self._rgi_sgn = RGI((self.lnl, self.lnl, self.cosg), self.sign,
                            bounds_error=False, fill_value=None)
        self.l_min, self.l_max = l_min, l_max
        self.chi_s = chi_s
        self.Dz = Dz

    def __call__(self, l1, l2, l3):
        l1 = np.asarray(l1, dtype=float)
        l2 = np.asarray(l2, dtype=float)
        l3 = np.asarray(l3, dtype=float)
        shape = np.broadcast(l1, l2, l3).shape
        l1b = np.broadcast_to(l1, shape)
        l2b = np.broadcast_to(l2, shape)
        l3b = np.broadcast_to(l3, shape)
        denom = 2.0 * l2b * l3b
        cos23 = np.where(denom > 0,
                         (np.asarray(l1b) ** 2 - l2b ** 2 - l3b ** 2)
                         / np.where(denom > 0, denom, 1.0), 0.0)
        cos23 = np.clip(cos23, -1.0, 1.0)
        lnl2 = np.log(np.clip(l2b, self.l_min, self.l_max))
        lnl3 = np.log(np.clip(l3b, self.l_min, self.l_max))
        pts = np.stack([lnl2.ravel(), lnl3.ravel(), cos23.ravel()], axis=-1)
        logabs = self._rgi_log(pts)
        sgn = np.sign(self._rgi_sgn(pts))
        out = sgn * np.exp(logabs)
        inside = ((l2b >= self.l_min) & (l2b <= self.l_max)
                  & (l3b >= self.l_min) & (l3b <= self.l_max)
                  & (l1b >= self.l_min) & (l1b <= self.l_max)).ravel()
        out = np.where(inside, out, 0.0)
        return out.reshape(shape)


def fmt_c(z):
    return f"({z.real:+.4e}{z.imag:+.4e}j)"


def frame_inv(res):
    return np.sqrt(sum(abs(res[m]) ** 2 for m in range(4)))


def main():
    print("=" * 84)
    print("DIAGNOSTIC A: radial swap through the IDENTICAL J_2 orientation kernel")
    print("=" * 84)

    # --- reference radial shape: B_kappa (exact Limber + interp) ---
    print("\n[1] building reference B_kappa radial shape ...")
    bk_exact = BkappaLimber(n_chi=160)
    refB = BkappaInterp(bk_exact, n_l=96, n_c=64)

    # --- canoes radial shape: W_can (general triangle) + interp ---
    print("[2] building canoes W_can (B_delta HIGH) radial shape ...")
    wc = WcanLimber(n_z=40)
    canB = WcanInterp(wc.grid, chi_s=bk_exact.chi_s, Dz=bk_exact.Dz,
                      n_l=96, n_c=64)

    # sanity: per-triangle ratio at a few equilateral L (should match the
    # documented +0.44 slope / R(L) from Wcan_vs_Bkappa_decomp.npz).
    print("\n  sanity equilateral R = W_can/B_kappa (vs documented +0.44 slope):")
    Ls = [60., 100., 200., 500., 1000.]
    Rs = []
    for L in Ls:
        r = wc._scalar(L, L, L) / float(bk_exact(L, L, L))
        Rs.append(r)
        print(f"    L={L:6.0f}  R={r:9.4f}")
    slope = np.polyfit(np.log(Ls), np.log(np.abs(Rs)), 1)[0]
    print(f"    log-log slope = {slope:+.4f}  (documented +0.44)")

    # --- the three configs at gamma=10' ---
    configs = [(10.0, 30.0), (10.0, 60.0), (10.0, 120.0)]
    PROD = dict(n_l=128, n_ang=128, l_max=3.0e4)

    print("\n[3] feeding BOTH radial shapes through the IDENTICAL J_2 assembly")
    print(f"    (natural_components, settings {PROD})")
    print("-" * 84)

    rows = []
    for (g_am, p_deg) in configs:
        g_rad = g_am * ARCMIN2RAD
        p_rad = np.radians(p_deg)
        res_ref = natural_components(refB, g_rad, p_rad,
                                     n_l=PROD["n_l"], n_ang=PROD["n_ang"],
                                     l_max=PROD["l_max"])
        res_can = natural_components(canB, g_rad, p_rad,
                                     n_l=PROD["n_l"], n_ang=PROD["n_ang"],
                                     l_max=PROD["l_max"])
        fi_ref = frame_inv(res_ref)
        fi_can = frame_inv(res_can)
        ratio_fi = fi_can / fi_ref
        # per-component magnitude ratios (channel split)
        comp_ratio = np.array([abs(res_can[m]) / max(abs(res_ref[m]), 1e-300)
                               for m in range(4)])
        rows.append(dict(g=g_am, phi=p_deg, fi_ref=fi_ref, fi_can=fi_can,
                         ratio_fi=ratio_fi, comp_ratio=comp_ratio,
                         res_ref=res_ref, res_can=res_can))
        print(f"\n  gamma={g_am:.0f}' phi={p_deg:.0f}deg")
        print(f"    fi_ref(refB through J2) = {fi_ref:.4e}")
        print(f"    fi_can(canB through J2) = {fi_can:.4e}")
        print(f"    ratio_fi = canB/refB    = {ratio_fi:.4f}")
        print(f"    per-comp |canB/refB|    = "
              f"[{comp_ratio[0]:.3f}, {comp_ratio[1]:.3f}, "
              f"{comp_ratio[2]:.3f}, {comp_ratio[3]:.3f}]")

    # --- the phi-trend (frame-invariant) and comparison to the target ---
    phis = np.array([r["phi"] for r in rows])
    radial_only = np.array([r["ratio_fi"] for r in rows])  # [phi30, phi60, phi120]
    target = np.array([0.68, 0.43, 2.03])                  # documented J_2 residual

    # normalize both to median (the target IS already the residual-after-median;
    # radial_only is a raw canB/refB ratio whose median scale should match too)
    med = np.median(radial_only)
    radial_norm = radial_only / med

    print("\n" + "=" * 84)
    print("RESULT: phi-trend (frame-invariant) at gamma=10'")
    print("=" * 84)
    print(f"  phi (deg)                 : {phis}")
    print(f"  radial-only canB/refB raw : {radial_only}")
    print(f"  radial-only / median      : {radial_norm}    (median={med:.4f})")
    print(f"  TARGET J_2 residual       : {target}")

    # residual after dividing the radial-only trend OUT of the target
    leftover = target / radial_norm
    leftover_norm = leftover / np.median(leftover)
    print(f"\n  leftover = target / radial_norm        : {leftover}")
    print(f"  leftover / median(leftover)            : {leftover_norm}")
    print(f"  leftover spread (max/min)              : "
          f"{leftover_norm.max() / leftover_norm.min():.4f}")

    # also: does the RAW radial-only ratio reproduce the target directly?
    repro = radial_norm / target
    print(f"\n  radial_norm / target (==1 perfect)     : {repro}")
    print(f"  reproduction spread (max/min)          : "
          f"{repro.max() / repro.min():.4f}")

    # frame-invariant target trend (full) for the record
    print("\n  (the full target phi-trend [0.68,0.43,2.03] reorders the FINAL "
          "residual_j2 [phi30,phi60,phi120])")

    np.savez(
        OUT / "J2_radial_swap.npz",
        phi_deg=phis,
        fi_ref=np.array([r["fi_ref"] for r in rows]),
        fi_can=np.array([r["fi_can"] for r in rows]),
        ratio_fi=radial_only,
        ratio_fi_norm=radial_norm,
        comp_ratio=np.array([r["comp_ratio"] for r in rows]),
        target=target,
        leftover=leftover,
        leftover_norm=leftover_norm,
        repro=repro,
        equilateral_L=np.array(Ls),
        equilateral_R=np.array(Rs),
        equilateral_slope=slope,
        note=("DIAGNOSTIC A: refB=B_kappa and canB=W_can fed through the IDENTICAL "
              "natural_components J_2 assembly at gamma=10', phi=30/60/120. "
              "ratio_fi=canB/refB frame-invariant; radial_norm=ratio_fi/median; "
              "target=[0.68,0.43,2.03] documented J_2 residual; leftover=target/"
              "radial_norm quantifies any non-radial (angular/channel) residual."),
    )
    print(f"\nsaved -> {OUT / 'J2_radial_swap.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
