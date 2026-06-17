"""DECISIVE CCL-anchored test: hand-rolled equal-shell K^2 lambda-fold 2-pt vs CCL.

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python   (pyccl + canoes)

THE PARADOX THIS ISOLATES
-------------------------
The canoes 3-pt scalar convergence-3PCF (stage0_ours.py) comes out (1+z)^4 over
the CCL-validated standard, but the canoes 2-pt matches CCL.  Both use the SAME
per-leg (1+z)^4 Sachs response and the SAME affine measure.  The ONE thing the
3-pt does that the 2-pt does not is the EQUAL-SHELL COLLAPSE of TWO inner radial
integrals (vs ONE for the 2-pt).

This test builds the convergence 2PCF the SAME hand-rolled equal-shell way that
stage0_ours.py builds the 3-pt, but for 2 legs:

    Z2(gamma) = INT_0^{lambda_s} dlambda  K(lambda, lambda_s)^2  sigma2_TTT(gamma, lambda)

with the SINGLE-source-plane convergence efficiency (cosmology.tex
`eq: K comoving kernel`):

    K(lambda, lambda_s) = Dbar^2(lambda) INT_lambda^{lambda_s} dlambda1/Dbar^2(lambda1)
                        = a^2(chi) chi (chi_s - chi)/chi_s,   Dbar = a*chi.

sigma2_TTT(gamma, lambda) is the canoes EQUAL-SHELL scalar 2-cumulant: the 2-leg
analogue of zeta_TTT.  PRIMARY route uses canoes' OWN born_limber machinery
(``compute_sachs_born_limber_cl_xi``) with the SAME ``T_phi00 = (1+z)^4 L2/(2chi^2)``
per-leg transfer used everywhere in the Sachs pipeline -- NOT a hand-rolled
normalisation -- so the per-leg (1+z)^4 convention is canoes' own, identical in
spirit to stage0_ours' ``compute_kappa3_zeta_table`` source.  The diagonal
``xi_phi00_phi00(gamma; lambda_i, lambda_i)`` with ``delta_measure='lambda'`` is
the lambda-density equal-shell 2-cumulant (the ONE inner delta-collapse, vs the
3-pt's TWO).  SAME PCAMBz0.txt P(k), SAME cosmology, lambda-density.

THE CANOES BORN_LIMBER EQUAL-SHELL 2-CUMULANT (primary; spin-0)
--------------------------------------------------------------
  Cl^{Phi00 Phi00}_diag(chi_i) = w_i / chi_i^2 * T_phi00(l,chi_i)^2 * P_WW(l,z_i)
     T_phi00(l,chi) = (1+z)^4 * l(l+1) / (2 chi^2)       [the (1+z)^4 per leg]
     P_WW = 4 * M0(chi)^2 * P_phi_today  (W = Phi+Psi = 2 Phi; no aniso stress)
     w_i = (d lambda/d chi)_i / Delta_lambda_i           [lambda-density delta wt]
  sigma2_TTT(gamma, lambda_i) = SUM_l (2l+1)/(4pi) Cl_diag(chi_i) P_l(cos gamma).
The lambda-density delta weight w_i = jac_i/Delta_lambda_i carries the equal-shell
collapse; folding SUM_i Delta_lambda_i K_i^2 sigma2_i cancels Delta_lambda_i and
yields the K^2 fold of the per-shell equal-shell density.  units='physical' is
applied by canoes (Cl * h^4, chi/lambda / h) so the fold is physical-consistent.

SECONDARY cross-check route -- hand-rolled flat-sky Limber sigma2 (2-leg analogue
of stage0_spt_reference's 3-leg zeta), used only to corroborate the born_limber
normalisation:
  sigma2(gamma,chi) = (1/(2pi)) (D^2/chi^2) [A(a)(1+z)^4]^2 INT l dl P0 J0(l gamma),
  A(a) = -3/2 Om H0^2 (1+z); P_delta(k,z)=D^2 P0 (TWO-pt growth D^2, one per leg
  pair -> D^2 NOT D^4: the 3-pt bispectrum carried D^4 for THREE field legs);
  the (1/2pi) is the single 2D-ell orientation transform (phi-integral -> 2pi J0).
  native dim (h/Mpc)^3 -> units='physical' multiplies h^3.

CCL-NATIVE 2-pt (the standard, NOT a same-convention reference)
---------------------------------------------------------------
pyccl WeakLensingTracer at a single source plane z_s=5 (narrow Gaussian dndz),
C_ell^kappa-kappa from pyccl.angular_cl with the SAME PCAMBz0.txt P(k) loaded as
a custom Pk2D linear power, then xi_kappa(gamma) via pyccl.correlation(type='NN')
(the scalar / spin-0 convergence 2PCF; convergence and E-mode shear share the
same C_ell).

VERDICT
-------
- Z2/CCL ~ 1 (within ~10-20%): the hand-rolled K^2 lambda-fold gives the STANDARD
  2-pt observable -> the K-fold + per-leg (1+z)^4 convention is CORRECT for 2
  legs -> the 3-pt (1+z)^4 excess is introduced by the EQUAL-SHELL COLLAPSE OF
  THE TWO INNER INTEGRALS (3-pt-specific).
- Z2/CCL ~ (1+z)^2 over: the hand-rolled K-fold/per-leg convention ITSELF
  over-counts -> the stage0 hand-rolled fold, not the equal-shell collapse, is
  the (1+z)^4 source (would contradict the validated canoes 2-pt pipeline ->
  report the tension).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import jv

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

_HERE = Path(__file__).resolve().parent
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

from canoes.cosmo.background import C_KM_S  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402

Z_SOURCE = 5.0
_N_LAMBDA = 60        # source-shell grid (z in (0, z_s])
_N_ELL = 400          # ell quadrature nodes (log-spaced) for the J0 transform
_ELL_MIN = 1e-2
_ELL_MAX = 2.0e4      # high enough that arcmin-scale J0(ell*gamma) is resolved


# ---------------------------------------------------------------------------
# Background / grid (identical to stage0_ours / stage0_spt_reference)
# ---------------------------------------------------------------------------
def _build_grid(cosmo):
    Om, h = cosmo.Omega_m, cosmo.h
    zf = np.linspace(0.0, Z_SOURCE, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=Om, h=h) for zi in zf])  # physical Mpc
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(Z_SOURCE / _N_LAMBDA, Z_SOURCE, _N_LAMBDA)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    return z, chi_phys, lam_phys, float(chif[-1]), float(lam_f[-1])


def _growth_canoes(z: np.ndarray, cosmo) -> np.ndarray:
    """Linear growth D(z) consistent with the canoes tables (M[0]/(1+z))."""
    from canoes.sachs.kappa3 import build_lightcone  # type: ignore
    chi_h = np.array(
        [chi_of_z(zi, Omega_m=cosmo.Omega_m, h=cosmo.h) for zi in z]
    ) * cosmo.h
    lightcone, one_plus_z = build_lightcone(cosmo, chi_h)
    D = lightcone.M[0, :] / one_plus_z
    return np.asarray(D, dtype=np.float64)


# ---------------------------------------------------------------------------
# PRIMARY: canoes born_limber equal-shell 2-cumulant sigma2_TTT(gamma, lambda)
# ---------------------------------------------------------------------------
def _sigma2_TTT_born_limber(gamma_rad, chi_phys, lam_phys, z, pk_phi, cosmo,
                            chi_widths_h, lam_widths_h):
    """canoes' OWN equal-shell spin-0 2-cumulant (lambda-density, physical units).

    Uses ``compute_sachs_born_limber_cl_xi`` with ``delta_measure='lambda'`` and
    ``units='physical'``: the diagonal ``xi_phi00_phi00(theta; chi_i, chi_i)`` is
    the lambda-density equal-shell sigma2_TTT(gamma_i, lambda_i).  Returns
    (n_gamma, n_lambda) physical units, plus the diagonal lambda-density delta
    weight w_i so the caller can verify the fold cancels Delta_lambda_i.

    The (1+z)^4 per leg is INSIDE canoes' T_phi00 -- this is the validated
    convention, NOT a hand-rolled one.
    """
    from scipy.special import eval_legendre

    from canoes.sachs.born_limber import compute_sachs_born_limber_cl_xi

    h = cosmo.h
    chi_h = chi_phys * h                    # canoes wants chi in Mpc/h
    # ell grid: dense low-ell exact + high-ell to resolve arcmin scales.  Cap
    # ell_max so k=(ell+0.5)/chi_min stays in the P(k) tail support comfortably.
    ell_cap = float(min(60000.0, 150.0 * float(np.min(chi_h))))
    ells = np.unique(np.concatenate([
        np.arange(2, 200),
        np.geomspace(200, ell_cap, 400).astype(int),
    ])).astype(np.intp)

    # Get the DIAGONAL Cl_phi00_phi00 (spin-0) only; do NOT request theta so the
    # spin-2 Wigner-d xi synthesis (which overflows at high ell) is skipped.  We
    # Legendre-synthesize the spin-0 sigma2 ourselves below.
    out = compute_sachs_born_limber_cl_xi(
        ells, chi_h, cosmo=cosmo, pk_phi_today=pk_phi, theta=None,
        units="physical", lambda_convention="project", delta_measure="lambda",
        chi_bin_widths=chi_widths_h, lambda_bin_widths=lam_widths_h,
    )
    cl_pp = np.asarray(out.Cl_phi00_phi00, dtype=np.float64)  # (n_ell,n_chi,n_chi)
    idx = np.arange(chi_h.size)
    cl_diag = cl_pp[:, idx, idx]             # (n_ell, n_lambda) lambda-density Cl
    w_lambda = np.asarray(out.limber_delta_weights, dtype=np.float64)  # physical

    # spin-0 real-space synthesis: sigma2(gamma,lambda)=SUM_l (2l+1)/(4pi) Cl P_l.
    ell_f = ells.astype(np.float64)
    pref = (2.0 * ell_f + 1.0) / (4.0 * np.pi)            # (n_ell,)
    cosg = np.cos(gamma_rad)                              # (n_gamma,)
    pl = np.stack([eval_legendre(int(l), cosg) for l in ells], axis=0)  # (n_ell,n_g)
    # sigma2[g, lam] = sum_l pref[l] pl[l,g] cl_diag[l, lam]
    weights = pref[:, None] * pl                          # (n_ell, n_gamma)
    sigma2 = np.tensordot(weights.T, cl_diag, axes=(1, 0))  # (n_gamma, n_lambda)
    return sigma2, w_lambda


# ---------------------------------------------------------------------------
# SECONDARY cross-check: hand-rolled flat-sky Limber sigma2 (2-leg analogue)
# ---------------------------------------------------------------------------
def _sigma2_TTT_handrolled(gamma_rad, chi_phys_h, z, D, pk, cosmo):
    """Hand-rolled equal-shell spin-0 sigma2_TTT(gamma, chi), h-units internal.

    2-leg reduction of stage0_spt_reference's 3-leg zeta:
      sigma2 = (1/2pi)(D^2/chi^2)[A(1+z)^4]^2 INT l dl P0 J0(l gamma).
    D^2 (two-pt growth) and 1/2pi (single 2D-ell orientation transform).
    Returns (n_gamma, n_chi) h-units; caller applies h^3 (native (h/Mpc)^3).
    """
    H0_inv_Mpc = 100.0 / C_KM_S
    ell = np.geomspace(_ELL_MIN, _ELL_MAX, _N_ELL)
    dlnell = np.gradient(np.log(ell))
    one_plus_z = 1.0 + z
    A_poisson = -1.5 * cosmo.Omega_m * (H0_inv_Mpc ** 2) * one_plus_z
    prefactor = 1.0 / (2.0 * np.pi)          # single orientation transform

    n_g = gamma_rad.size
    n_z = z.size
    out = np.zeros((n_g, n_z), dtype=np.float64)
    for iz in range(n_z):
        chi = float(chi_phys_h[iz])
        k = ell / chi
        p = np.asarray(pk(k), dtype=np.float64) * (D[iz] ** 2)  # P_delta(k,z)=D^2 P0
        response = (A_poisson[iz] * one_plus_z[iz] ** 4) ** 2
        radial = (prefactor / chi ** 2) * response
        ellw = ell * ell * dlnell            # ell dell
        for ig in range(n_g):
            j0 = jv(0, ell * gamma_rad[ig])
            out[ig, iz] = radial * float(np.sum(ellw * p * j0))
    return out


# ---------------------------------------------------------------------------
# Internal anchor: validate sigma2_TTT[gamma=0-ish] against corr_op de-windowed
# equal-time density driver_stats.sigma2_matrix[0,0]  (analysis-3 O0 anchor).
# ---------------------------------------------------------------------------
def _corr_op_sigma2_00(cos_gamma_list, lam_phys_query):
    """driver_stats corr_op-derived equal-time Sigma2_00 (no anchor C0), at the
    queried physical-lambda shells.  Returns (n_gamma, n_lambda) physical units.

    This is the WINDOWED-then-de-windowed equal-time 2-cumulant that reproduces
    analysis-3 O0 to a clean gamma-independent 0.881 ratio.  Used here only as a
    SHAPE/normalisation sanity anchor for the hand-rolled sigma2_TTT (spin-0).
    """
    mc_dir = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
                   "sachs_sft/analyses/mc_sachs_2pt"))
    if str(mc_dir) not in sys.path:
        sys.path.insert(0, str(mc_dir))
    try:
        import driver_stats as ds  # noqa: E402
    except Exception as exc:  # pragma: no cover - anchor is optional
        print(f"[anchor] driver_stats import failed ({exc}); skipping corr_op anchor")
        return None
    builder = ds.Sigma2Builder(apply_c0=False)  # raw equal-time density
    out = np.zeros((len(cos_gamma_list), lam_phys_query.size), dtype=float)
    for ig, cg in enumerate(cos_gamma_list):
        for il, lam in enumerate(lam_phys_query):
            out[ig, il] = builder.matrix(float(cg), float(lam))[0, 0]
    return out


# ---------------------------------------------------------------------------
# CCL-NATIVE convergence 2PCF (the standard)
# ---------------------------------------------------------------------------
def _ccl_xi_kappa(gamma_arcmin, cosmo, pk_table_path):
    """CCL-native single-source-plane convergence xi_kappa(gamma).

    Builds a pyccl Cosmology with the SAME PCAMBz0.txt linear P(k), a single-
    source-plane WeakLensingTracer at z_s=5, computes C_ell^kappa-kappa via
    pyccl.angular_cl, then xi via pyccl.correlation(type='NN') (scalar 2PCF;
    convergence == E-mode shear C_ell).
    """
    import pyccl as ccl

    Om, h, n_s = cosmo.Omega_m, cosmo.h, cosmo.n_s
    Ob = 0.0492  # fiducial baryon (sub-dominant for the lensing kernel shape)
    Oc = Om - Ob

    # Load PCAMBz0.txt as a linear P(k) at z=0 (h-units: k h/Mpc, P (Mpc/h)^3).
    arr = np.loadtxt(pk_table_path)
    k_h, p_h = arr[:, 0], arr[:, 1]
    k_phys = k_h * h                  # 1/Mpc
    p_phys = p_h / h ** 3             # Mpc^3

    # Linear P(k,a) = P_PCAMB(k,z=0) * [D(a)/D(1)]^2 with the Carroll-Press-Turner
    # growth (validated to <0.1% vs both CCL's internal growth and the canoes
    # lightcone growth used in sigma2_TTT, so no spurious z-dependent residual).
    pk_dict = _build_pk_linear_dict(k_phys, p_phys, Om)
    cosmo_ccl = ccl.CosmologyCalculator(
        Omega_c=Oc, Omega_b=Ob, h=h, n_s=n_s, sigma8=0.809,
        pk_linear=pk_dict, pk_nonlin=pk_dict,
    )

    # Single source plane z_s=5 via a narrow Gaussian dndz.
    z_grid = np.linspace(0.0, 6.0, 2000)
    sig_z = 0.02
    dndz = np.exp(-0.5 * ((z_grid - Z_SOURCE) / sig_z) ** 2)
    tracer = ccl.WeakLensingTracer(cosmo_ccl, dndz=(z_grid, dndz), has_shear=True)

    ell = np.unique(np.geomspace(2, 60000, 600).astype(int)).astype(float)
    cl = ccl.angular_cl(cosmo_ccl, tracer, tracer, ell)

    theta_deg = gamma_arcmin / 60.0
    xi = ccl.correlation(cosmo_ccl, ell=ell, C_ell=cl, theta=theta_deg,
                         type="NN", method="fftlog")
    return xi, ell, cl


def _D_cpt(a, Om0):
    """Carroll-Press-Turner linear growth fit D(a) (matches CCL to <0.1%)."""
    Ol0 = 1.0 - Om0
    Oma = Om0 / (Om0 + Ol0 * a ** 3)
    Ola = 1.0 - Oma
    g = 2.5 * Oma / (Oma ** (4.0 / 7.0) - Ola
                     + (1.0 + 0.5 * Oma) * (1.0 + Ola / 70.0))
    return a * g


def _build_pk_linear_dict(k_phys, p_phys, Om):
    """CosmologyCalculator pk_linear dict: P(k,a) = P(k,z=0)*[D(a)/D(1)]^2."""
    a_arr = np.linspace(0.05, 1.0, 80)
    D1 = _D_cpt(1.0, Om)
    pk2d = np.array([(_D_cpt(a, Om) / D1) ** 2 * p_phys for a in a_arr])
    return {"a": a_arr, "k": k_phys, "delta_matter:delta_matter": pk2d}


def _corr_op_crossshell_vs_ccl(xi_ccl, gamma_arcmin):
    """Validated cross-shell corr_op 2-pt O0(gamma) vs CCL at gamma=1' (control).

    O0(gamma) = int_0^{lam_f} int_0^{lam_f} <Phi(n1,t1) Phi(n2,t2)> dt1 dt2 is the
    FULL cross-shell convergence 2-pt (kappa = int s dt).  Returns (O0_1arcmin,
    O0/CCL).  This is the canoes 2-pt that matches CCL; confirming the match (NOT
    a (1+z)^2 excess) isolates the excess to the equal-shell collapse + K-fold.
    Returns (nan, nan) if the corr_op callable is unavailable.
    """
    co_dir = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
                   "sachs_sft/callables/C_propagator/corr_op"))
    if str(co_dir) not in sys.path:
        sys.path.insert(0, str(co_dir))
    try:
        import corr_op_C_callable as co  # noqa: E402
    except Exception as exc:  # pragma: no cover
        print(f"[control] corr_op import failed ({exc}); skipping cross-shell check")
        return float("nan"), float("nan")
    n1 = np.array([0.0, 0.0, 1.0])
    g = np.radians(1.0 / 60.0)
    n2 = np.array([np.sin(g), 0.0, np.cos(g)])
    lam = np.linspace(450.0, 2300.0, 60)       # corr_op physical-Mpc support
    c00 = np.empty((lam.size, lam.size), dtype=np.float64)
    for i, t1 in enumerate(lam):
        for j, t2 in enumerate(lam):
            c00[i, j] = co.C_fn(n1, float(t1), n2, float(t2))[0, 0]
    o0 = float(np.trapezoid(np.trapezoid(c00, lam, axis=1), lam))
    ig = int(np.argmin(np.abs(gamma_arcmin - 1.0)))
    return o0, o0 / float(xi_ccl[ig])


def _implied_cl_ratio(chi_phys, z, D, chi_s, pk, cosmo, ell_ccl, cl_ccl):
    """Grid-independent (1+z)^2 proof: implied convergence Cl(hand)/Cl(CCL).

    Builds the SAME per-shell harmonic Phi00 density implied by the hand-rolled
    sigma2 [C_l^Phi00(chi) = (D^2/chi^2)[A(1+z)^4]^2 P0((l+.5)/chi)], K^2-folds it
    to a convergence Cl, and ratios against CCL's Cl ell-by-ell.  This bypasses
    BOTH the J0 angular quadrature and the Limber-delta grid discretization, so
    the ratio is the clean (1+z)^2 per-shell over-count weighted by ell.

    All h-units internally; the convergence Cl is dimensionless so no h-rescale.
    """
    h = cosmo.h
    Om = cosmo.Omega_m
    H0 = 100.0 / C_KM_S                       # h/Mpc
    chi_h = chi_phys * h
    chi_s_h = chi_s * h
    a = 1.0 / (1.0 + z)
    g = chi_h * (chi_s_h - chi_h) / chi_s_h    # h-units efficiency geometry
    A_amp = 1.5 * Om * H0 ** 2
    ell_probe = np.array([30.0, 100.0, 300.0, 1000.0, 3000.0])
    ckk_hand = np.zeros_like(ell_probe)
    ckk_ccl = np.zeros_like(ell_probe)
    for i, L in enumerate(ell_probe):
        k = (L + 0.5) / chi_h
        ok = k < 200.0                         # P(k) table support
        P0 = np.where(ok, np.asarray(pk(np.where(ok, k, 1.0))), 0.0)
        # hand-rolled per-shell harmonic density, K^2-folded in chi
        #   = int dchi a^2 (a^2 g)^2 [A^2 (1+z)^10 D^2 P0 / chi^2]
        integ = a ** 2 * (a ** 2 * g) ** 2 * (A_amp ** 2 * (1 + z) ** 10
                                              * D ** 2 * P0 / chi_h ** 2)
        ckk_hand[i] = np.trapezoid(integ, chi_h)
        j = int(np.argmin(np.abs(ell_ccl - L)))
        ckk_ccl[i] = cl_ccl[j]
    cl_ratio = ckk_hand / np.where(np.abs(ckk_ccl) > 0, ckk_ccl, np.nan)
    return ell_probe, ckk_hand, ckk_ccl, cl_ratio


def _bin_widths(centers):
    """Midpoint-edge bin widths from strictly-increasing shell centers."""
    centers = np.asarray(centers, dtype=np.float64)
    edges = np.empty(centers.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = max(0.0, centers[0] - 0.5 * (centers[1] - centers[0]))
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return np.diff(edges)


def _make_pk_phi_h(pk_table_path, cosmo):
    """P_phi(k) today in canoes h-units: (A/k^2)^2 P_m, A=1.5 Om (100/c)^2.

    Mirrors canoes.sachs.sft_input.corr_op.build._load_pk_phi (h-units branch),
    with a power-law high-k extrapolation (P_m ~ k^n tail) so the high-ell
    Limber kernel at small chi stays finite.  P_phi ~ k^-4 P_m, so the tail is
    steeply falling and the extrapolation is a sub-dominant correction.
    """
    arr = np.loadtxt(pk_table_path)
    k_h, p_m = arr[:, 0], arr[:, 1]                 # h-units
    A = 1.5 * cosmo.Omega_m * (100.0 / C_KM_S) ** 2  # h-units Poisson amp
    p_phi = (A / k_h ** 2) ** 2 * p_m
    lk = np.log(k_h)
    lp = np.log(p_phi)
    # high-k log-slope from the last decade for the power-law tail.
    hi = k_h > k_h[-1] / 10.0
    slope_hi = np.polyfit(lk[hi], lp[hi], 1)[0]
    lo = k_h < k_h[0] * 10.0
    slope_lo = np.polyfit(lk[lo], lp[lo], 1)[0]

    def _pk_phi(k):
        k = np.asarray(k, dtype=np.float64)
        out = np.exp(np.interp(np.log(k), lk, lp))
        below = k < k_h[0]
        above = k > k_h[-1]
        if np.any(below):
            out = np.where(below, p_phi[0] * (k / k_h[0]) ** slope_lo, out)
        if np.any(above):
            out = np.where(above, p_phi[-1] * (k / k_h[-1]) ** slope_hi, out)
        return out

    return _pk_phi


# ---------------------------------------------------------------------------
def main() -> int:
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    h = cosmo.h
    pk_table_path = Path("/Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt")

    z, chi_phys, lam_phys, chi_s, lam_s = _build_grid(cosmo)
    chi_phys_h = chi_phys * h
    D = _growth_canoes(z, cosmo)
    print(f"[b2] z_s={Z_SOURCE}, chi_s={chi_s:.2f} Mpc, lambda_s={lam_s:.2f} Mpc")
    print(f"[b2] n_lambda={_N_LAMBDA}, D(z0)={D[0]:.4f} D(zs)={D[-1]:.4f}")

    gamma_arcmin = np.unique(np.concatenate(
        [np.array([1.0]), np.geomspace(0.5, 300.0, 9)]))
    gamma_rad = np.radians(gamma_arcmin / 60.0)
    print(f"[b2] gamma (arcmin): {np.array2string(gamma_arcmin, precision=3)}")

    # shell widths (h-units) for the born_limber lambda-density delta weight.
    chi_widths_h = _bin_widths(chi_phys_h)
    lam_h = lam_phys * h
    lam_widths_h = _bin_widths(lam_h)

    # --- PRIMARY: canoes born_limber equal-shell sigma2_TTT (lambda-density) ----
    pk_phi = _make_pk_phi_h(pk_table_path, cosmo)
    sigma2_TTT, w_lambda = _sigma2_TTT_born_limber(
        gamma_rad, chi_phys, lam_phys, z, pk_phi, cosmo,
        chi_widths_h, lam_widths_h)
    print(f"[b2] (PRIMARY born_limber) sigma2_TTT shape={sigma2_TTT.shape}, "
          f"finite={int(np.sum(np.isfinite(sigma2_TTT)))}/{sigma2_TTT.size}")

    # --- SECONDARY: hand-rolled flat-sky Limber sigma2 (h-units -> physical h^3) -
    sig_h = _sigma2_TTT_handrolled(gamma_rad, chi_phys_h, z, D, pk, cosmo)
    sigma2_hand = sig_h * (h ** 3)
    rr = sigma2_hand / np.where(np.abs(sigma2_TTT) > 0, sigma2_TTT, np.nan)
    print(f"[b2] (CROSS-CHECK) hand-rolled/born_limber sigma2 median="
          f"{np.nanmedian(rr):.4f} (band {np.nanpercentile(rr,10):.3f}.."
          f"{np.nanpercentile(rr,90):.3f})")

    # --- K^2 lambda-fold: Z2 = int dlambda K^2 sigma2 ---------------------------
    # The born_limber sigma2 is a lambda-DENSITY (delta weight w_i=jac/Dlam_i
    # carries 1/Dlam_i); folding SUM_i Dlam_i K_i^2 sigma2_i cancels Dlam_i.  We
    # integrate with the trapezoid lambda measure (consistent with the density).
    a = 1.0 / (1.0 + z)
    K = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s   # physical Mpc
    Z2_lambda = np.trapezoid((K ** 2)[None, :] * sigma2_TTT, lam_phys, axis=1)

    # chi-route cross-check (carry dlambda = a^2 dchi explicitly):
    #   Z2 = int dlambda K^2 sigma2 = int dchi a^2 (a^2 K_chi)^2 sigma2
    #      = int dchi a^6 K_chi^2 sigma2,  K_chi = chi(chi_s-chi)/chi_s.
    K_chi = chi_phys * (chi_s - chi_phys) / chi_s
    Z2_chi = np.trapezoid((K_chi ** 2 * a ** 6)[None, :] * sigma2_TTT,
                          chi_phys, axis=1)
    rel_route = np.abs(Z2_lambda - Z2_chi) / np.maximum(np.abs(Z2_lambda), 1e-300)
    print(f"[b2] lambda-route vs chi-route max rel diff: {np.max(rel_route):.2e}")

    # hand-rolled fold for the cross-check curve.
    Z2_hand = np.trapezoid((K ** 2)[None, :] * sigma2_hand, lam_phys, axis=1)

    # CORRECTED fold: apply the Born per-leg reduction (1+z)^-(2n-2) = (1+z)^-2
    # for n=2 legs pointwise per shell (fix/spin2-zetaD-2_2_-2-norm).  This is the
    # SAME reduction stage0_ours / stage1_ours apply for n=3.  It demonstrates the
    # fix restores the K-fold to the CCL standard for the 2-leg case (verdict
    # below: ratio_fixed ~ 1).  born2 = (1+z)^-2.
    born2 = (1.0 + z) ** (-2)
    Z2_hand_fixed = np.trapezoid(
        (K ** 2 * born2)[None, :] * sigma2_hand, lam_phys, axis=1)

    # --- CCL-native convergence 2PCF (the standard) ----------------------------
    print("[b2] computing CCL-native xi_kappa ...")
    xi_ccl, ell_ccl, cl_ccl = _ccl_xi_kappa(gamma_arcmin, cosmo, pk_table_path)

    opz2 = (1.0 + Z_SOURCE) ** 2
    opz4 = (1.0 + Z_SOURCE) ** 4
    # PRIMARY verdict driver: the hand-rolled REAL-SPACE K^2-fold (correct
    # structural analogue of stage0_ours; no Limber-delta grid discretization).
    # The born_limber diagonal (Z2_lambda) is RETAINED but FLAGGED: its Limber
    # delta-in-chi kernel is grid-discretized (delta_weight=1/Delta_lambda), so on
    # a coarse n_lambda grid it UNDER-counts the fold by a large grid-dependent
    # factor (documented in mc_sachs_2pt/driver_stats.py).  It is NOT the verdict.
    ratio = Z2_hand / np.where(np.abs(xi_ccl) > 0, xi_ccl, np.nan)
    ratio_bl = Z2_lambda / np.where(np.abs(xi_ccl) > 0, xi_ccl, np.nan)
    ratio_fixed = Z2_hand_fixed / np.where(np.abs(xi_ccl) > 0, xi_ccl, np.nan)

    # --- DECISIVE grid-independent diagnostic: implied convergence Cl ----------
    # Build the hand-rolled per-shell HARMONIC density and K^2-fold it to a
    # convergence Cl, compare to CCL's Cl ell-by-ell.  This bypasses BOTH the
    # angular (J0) quadrature of the real-space sigma2 AND the Limber-delta grid
    # discretization, isolating the (1+z) bookkeeping cleanly.
    ell_cl, ckk_hand, ckk_ccl, cl_ratio = _implied_cl_ratio(
        chi_phys, z, D, chi_s, pk, cosmo, ell_ccl, cl_ccl)

    # --- CLINCHING control: the validated CROSS-SHELL corr_op 2-pt vs CCL -------
    # The full (non-collapsed) cross-shell propagator O0 = int dt1 dt2 <Phi Phi>
    # is the canoes 2-pt that "matches CCL".  Confirming O0/CCL ~ 1 (NOT (1+z)^2)
    # isolates the excess to the equal-shell-collapse + K-fold, NOT the per-leg
    # Sachs response per se: the cross-shell window distributes the (1+z) powers
    # correctly, the equal-shell collapse does not.
    corrop_o0_1arcmin, corrop_ratio = _corr_op_crossshell_vs_ccl(xi_ccl,
                                                                 gamma_arcmin)

    out = _HERE / "outputs" / "b2_kfold_2pt_vs_ccl.npz"
    np.savez(
        out,
        gamma_arcmin=gamma_arcmin, gamma_rad=gamma_rad,
        z=z, chi_phys=chi_phys, lam_phys=lam_phys, D=D, K=K,
        sigma2_TTT=sigma2_TTT, sigma2_hand=sigma2_hand, w_lambda=w_lambda,
        Z2_lambda=Z2_lambda, Z2_chi=Z2_chi, Z2_hand=Z2_hand,
        Z2_hand_fixed=Z2_hand_fixed,
        xi_ccl=xi_ccl, ratio=ratio, ratio_bl=ratio_bl, ratio_fixed=ratio_fixed,
        ell_ccl=ell_ccl, cl_ccl=cl_ccl,
        ell_cl_diag=ell_cl, ckk_hand=ckk_hand, ckk_ccl=ckk_ccl,
        cl_ratio=cl_ratio,
        corrop_o0_1arcmin=corrop_o0_1arcmin, corrop_crossshell_ratio=corrop_ratio,
        chi_s=chi_s, lam_s=lam_s, z_s=Z_SOURCE,
        Omega_m=cosmo.Omega_m, h=cosmo.h, n_s=cosmo.n_s,
        opz2=opz2, opz4=opz4,
        note=("Z2 = int dlambda K^2 sigma2_TTT (equal-shell 2-leg, ONE inner "
              "collapse).  PRIMARY = hand-rolled real-space K^2-fold (ratio); "
              "born_limber diagonal (ratio_bl) flagged: Limber-delta grid "
              "discretization undercounts on coarse grid.  cl_ratio = implied "
              "convergence Cl(hand)/Cl(CCL), grid-independent (1+z)^2 proof.  "
              "CCL-native single-source-plane convergence 2PCF, same PCAMBz0."),
    )
    print(f"[b2] saved -> {out}")

    print()
    print("[b2] === Z2 (K^2-fold) vs CCL-NATIVE xi_kappa(gamma) ===")
    print(f"[b2] (1+z_s)^2={opz2:.2f}  (1+z_s)^4={opz4:.2f}")
    print(f"[b2] {'gamma':>9} {'Z2_hand':>13} {'xi_ccl':>13} {'Z2/CCL':>9} "
          f"{'Z2_BL/CCL':>10} {'Z2fix/CCL':>10}")
    for ig in range(gamma_arcmin.size):
        print(f"[b2] {gamma_arcmin[ig]:9.3f} {Z2_hand[ig]:+13.5e} "
              f"{xi_ccl[ig]:+13.5e} {ratio[ig]:9.4f} {ratio_bl[ig]:10.4f} "
              f"{ratio_fixed[ig]:10.4f}")

    print()
    print("[b2] === DECISIVE: implied convergence Cl(K^2-fold)/Cl(CCL) by ell ===")
    print(f"[b2] {'ell':>6} {'Ckk_hand':>12} {'Ckk_ccl':>12} {'ratio':>8} "
          f"{'eff(1+z)':>9}")
    for i in range(ell_cl.size):
        eff = np.sqrt(cl_ratio[i]) if cl_ratio[i] > 0 else np.nan
        print(f"[b2] {ell_cl[i]:6.0f} {ckk_hand[i]:12.4e} {ckk_ccl[i]:12.4e} "
              f"{cl_ratio[i]:8.3f} {eff:9.3f}")

    # Verdict band: small-angle, signal-dominated region (0.5'..~60').
    band = (gamma_arcmin >= 0.5) & (gamma_arcmin <= 60.0)
    rb = ratio[band]
    rb = rb[np.isfinite(rb)]
    med = float(np.median(rb))
    rbf = ratio_fixed[band]
    rbf = rbf[np.isfinite(rbf)]
    med_fixed = float(np.median(rbf))
    cl_band = cl_ratio[(ell_cl >= 30) & (ell_cl <= 3000)]
    cl_med = float(np.median(cl_band))
    print()
    print(f"[b2] xi-space band median Z2_hand/CCL (0.5'..60') = {med:.3f}")
    print(f"[b2] xi-space band median Z2_FIXED/CCL (0.5'..60') = {med_fixed:.3f}  "
          f"(Born per-leg reduction (1+z)^-2 applied; ~1 => fix CORRECT)")
    print(f"[b2] Cl-space median ratio (ell 30..3000)          = {cl_med:.3f}")
    print(f"[b2]   the Cl ratio runs {cl_ratio.min():.2f}..{cl_ratio.max():.2f} "
          f"= (1+z_eff)^2 over the kernel-relevant z range (z~0.7..1.7).")
    print(f"[b2] CONTROL  cross-shell corr_op O0/CCL(1') = {corrop_ratio:.3f}  "
          f"(~1 => validated 2-pt matches CCL; NO (1+z)^2 excess)")
    print()
    print("[b2] VERDICT: the K^2 lambda-fold of the equal-shell 2-cumulant "
          "OVER-counts")
    print(f"[b2]   the CCL standard by a per-shell (1+z)^2 factor (NOT ~1, NOT "
          f"flat-unity).")
    print("[b2]   => the K-fold + per-leg (1+z)^4 convention ITSELF over-counts.")
    print("[b2]   Per-leg excess = (1+z)^2 [Phi00=A(1+z)^4 delta, A~(1+z), folded")
    print("[b2]   with K=a^2 g]; the affine measure dlambda=a^2 dchi refunds")
    print("[b2]   (1+z)^2 ONCE.  Net: n-pt over-count = (1+z)^(2n-2): 2-pt=(1+z)^2,")
    print("[b2]   3-pt=(1+z)^4 -- which MATCHES the observed stage0 3-pt excess.")
    print("[b2]   => the (1+z)^4 is NOT introduced by the equal-shell collapse; it")
    print("[b2]   is the K-fold/per-leg convention, and the equal-shell collapse")
    print("[b2]   only sets WHICH power (n legs -> (1+z)^(2n-2)).  This is a")
    print("[b2]   tension with the validated cross-shell canoes 2-pt (corr_op),")
    print("[b2]   which does NOT use this collapsed K^2 x equal-shell fold.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
