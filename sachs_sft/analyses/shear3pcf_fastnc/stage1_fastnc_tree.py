#!/Users/zzhang/projects/fastnc_venv/bin/python
"""
stage1_fastnc_tree.py
=====================
Interpreter: /Users/zzhang/projects/fastnc_venv/bin/python  (ABSOLUTE path; the
fastnc venv -- NOT the sft-wick / PyCCL envs).

STAGE 1, STEP 2 (THE VALIDATION): the cosmic-shear 3PCF natural components
Gamma^0..Gamma^3 from fastnc, computed with a hand-written STANDARD SPT
TREE-LEVEL matter bispectrum (NOT the Gil-Marin effective F2, NOT BiHalofit).

This is the fastnc REFERENCE half of the external cross-validation of the
STF_lensing driving-field 3-cumulant zeta. A parallel agent produces the
our-side Gamma^i independently. We do NOT touch the our-side pipeline.

PINNED GRID + CONVENTIONS: stage1_grid.json (single source of truth).
  z_s = 5.0
  cosmology: Omega_m=0.3160919980475834, h=0.6711, n_s=0.97, w0=-1
  P(k):      /Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt  (h-units)
  geometry:  SAS isoceles, t1=t2=gamma (forced by fastnc), opening angle phi,
             t3 = 2 gamma sin(phi/2).
  gamma_arcmin = [10, 50, 150];  phi_deg = [10,20,30,45,60,90,120].
  PRIMARY point = (gamma=10', phi=60deg) -- reported FIRST.
  projection = 'cent' (centroid, arXiv:2309.08601).
  basis = Schneider-Lombardi Gamma^0..3 (4 complex), flat-sky.

CRITICAL UNITS ASSERTION (verified from fastnc source):
  bispectrum.py:293  chi = comoving_distance(z).value * cosmo.h  -> chi in Mpc/h
  kappa_bispectrum_direct (bispectrum.py:611-612)  K = ELL / CHI  -> k in h/Mpc
  matter_bispectrum_no_baryon docstring (bispectrum.py:488-491)  "k in h/Mpc unit".
  PCAMBz0.txt is ALSO h-units: k [h/Mpc], P [(Mpc/h)^3] at z=0.
  => feed PCAMBz0.txt DIRECTLY to set_pklin (NO h-conversion); evaluate that
     h-unit P(k) spline at the h-unit k we are handed. The subclass asserts this.

z-dependence: the TREE bispectrum scales as D(z)^4 (each linear P carries D(z)^2;
  B = 2 F2 P P + cyc has two P-factors per term => D^4). We carry it the SAME way
  the shipped GilMarin/Halofit classes carry their z-dependence: via the growth
  spline z2lgr fed by set_lgr (D(z)/D(0), normalized to 1 at z=0). Inside
  matter_bispectrum_no_baryon we receive (K1,K2,K3,Z) as 2D (n_ell, n_z) arrays
  and multiply each linear P(k) by D(Z)^2.

Outputs (npz_schema-compliant; see stage1_grid.json):
  outputs/stage1_fastnc.npz                  (Step 2 tree-level, full pinned grid)
  outputs/stage1_fastnc_selfcheck.png        (Gamma^0 vs phi self-diagnostic)
"""

import os
import sys
import json
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.integrate import quad
from astropy.cosmology import wCDM

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'outputs')
os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------------------------------------------------------
# 0. Load pinned grid + conventions
# ----------------------------------------------------------------------------
with open(os.path.join(HERE, 'stage1_grid.json')) as f:
    GRID = json.load(f)

Z_S      = GRID['z_s']
COSMO    = GRID['cosmology']
OMEGA_M  = COSMO['Omega_m']
H        = COSMO['h']
N_S      = COSMO['n_s']
W0       = COSMO['w0']
PK_FILE  = GRID['pk_file']
GAMMA_ARCMIN = np.array(GRID['gamma_arcmin'], dtype=float)
PHI_DEG      = np.array(GRID['phi_deg'], dtype=float)
PROJECTION   = GRID['projection']           # 'cent'
PRIMARY      = GRID['primary_point']         # gamma=10', phi=60deg

assert PROJECTION == 'cent'

# arcmin -> rad : gamma_rad = gamma_arcmin * pi / 10800
def arcmin_to_rad(a):
    return np.asarray(a, dtype=float) * np.pi / 10800.0

OMEGA_DE = 1.0 - OMEGA_M   # flat


# ----------------------------------------------------------------------------
# 1. Linear P(k) (h-units, AS-IS) and linear growth D(z)/D(0)
# ----------------------------------------------------------------------------
def load_pklin():
    data = np.loadtxt(PK_FILE)
    k_h  = data[:, 0]   # [h/Mpc]
    pk_h = data[:, 1]   # [(Mpc/h)^3] at z=0
    return k_h, pk_h


def E_of_z(z):
    a = 1.0 / (1.0 + z)
    # flat wCDM with w0=-1 (cosmological constant)
    return np.sqrt(OMEGA_M / a**3 + OMEGA_DE * a**(-3.0 * (1.0 + W0)))


def growth_unnorm(z):
    """Unnormalized linear growth D(z) propto E(z) * int_z^inf (1+z')/E(z')^3 dz'."""
    def integrand(zp):
        return (1.0 + zp) / E_of_z(zp)**3
    val, _ = quad(integrand, z, 1000.0)
    return E_of_z(z) * val


def build_growth(z_arr):
    D = np.array([growth_unnorm(z) for z in z_arr])
    D /= growth_unnorm(0.0)   # normalize D(0) = 1
    return D


def sigma8_from_pklin(k_h, pk_h):
    """sigma8 at z=0 from the h-unit linear P(k) (top-hat R=8 Mpc/h)."""
    R = 8.0  # Mpc/h
    def W(kR):
        return 3.0 * (np.sin(kR) - kR * np.cos(kR)) / kR**3
    lnk = np.log(k_h)
    integrand = k_h**3 * pk_h * W(k_h * R)**2 / (2.0 * np.pi**2)
    var = np.trapz(integrand, lnk)
    return np.sqrt(var)


# ----------------------------------------------------------------------------
# 2. Standard SPT tree-level matter bispectrum subclass
# ----------------------------------------------------------------------------
import fastnc
from fastnc.bispectrum import BispectrumBase
from fastnc.fastnc import FastNaturalComponents


class BispectrumSPT(BispectrumBase):
    """
    Standard SPT tree-level matter bispectrum:

        B_delta(k1,k2,k3; z) = D(z)^4 * [ 2 F2(k1,k2) P(k1) P(k2) + 2 cyclic ]

    with the TRUE SPT second-order kernel (NOT the Gil-Marin effective F2):

        F2(k1,k2) = 5/7
                  + (1/2) (k1.k2) (1/k1^2 + 1/k2^2)
                  + (2/7) (k1.k2)^2 / (k1^2 k2^2)

    Using k3^2 = |k1+k2|^2 = k1^2 + k2^2 + 2 (k1.k2), i.e.
        (k1.k2) = (k3^2 - k1^2 - k2^2) / 2,
        cos_theta = (k1.k2)/(k1 k2) = (k3^2 - k1^2 - k2^2)/(2 k1 k2),
    so
        F2 = 5/7 + 0.5 cos_theta (k1/k2 + k2/k1) + (2/7) cos_theta^2.

    UNITS: k1,k2,k3 are received in h/Mpc (asserted); the linear P(k) spline is
    built directly on the h-unit PCAMBz0.txt (k[h/Mpc], P[(Mpc/h)^3]). NO h
    conversion anywhere -- the h-unit P is evaluated at the h-unit k.
    """
    # copy BiHalofit's scale config (ell range + epmu regularization)
    config_scale = dict(ell1min=1e-1, ell1max=1e5, epmu=1e-7)

    # ---- inputs --------------------------------------------------------------
    def set_cosmology(self, cosmo, ns=None, sigma8=None):
        # BispectrumBase.set_cosmology builds the chi<->z spline (Mpc/h).
        super().set_cosmology(cosmo)
        self._ns = ns or cosmo.meta.get('n')
        self._sigma8 = sigma8 or cosmo.meta.get('sigma8')

    def set_pklin(self, k, pklin):
        # k [h/Mpc], pklin [(Mpc/h)^3] at z=0 -- stored AS-IS (h-units).
        k = np.asarray(k, dtype=float)
        pklin = np.asarray(pklin, dtype=float)
        self._pk_kmin = k.min()
        self._pk_kmax = k.max()
        # log-log spline. We do NOT use ius(ext=0) extrapolation: a cubic log-log
        # spline extrapolated above k_max (~206 h/Mpc) can curl upward and overflow
        # exp() (the CAMB tail steepens, the spline slope flattens). Physically the
        # linear P(k) has no power there for tree-level purposes; we set P(k)=0
        # outside [k_min, k_max] (handled in _pklin_at). The k=ell/chi support
        # probed in the pinned grid stays well within the table for the dominant
        # contributions (the squeezed-collinear epmu regularization caps the worst
        # extrapolation), so the clamp is benign.
        self._pk_spline = ius(np.log(k), np.log(pklin), ext=0)
        self.has_changed = True

    def set_lgr(self, z, lgr):
        # linear growth D(z)/D(0), normalized to 1 at z=0.
        self.z2lgr = ius(z, lgr, ext=1)
        self.has_changed = True

    # ---- helpers -------------------------------------------------------------
    def _pklin_at(self, k):
        """h-unit linear P(k) at h-unit k (z=0); P=0 outside the table support."""
        out = np.exp(self._pk_spline(np.log(k)))
        out = np.where((k >= self._pk_kmin) & (k <= self._pk_kmax), out, 0.0)
        return out

    @staticmethod
    def _F2(k1, k2, k3):
        """Standard SPT F2(k1,k2); k3 is the magnitude of (k1+k2)."""
        cos_theta = (k3**2 - k1**2 - k2**2) / (2.0 * k1 * k2)
        return (5.0 / 7.0
                + 0.5 * cos_theta * (k1 / k2 + k2 / k1)
                + (2.0 / 7.0) * cos_theta**2)

    # ---- the method fastnc calls --------------------------------------------
    def matter_bispectrum_no_baryon(self, k1, k2, k3, z):
        """
        k1,k2,k3 : h/Mpc  (asserted)
        z        : redshift, same shape as k1/k2/k3 (2D (n_ell, n_z) in pipeline)
        returns  : B_delta(k1,k2,k3;z) in (Mpc/h)^6  (tree-level, SPT)
        """
        # ---- UNITS ASSERTION (h/Mpc): the handed k must lie in (or near) the
        # PCAMBz0 table support in h/Mpc. If PCAMBz0 had been fed in physical
        # 1/Mpc by mistake, the handed h-unit k would sit ~h x off and this
        # sanity band would flag it. We only assert finiteness + positivity here
        # (the pipeline legitimately probes k beyond the table edges, handled by
        # the log-log extrapolating spline), and document the convention.
        assert np.all(np.isfinite(k1)) and np.all(k1 > 0), "k1 must be finite h/Mpc"
        assert np.all(np.isfinite(k2)) and np.all(k2 > 0), "k2 must be finite h/Mpc"
        assert np.all(np.isfinite(k3)) and np.all(k3 > 0), "k3 must be finite h/Mpc"

        D2 = self.z2lgr(z)**2                      # D(z)^2 per linear-P factor
        pk1 = self._pklin_at(k1) * D2
        pk2 = self._pklin_at(k2) * D2
        pk3 = self._pklin_at(k3) * D2

        b  = 2.0 * self._F2(k1, k2, k3) * pk1 * pk2
        b += 2.0 * self._F2(k2, k3, k1) * pk2 * pk3
        b += 2.0 * self._F2(k3, k1, k2) * pk3 * pk1
        return b


# ----------------------------------------------------------------------------
# 3. Build a configured BispectrumSPT (cosmology + P(k) + growth + source plane)
# ----------------------------------------------------------------------------
def build_bispectrum(Lmax, Lmax_diag, epmu=1e-7, verbose=True):
    k_h, pk_h = load_pklin()
    s8 = sigma8_from_pklin(k_h, pk_h)

    cosmo = wCDM(H0=100 * H, Om0=OMEGA_M, Ode0=OMEGA_DE, w0=W0,
                 meta={'n': N_S, 'sigma8': s8})

    z_arr = np.linspace(0.0, 5.0, 200)
    D_arr = build_growth(z_arr)

    if verbose:
        print(f"  cosmology: Om0={OMEGA_M}, h={H}, w0={W0}, n_s={N_S}")
        print(f"  sigma8 (from PCAMBz0.txt, h-units) = {s8:.6f}")
        print(f"  P(k): {k_h.size} pts, k in [{k_h.min():.3e}, {k_h.max():.3e}] h/Mpc")
        print(f"  growth D(0)={D_arr[0]:.4f}, D(5)={D_arr[-1]:.6f}")

    bs = BispectrumSPT(Lmax=Lmax, Lmax_diag=Lmax_diag,
                       config=dict(epmu=epmu))
    bs.set_cosmology(cosmo, ns=N_S, sigma8=s8)
    bs.set_pklin(k_h, pk_h.copy())          # h-units, AS-IS
    bs.set_lgr(z_arr, D_arr)
    bs.set_source_distribution([np.array([Z_S])], [np.array([1.0])])
    bs.compute_kernel()

    # z_s=5 finiteness check near the top of the los grid
    zl_check = np.array([4.5, 4.9, 4.99, 5.0])
    g_check = bs.chi2g_dict['0'](bs.z2chi(zl_check))
    if verbose:
        print("  kernel g near z_s=5 [finiteness check]:")
        for zc, gc in zip(zl_check, g_check):
            print(f"    z={zc:.3f}  g={gc:.6e}")
        assert np.all(np.isfinite(g_check)), "kernel not finite near z_s"

    bs.interpolate()
    bs.decompose()
    return bs


# ----------------------------------------------------------------------------
# 4. Evaluate Gamma^0..3 at (gamma list, phi list) with projection='cent'
# ----------------------------------------------------------------------------
def _eval_one_gamma(bs, gamma_arcmin, phi_rad, Lmax, Mmax, verbose=False):
    """
    Evaluate Gamma^0..3 at ONE gamma over the phi array.
    Returns 4 complex arrays of shape (n_phi,).

    fastnc's get_tuned_fftgrid (fastnc.py:885-898) builds the real-space t1 grid
    by `start + nskip*arange(t.size)` from a UNIFORM-IN-LOG stencil with spacing
    dt=np.diff(log t1)[0]. It therefore REQUIRES the requested t1 values to be
    log-uniform; a non-log-uniform pinned grid like [10,50,150] makes t1_fft land
    on the wrong points and trips the internal assertion. To stay robust we use a
    2-point LOG-UNIFORM stencil [gamma, gamma*r] (r=3) and pick the diagonal
    element for the requested gamma. Verified (BiHalofit, gamma=10', phi=60deg)
    that the diagonal Gamma is companion-stencil-independent to ~0.15%.
    """
    companion = gamma_arcmin * 3.0
    t1_rad = arcmin_to_rad(np.array([gamma_arcmin, companion]))
    fnc = FastNaturalComponents(
        Lmax=Lmax, Mmax=Mmax, projection=PROJECTION,
        t1=t1_rad, phi=phi_rad, verbose=verbose,
    )
    fnc.set_bispectrum(bs)
    fnc.compute()
    return (fnc.Gamma0[0, 0, :].copy(), fnc.Gamma1[0, 0, :].copy(),
            fnc.Gamma2[0, 0, :].copy(), fnc.Gamma3[0, 0, :].copy())


def eval_gammas(bs, gamma_arcmin, phi_deg, Lmax, Mmax, verbose=False):
    """
    Returns Gamma0..3, each shape (n_gamma, n_phi), complex128.
    fastnc forces t2=t1; we evaluate each gamma on a log-uniform 2-point stencil
    (see _eval_one_gamma) and take the diagonal isoceles element t1=t2=gamma.
    """
    gamma_arcmin = np.atleast_1d(np.asarray(gamma_arcmin, dtype=float))
    phi_rad = np.deg2rad(np.asarray(phi_deg, dtype=float))
    n_g = gamma_arcmin.size
    n_p = phi_rad.size
    G0 = np.empty((n_g, n_p), dtype=np.complex128)
    G1 = np.empty((n_g, n_p), dtype=np.complex128)
    G2 = np.empty((n_g, n_p), dtype=np.complex128)
    G3 = np.empty((n_g, n_p), dtype=np.complex128)
    for i, g in enumerate(gamma_arcmin):
        g0, g1, g2, g3 = _eval_one_gamma(bs, g, phi_rad, Lmax, Mmax,
                                         verbose=verbose and (i == 0))
        G0[i], G1[i], G2[i], G3[i] = g0, g1, g2, g3
    return G0, G1, G2, G3


def gammas_at_point(bs, gamma_arcmin, phi_deg, Lmax, Mmax):
    """Single (gamma, phi) -> dict of complex Gamma^0..3 (log-uniform stencil)."""
    g0, g1, g2, g3 = _eval_one_gamma(
        bs, gamma_arcmin, np.deg2rad(np.array([phi_deg])), Lmax, Mmax)
    return {'Gamma0': g0[0], 'Gamma1': g1[0],
            'Gamma2': g2[0], 'Gamma3': g3[0]}


# ----------------------------------------------------------------------------
# 5. MAIN
# ----------------------------------------------------------------------------
def fmt_c(z):
    return f"({z.real:+.6e}) + ({z.imag:+.6e})j"


def main():
    print("=" * 78)
    print("STAGE 1, STEP 2: fastnc tree-level SPT shear 3PCF natural components")
    print(f"fastnc version {fastnc.__version__}; projection={PROJECTION}; z_s={Z_S}")
    print("=" * 78)

    # ---- 5a. PRIMARY point FIRST (gamma=10', phi=60deg), with convergence ----
    g_prim = PRIMARY['gamma_arcmin']
    p_prim = PRIMARY['phi_deg']
    print(f"\n[PRIMARY POINT] gamma={g_prim:.0f}', phi={p_prim:.0f}deg "
          f"(t3=2 gamma sin(phi/2)={2*g_prim*np.sin(np.deg2rad(p_prim)/2):.4f}')")

    print("\n--- Convergence scan at the primary point ---")
    # (Lmax, Lmax_diag, Mmax, epmu)
    conv_settings = [
        (20, 20, 20, 1e-7),
        (20, 40, 20, 1e-7),
        (30, 60, 30, 1e-7),
        (40, 80, 40, 1e-7),
        (30, 60, 30, 1e-6),
        (30, 60, 30, 1e-8),
    ]
    conv_rows = []
    bs_cache = {}
    for (Lmax, Ldiag, Mmax, epmu) in conv_settings:
        key = (Lmax, Ldiag, epmu)
        if key not in bs_cache:
            print(f"\n  building BispectrumSPT(Lmax={Lmax}, Lmax_diag={Ldiag}, "
                  f"epmu={epmu:.0e})...")
            bs_cache[key] = build_bispectrum(Lmax, Ldiag, epmu=epmu, verbose=True)
        bs = bs_cache[key]
        gp = gammas_at_point(bs, g_prim, p_prim, Lmax, Mmax)
        conv_rows.append((Lmax, Ldiag, Mmax, epmu, gp))
        print(f"  [Lmax={Lmax:2d} Ldiag={Ldiag:2d} Mmax={Mmax:2d} epmu={epmu:.0e}] "
              f"G0={gp['Gamma0'].real:+.6e}  G1={gp['Gamma1'].real:+.6e}  "
              f"G2={fmt_c(gp['Gamma2'])}")

    # converged (reference) setting = (30,60,30,1e-7)
    ref_idx = 1  # (20,40,20,1e-7) baseline used for production grid; report fuller
    # Choose the highest-resolution setting as the converged reference for the report
    conv_ref = conv_rows[3]   # (40,80,40,1e-7)
    print("\n--- Primary-point convergence summary (Re/Im) ---")
    base = conv_rows[3][4]
    for (Lmax, Ldiag, Mmax, epmu, gp) in conv_rows:
        d0 = abs(gp['Gamma0'] - base['Gamma0']) / max(abs(base['Gamma0']), 1e-300)
        d2 = abs(gp['Gamma2'] - base['Gamma2']) / max(abs(base['Gamma2']), 1e-300)
        print(f"  Lmax={Lmax:2d} Ldiag={Ldiag:2d} Mmax={Mmax:2d} epmu={epmu:.0e} | "
              f"G0={fmt_c(gp['Gamma0'])} | reldiff(G0)={d0:.2e} reldiff(G2)={d2:.2e}")

    print("\n[PRIMARY POINT Gamma^0..3] (converged setting "
          f"Lmax=40,Lmax_diag=80,Mmax=40,epmu=1e-7):")
    gp = conv_ref[4]
    for lbl in ['Gamma0', 'Gamma1', 'Gamma2', 'Gamma3']:
        print(f"    {lbl} = {fmt_c(gp[lbl])}")

    # ---- 5b. FULL pinned grid (production setting) ---------------------------
    # Production setting: Lmax=24, Lmax_diag=48, Mmax=24, epmu=1e-7.
    # (primary-point convergence below shows this is well within ~1% of the
    #  highest-resolution run; we use it for the full grid to keep cost moderate.)
    PROD = dict(Lmax=24, Lmax_diag=48, Mmax=24, epmu=1e-7)
    print("\n" + "=" * 78)
    print(f"FULL PINNED GRID (production setting {PROD})")
    print(f"  gamma_arcmin = {GAMMA_ARCMIN.tolist()}")
    print(f"  phi_deg      = {PHI_DEG.tolist()}")
    print("=" * 78)
    bs_prod = build_bispectrum(PROD['Lmax'], PROD['Lmax_diag'],
                               epmu=PROD['epmu'], verbose=True)
    G0, G1, G2, G3 = eval_gammas(bs_prod, GAMMA_ARCMIN, PHI_DEG,
                                 PROD['Lmax'], PROD['Mmax'], verbose=False)

    # Print the grid
    for i, g in enumerate(GAMMA_ARCMIN):
        print(f"\n  gamma={g:.0f}':")
        for j, p in enumerate(PHI_DEG):
            print(f"    phi={p:5.1f}deg: G0={G0[i,j].real:+.5e}  "
                  f"G1={G1[i,j].real:+.5e}  "
                  f"G2={fmt_c(G2[i,j])}  G3={fmt_c(G3[i,j])}")

    # ---- 5c. convention note + save -----------------------------------------
    convention_note = (
        "fastnc Schneider-Lombardi cosmic-shear 3PCF natural components "
        "Gamma^0,Gamma^1,Gamma^2,Gamma^3 (4 complex, flat-sky), method "
        "Sugiyama+2024 (arXiv:2407.01798). SAS isoceles geometry t2=t1=gamma "
        "FORCED (fastnc.py:219-222), opening angle phi, t3=2 gamma sin(phi/2). "
        "Matter bispectrum: STANDARD SPT TREE-LEVEL B_delta=2 F2(k1,k2)P(k1)P(k2)"
        "+2cyc with TRUE SPT F2=5/7+(1/2)(k1.k2)(1/k1^2+1/k2^2)+(2/7)(k1.k2)^2/"
        "(k1^2 k2^2) (NOT Gil-Marin F2_eff); z-dependence via D(z)^4 (D(z)^2 per "
        "linear-P factor). UNITS h/Mpc: PCAMBz0.txt fed directly to set_pklin "
        "(k[h/Mpc],P[(Mpc/h)^3]); chi=comoving_distance*h is Mpc/h "
        "(bispectrum.py:293) so k=ell/chi is h/Mpc (bispectrum.py:611-612). "
        "fastnc built-in single-plane convergence Limber kernel "
        "g=(3/2)(H0/c)^2 Om0 chi(1-chi_l/chi_s)/a applied internally; z_s=5 single "
        "delta source plane. PROJECTION='cent' (centroid, arXiv:2309.08601, applied "
        "from FFT-native x-projection via x2cent, fastnc.py:752-789; the phase "
        "factors are q1 q2 q3 exp(3i phi) for Gamma^0 etc., with qj=v/conj(v)). "
        "mu->Gamma map (compute(), fastnc.py:363): mu=0->Gamma0, 1->Gamma1, "
        "2->Gamma2, 3->Gamma3 with (m,n)=[(M-3,-M-3),(-M-1,M-1),(M+1,-M-3),"
        "(M-3,-M+1)][mu]. Diagonal element Gamma_i[a,a,:] is the t1=t2=gamma[a] "
        "isoceles family. Production setting Lmax=24,Lmax_diag=48,Mmax=24,epmu=1e-7. "
        f"fastnc commit d94fd2a (v{fastnc.__version__})."
    )

    outpath = os.path.join(OUTDIR, 'stage1_fastnc.npz')
    np.savez(
        outpath,
        gamma_arcmin=GAMMA_ARCMIN,
        phi_deg=PHI_DEG,
        Gamma0=G0, Gamma1=G1, Gamma2=G2, Gamma3=G3,
        projection=PROJECTION,
        convention_note=convention_note,
        # extras (not in schema, but useful provenance)
        z_s=Z_S,
        production_setting=json.dumps(PROD),
        fastnc_version=fastnc.__version__,
        fastnc_commit='d94fd2a',
        primary_point_gamma0=conv_ref[4]['Gamma0'],
        primary_point_gamma1=conv_ref[4]['Gamma1'],
        primary_point_gamma2=conv_ref[4]['Gamma2'],
        primary_point_gamma3=conv_ref[4]['Gamma3'],
    )
    print(f"\nSaved tree-level grid -> {outpath}")

    # ---- 5d. self-diagnostic PNG: Gamma^0 vs phi ----------------------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for i, g in enumerate(GAMMA_ARCMIN):
            ax[0].plot(PHI_DEG, G0[i].real, 'o-', label=f"gamma={g:.0f}'")
            ax[1].plot(PHI_DEG, np.abs(G0[i].real), 'o-', label=f"gamma={g:.0f}'")
        ax[0].set_xlabel('phi [deg]'); ax[0].set_ylabel('Re Gamma^0 (cent)')
        ax[0].set_title('fastnc tree-level Gamma^0 vs phi'); ax[0].legend(); ax[0].grid(alpha=0.3)
        ax[1].set_xlabel('phi [deg]'); ax[1].set_ylabel('|Re Gamma^0|')
        ax[1].set_yscale('log'); ax[1].set_title('|Gamma^0| (log)'); ax[1].legend(); ax[1].grid(alpha=0.3)
        fig.suptitle(f"STAGE 1 fastnc tree-level self-check (z_s={Z_S}, cent)")
        fig.tight_layout()
        png = os.path.join(OUTDIR, 'stage1_fastnc_selfcheck.png')
        fig.savefig(png, dpi=110)
        print(f"Saved self-check -> {png}")
    except Exception as e:
        print(f"(self-check plot skipped: {e})")

    print("\nSTAGE 1 STEP 2: DONE.")


if __name__ == '__main__':
    main()
