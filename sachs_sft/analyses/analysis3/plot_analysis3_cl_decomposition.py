"""Angular power spectra C_ell (the ell-space analogue of Fig. 11 / the 2PCF
NLO decomposition): the four lensing 2-point observables decomposed into
Order-0 / FF / FK / Full, obtained by a forward CURVED-SKY (full-sky) transform
of the SAME 2PCF data that Fig. 11 (analysis3_NLO_FFFK) shows.

The exact full-sky relation between a 2PCF and its angular power spectrum is

  C_ell = 2 pi integral_{-1}^{1} d(cos theta) xi(theta) d^ell_{m n}(theta),

with the reduced Wigner d-functions (Ng & Liu 1999; Kamionkowski+1998):
  kappa-kappa  -> d^ell_{0, 0} = P_ell   (C_ell^{kappa kappa})
  xi_+         -> d^ell_{2, 2}           (C_ell^{EE} + C_ell^{BB})
  xi_-         -> d^ell_{2,-2}           (C_ell^{EE} - C_ell^{BB})
  kappa-gamma  -> d^ell_{2, 0}           (C_ell^{kappa E})

The flat-sky Hankel kernels J_0 / J_4 / J_2 used previously are only the
small-angle / high-ell limits  d^ell_{m n}(theta) -> J_{m-n}(ell theta).  They
are INACCURATE at low ell / large angles -- exactly where the FK enhancement
sits -- so the curved-sky kernel is used here.  (The flat-sky transform is
retained as ``forward_hankel`` for the fullsky-vs-flatsky probe.)

The 2PCF is sampled on [0.5', 83 deg], so interpolation and finite angular
coverage affect the recovered spectra. Only the scalar channel uses endpoint
subtraction to estimate its constant large-separation contribution. This is not
an exact subtraction of the disconnected moment. Spin channels retain their
endpoint values because a constant is not generally orthogonal to their
Wigner-d kernels.

The built-in Wigner-d self-test checks the kernel norms. A separate pure-E
round trip in r1_sft061/pure_e_transform_check.py measures transform residuals
for the accepted Order-0 spectrum on this angular grid. Those residuals are
specific to that input and are not a universal error bound for FK spectra.

Run with the PyCCL interpreter (needs pyccl + scipy + matplotlib):
  <PyCCL python> plot_analysis3_cl_decomposition.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import eval_legendre, jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _plot_style import (  # noqa: E402
    PALETTE, annotate_sign_legend, apply_rcparams, plot_signed_line,
    plot_signed_markers,
)

RUNS = HERE.parent.parent / "sftwick_outputs" / "2PCF"
O0_NPZ = RUNS / "C_corr_op_O0" / "xi_C_corr_op_O0.npz"
FF_NPZ = RUNS / "C_corr_op_K_limber_FF" / "xi_C_corr_op_K_limber_FF.npz"
# The manuscript's FK: ell_max = 15360 with the permutation-aware vertex (2026-08-26).
# The June cut1000 sweep in C_corr_op_K_limber_FK/ is superseded and about 4.5x low at 0.5'.
FK_NPZ = RUNS / "C_corr_op_K_limber_FK_cut15360_permfix" / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz"
OUT_STEM = HERE / "outputs" / "analysis3_cl_O0_FF_FK"

sys.path.insert(0, str(HERE.parents[2] / "reproduce"))
from product_paths import resolve_product  # noqa: E402

O0_NPZ = resolve_product("order0", O0_NPZ)
FF_NPZ = resolve_product("ff", FF_NPZ)
FK_NPZ = resolve_product("fk", FK_NPZ)
LOAD_KW: dict[str, Any] = {"allow_pickle": True}  # trusted local sft-wick output

COLOR_O0, COLOR_FF, COLOR_FK = PALETTE[0], PALETTE[1], PALETTE[2]
MARKER_O0, MARKER_FF, MARKER_FK = "o", "^", "s"

# observable -> (title, combos, (m, n)); same component combos as the 2PCF
# figure, with the curved-sky Wigner-d kernel d^ell_{m n} (flat-sky limit J_{m-n}).
OBSERVABLE_SPECS = (
    ("kappa-kappa", r"$C_\ell^{\kappa\kappa}$  (from $\xi_\kappa$)",
     [((0, 0), +1.0)], (0, 0)),
    ("xi_plus", r"$C_\ell^{EE}+C_\ell^{BB}$  (from $\xi_+$)",
     [((1, 1), +1.0), ((2, 2), +1.0)], (2, 2)),
    ("xi_minus", r"$C_\ell^{EE}-C_\ell^{BB}$  (from $\xi_-$)",
     [((1, 1), +1.0), ((2, 2), -1.0)], (2, -2)),
    ("kappa_gamma_t", r"$C_\ell^{\kappa E}$  (from $\xi_{\kappa\gamma_t}$)",
     [((0, 1), -1.0)], (2, 0)),
)

Z_SOURCE = 5.0
# ell range [3, 1500]: the source 2PCF spans gamma in [0.5', 5000'].  The
# curved-sky transform must be evaluated at INTEGER ell (Legendre / Wigner-d);
# ~30 log-spaced multipoles, rounded to unique integers, match the 2PCF figure.
# Finite angular coverage and interpolation can bias the recovered spectra,
# including E/B separation. The live transform applies no angular apodisation.
ELL = np.array(sorted({int(round(v)) for v in np.geomspace(3.0, 1500.0, 30)}),
               dtype=float)


def _gamma_arcmin(xi, yi):
    xi = np.asarray(xi, float) / np.linalg.norm(xi)
    yi = np.asarray(yi, float) / np.linalg.norm(yi)
    return np.degrees(np.arccos(np.clip(np.dot(xi, yi), -1, 1))) * 60.0


def load_sweep_order(path: Path, order: int):
    with np.load(path, **LOAD_KW) as d:
        a, b, o = np.asarray(d["a"]), np.asarray(d["b"]), np.asarray(d["order"])
        v, x, y = np.asarray(d["value"], float), d["x"], d["y"]
    grouped: dict[tuple[int, int], np.ndarray] = {}
    g_ref = None
    for (pa, pb) in {(int(ai), int(bi)) for ai, bi in zip(a, b)}:
        m = (a == pa) & (b == pb) & (o == order)
        if not m.any():
            continue
        idx = np.flatnonzero(m)
        g = np.array([_gamma_arcmin(x[i], y[i]) for i in idx])
        s = np.argsort(g)
        if g_ref is None:
            g_ref = g[s]
        grouped[(pa, pb)] = v[m][s]
    if g_ref is None:
        raise ValueError(f"no rows at order={order} in {path}")
    return g_ref, grouped


def _combine(grouped, combos):
    out = None
    for (pa, pb), sign in combos:
        if (pa, pb) not in grouped:
            return None
        term = sign * grouped[(pa, pb)]
        out = term if out is None else out + term
    return out


# Apodisation (Tukey/Hann taper in log-theta) applied to the 2PCF before the
# Hankel.  The finite angular range (gamma <= 5000') is a HARD boundary; for the
# near-flat FK piece this rings in C_ell at low ell with period ~pi/theta_max.
# A raised-cosine taper of the large-angle end removes the boundary discontinuity
# and smooths C_ell at small ell.  The small-angle end is left untapered so the
# high-ell baseline (which matches PyCCL) is unaffected.
_APOD_LO = 0.0    # taper fraction at small theta (high-ell side); 0 = none
_APOD_HI = 0.20   # taper fraction at large theta (low-ell side)


def _tukey_log(n, a_lo, a_hi):
    """Raised-cosine (Tukey) window over a length-n log-theta grid."""
    w = np.ones(n)
    nlo, nhi = int(round(a_lo * n)), int(round(a_hi * n))
    if nlo > 1:
        i = np.arange(nlo)
        w[:nlo] = 0.5 * (1.0 - np.cos(np.pi * i / (nlo - 1)))
    if nhi > 1:
        i = np.arange(nhi)
        w[n - nhi:] = 0.5 * (1.0 - np.cos(np.pi * (nhi - 1 - i) / (nhi - 1)))
    return w


def forward_hankel(xi_theta, gamma_arcmin, ell, order, n_fine=6000, apodise=True):
    """FLAT-SKY (small-angle) reference transform, RETAINED FOR THE PROBE ONLY.

      C_ell = 2 pi int dtheta theta J_order(ell theta) xi(theta).

    This is the high-ell / small-angle limit of ``forward_curved`` and is
    inaccurate at low ell.  The figure uses ``forward_curved`` instead.
    """
    theta = gamma_arcmin * (math.pi / 180.0 / 60.0)
    lt = np.log(theta)
    lt_f = np.linspace(lt.min(), lt.max(), n_fine)
    theta_f = np.exp(lt_f)
    xi_f = PchipInterpolator(lt, xi_theta)(lt_f)
    if apodise:
        xi_f = xi_f * _tukey_log(n_fine, _APOD_LO, _APOD_HI)
    du = np.gradient(lt_f)
    J = jv(order, np.outer(np.asarray(ell, float), theta_f))  # (N_ell, n_fine)
    return 2.0 * np.pi * (J @ (theta_f**2 * xi_f * du))


# --- curved-sky (full-sky) transform via reduced Wigner d-functions ---------

def wigner_d(ell_int, cos_theta, m, n):
    """Reduced Wigner d^ell_{m n}(theta), vectorised over cos_theta, for the
    requested integer multipoles.  Returns {ell: array}.

    (m, n) = (0, 0) is Legendre P_ell.  The spin cases (2, 2), (2, -2), (2, 0)
    are built by the stable upward ell-recurrence

      ell sqrt((L+1)^2-m^2) sqrt((L+1)^2-n^2) d^{L+1}
        = (2L+1) (L(L+1) x - m n) d^L
          - (L+1) sqrt(L^2-m^2) sqrt(L^2-n^2) d^{L-1},   x = cos theta,

    seeded at ell_min = max(|m|, |n|) = 2 with the closed forms
      d^2_{2, 2}  = cos^4(theta/2),
      d^2_{2,-2}  = sin^4(theta/2),
      d^2_{2, 0}  = sqrt(3/8) sin^2(theta).
    """
    x = np.asarray(cos_theta, float)
    want = sorted({int(round(l)) for l in np.atleast_1d(ell_int)})
    out: dict[int, np.ndarray] = {}
    if (m, n) == (0, 0):
        return {L: eval_legendre(L, x) for L in want}
    s2sq = np.clip((1.0 - x) / 2.0, 0.0, None)  # sin^2(theta/2)
    c2sq = np.clip((1.0 + x) / 2.0, 0.0, None)  # cos^2(theta/2)
    if (m, n) == (2, 2):
        seed = c2sq * c2sq
    elif (m, n) == (2, -2):
        seed = s2sq * s2sq
    elif (m, n) == (2, 0):
        seed = math.sqrt(3.0 / 8.0) * (1.0 - x * x)
    else:
        raise ValueError(f"unsupported Wigner indices (m, n)=({m}, {n})")
    lmin, lmax = max(abs(m), abs(n)), max(want)
    dlm1 = np.zeros_like(x)
    dl = np.array(seed, dtype=float)
    if lmin in want:
        out[lmin] = dl.copy()
    for L in range(lmin, lmax):
        a1 = (2 * L + 1) * (L * (L + 1) * x - m * n)
        t = (L * L - m * m) * (L * L - n * n)
        a2 = (L + 1) * math.sqrt(t) if t > 0 else 0.0
        denom = L * math.sqrt(((L + 1) ** 2 - m * m) * ((L + 1) ** 2 - n * n))
        dlp1 = (a1 * dl - a2 * dlm1) / denom
        dlm1, dl = dl, dlp1
        if (L + 1) in want:
            out[L + 1] = dl.copy()
    return out


def _selftest_wigner(tol=1e-9):
    """Orthonormality int_{-1}^{1} (d^ell_{m n})^2 d(cos theta) = 2/(2 ell + 1)
    This checks the seed and recurrence at the listed multipoles.

    (d^ell_{m n})^2 is a polynomial of degree 2*ell in x = cos theta, so
    Gauss-Legendre quadrature integrates it EXACTLY (to machine precision)."""
    test_ell = [5, 17, 60, 200]
    nodes, weights = np.polynomial.legendre.leggauss(max(test_ell) + 3)
    for (m, n) in ((2, 2), (2, -2), (2, 0)):
        d = wigner_d(test_ell, nodes, m, n)
        for L, arr in d.items():
            norm = float(weights @ (arr * arr))
            exp = 2.0 / (2.0 * L + 1.0)
            assert abs(norm / exp - 1.0) < tol, (
                f"Wigner-d ortho fail (m,n)=({m},{n}) L={L}: "
                f"{norm:.10e} vs {exp:.10e}")
    print("[wigner-d self-test PASS: orthonormality 2/(2l+1) to machine prec.]")


def build_curved_matrix(gamma_arcmin, ell_int, m, n, n_fine=20000, apodise=False):
    """Pre-build the linear operator D so that C_ell = D @ xi_fine, i.e.

      C_ell = 2 pi int dtheta sin(theta) xi(theta) d^ell_{m n}(theta),

    on a fine log-theta PCHIP grid (no extrapolation beyond the data range).
    Returns (lt_data, lt_fine, apod, D); reused for every 2PCF on the same grid.

    ``apodise`` (default off) is retained for the fullsky-vs-flatsky probe.
    The live figure does not apply a taper. Scalar endpoint subtraction is
    selected separately in ``forward_curved``.
    """
    theta = gamma_arcmin * (math.pi / 180.0 / 60.0)
    lt = np.log(theta)
    lt_f = np.linspace(lt.min(), lt.max(), n_fine)
    theta_f = np.exp(lt_f)
    x = np.cos(theta_f)
    sth = np.sin(theta_f)
    dln = np.gradient(lt_f)
    apod = _tukey_log(n_fine, _APOD_LO, _APOD_HI) if apodise else np.ones(n_fine)
    dblock = wigner_d(ell_int, x, m, n)
    meas = 2.0 * np.pi * sth * theta_f * dln   # 2 pi sin(theta) dtheta
    D = np.array([meas * dblock[int(round(L))] for L in ell_int])
    return lt, lt_f, apod, D


def forward_curved(xi_theta, setup, dc_subtract=True):
    """Transform interpolated xi with the operator from ``build_curved_matrix``.

    The legacy ``dc_subtract=True`` default is preserved for historical callers
    and is appropriate only for the scalar (m, n) = (0, 0) endpoint convention.
    Live callers set it explicitly and disable it for all spin kernels.

    A known constant contributes only to the scalar monopole on the full
    sphere. Subtracting the final sampled value estimates that constant and
    also removes any connected correlation remaining at that angle. It is
    therefore not an exact subtraction of the disconnected moment. A constant
    can contribute to spin spectra, so the scalar argument does not apply.
    """
    lt, lt_f, _apod, D = setup
    xi = np.asarray(xi_theta, float)
    if dc_subtract:
        xi = xi - xi[-1]            # data is gamma-sorted; xi[-1] = xi(theta_max)
    xi_f = PchipInterpolator(lt, xi)(lt_f)
    return D @ xi_f


def _pyccl_ckk():
    """PyCCL Order-0 C_ell^{kappa kappa} for the fiducial cosmology (validation)."""
    try:
        import pyccl as ccl
    except Exception as exc:  # noqa: BLE001
        print(f"[pyccl unavailable: {exc}]")
        return None
    cosmo = ccl.Cosmology(Omega_c=0.31609 - 0.0492, Omega_b=0.0492, h=0.6711,
                          n_s=0.97, sigma8=0.81, matter_power_spectrum="linear")
    tr = ccl.CMBLensingTracer(cosmo, z_source=Z_SOURCE)
    return np.asarray(ccl.angular_cl(cosmo, tr, tr, ELL), float)


def main() -> int:
    apply_rcparams()
    import matplotlib.pyplot as plt

    _selftest_wigner()
    g0, o0g = load_sweep_order(O0_NPZ, 0)
    gff, ffg = load_sweep_order(FF_NPZ, 2)
    gfk, fkg = load_sweep_order(FK_NPZ, 2)
    print(f"gamma grid: {g0.min():.3f}..{g0.max():.1f} arcmin ({g0.size} pts)")

    ckk_ccl = _pyccl_ckk()

    panel_rc = {"axes.labelsize": 18, "axes.titlesize": 16,
                "xtick.labelsize": 14.5, "ytick.labelsize": 14.5,
                "legend.fontsize": 14.5}
    saved = {k: plt.rcParams[k] for k in panel_rc}
    plt.rcParams.update(panel_rc)

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(12.8, 12.0))
    outer = GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.20,
                     top=0.92, bottom=0.06, left=0.075, right=0.975)
    pref = ELL * (ELL + 1.0) / (2.0 * np.pi)   # band-power D_ell
    legend_handles = None
    print("\n=== Analysis 3 NLO C_ell breakdown (O0 / FF / FK), curved-sky ===")
    for i, (name, title, combos, (m, n)) in enumerate(OBSERVABLE_SPECS):
        r, c = i // 2, i % 2
        # each observable cell = tall main C_ell panel + short ratio panel
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.06)
        ax = fig.add_subplot(inner[0])
        axr = fig.add_subplot(inner[1], sharex=ax)
        o0 = _combine(o0g, combos); ff = _combine(ffg, combos); fk = _combine(fkg, combos)
        if o0 is None or ff is None or fk is None:
            ax.text(0.5, 0.5, "(component_pair missing)", transform=ax.transAxes,
                    ha="center", va="center"); ax.set_title(title); continue
        # the FF/FK sweeps share the O0 gamma grid; guard then sum the 2PCF.
        assert np.allclose(gff, g0) and np.allclose(gfk, g0), "gamma grids differ"

        # Cumulative C_ell: Order-0, Order-0+FF, and the TOTAL field (O0+FF+FK).
        # The quadrature matrix is shared, while each summed real-space
        # curve is interpolated separately with shape-preserving PCHIP.
        setup = build_curved_matrix(g0, ELL, m, n)
        dc_subtract = (m, n) == (0, 0)
        Cl_o0 = forward_curved(o0, setup, dc_subtract=dc_subtract)
        Cl_o0ff = forward_curved(o0 + ff, setup, dc_subtract=dc_subtract)
        Cl_full = forward_curved(o0 + ff + fk, setup, dc_subtract=dc_subtract)

        # Diagnostic against the conventional CCL convergence spectrum.
        # Spin prefactors, Ricci operators and finite-angle recovery differ,
        # so equality of every panel with this scalar reference is not required.
        if ckk_ccl is not None:
            lowm = ELL <= 60
            vr = np.nanmedian(np.abs(Cl_o0[lowm]) / np.abs(ckk_ccl[lowm]))
        else:
            vr = float("nan")
        print(f"{name} (d^l_{{{m},{n}}}): D_ell^O0(peak)={np.nanmax(np.abs(pref*Cl_o0)):.3e}  "
              f"D_ell^full(peak)={np.nanmax(np.abs(pref*Cl_full)):.3e}  "
              f"O0/PyCCL |ratio| median(l<=60)={vr:.3f}")

        # --- main panel: band-power C_ell, Full vs O0+FF vs Order-0 ---
        # marker/colour convention matched to Fig. 11 for ALL curves
        # (Full -> signed line; Order-0 -> COLOR_O0/"o"; O0+FF -> COLOR_FF/"^").
        # PyCCL is used for the console validation only, not drawn here.
        plot_signed_line(ax, ELL, pref * Cl_full, label="Full (O0+FF+FK)")
        # Higher zorder keeps O0+FF visible where it overlaps Order-0.
        plot_signed_markers(ax, ELL, pref * Cl_o0ff, color=COLOR_FF, marker=MARKER_FF, label="O0 + FF", zorder=4)
        plot_signed_markers(ax, ELL, pref * Cl_o0, color=COLOR_O0, marker=MARKER_O0, label="Order-0", zorder=3)
        ax.set_title(title)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylim(1e-7, 1e-3)
        ax.set_xlim(ELL.min(), ELL.max())
        plt.setp(ax.get_xticklabels(), visible=False)  # ell labels on the ratio panel
        if c == 0:
            ax.set_ylabel(r"$\ell(\ell+1)\,C_\ell/2\pi$")

        # --- ratio panel: fractional Order-2 correction (full - O0)/O0 ---
        # percent units so the strip reads directly as "size of the correction"
        corr = 100.0 * (Cl_full - Cl_o0) / Cl_o0
        # Finite-angle recovery and the large-separation FK tail limit the
        # interpretation at low ell. Draw the strip faint below ell ~ 50.
        plot_signed_line(axr, ELL, corr, color=COLOR_FK, lw=1.6, alpha=0.95,
                         sign_marker_size=4.0, sign_marker_alpha=0.75,
                         faint_outside=(50.0, None))
        axr.axhline(100.0, color="0.55", lw=0.8, ls=":")  # correction = Order-0
        axr.set_xscale("log"); axr.set_yscale("log")
        axr.set_xlim(ELL.min(), ELL.max())
        axr.set_ylim(0.1, 100.0)
        axr.set_xlabel(r"$\ell$")
        axr.tick_params(labelsize=12)
        if c == 0:
            axr.set_ylabel(r"$\frac{\mathrm{full}-\mathrm{O0}}{\mathrm{O0}}\;[\%]$", fontsize=15)

        if legend_handles is None:
            handles, labels = ax.get_legend_handles_labels()
            seen, uh, ul = set(), [], []
            for h, l in zip(handles, labels):
                if l not in seen:
                    seen.add(l); uh.append(h); ul.append(l)
            legend_handles = (uh, ul)
            annotate_sign_legend(ax)

    if legend_handles is not None:
        h, l = legend_handles
        fig.legend(h, l, loc="upper center", ncol=len(h),
                   bbox_to_anchor=(0.5, 0.985), frameon=False, fontsize=15.5)
    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf"):
        fig.savefig(OUT_STEM.with_suffix(ext))
    plt.close(fig)
    plt.rcParams.update(saved)
    print(f"\n-> {OUT_STEM.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
