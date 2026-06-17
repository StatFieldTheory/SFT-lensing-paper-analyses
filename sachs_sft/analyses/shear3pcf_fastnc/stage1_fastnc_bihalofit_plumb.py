#!/Users/zzhang/projects/fastnc_venv/bin/python
"""
stage1_fastnc_bihalofit_plumb.py
================================
Interpreter: /Users/zzhang/projects/fastnc_venv/bin/python  (ABSOLUTE; fastnc venv)

STAGE 1, STEP 1 (PLUMBING ONLY): order-of-magnitude geometry/units sanity at
LARGE gamma using the SHIPPED nonlinear BispectrumHalofit (BiHalofit). This is
NOT the physics validation (Step 2 / stage1_fastnc_tree.py is). At large gamma
(>= ~100') tree-level ~ nonlinear to ~10-20%, so the orchestrator can use these
BiHalofit numbers as a loose (~20%) plumbing check of the geometry, projection,
units, and overall magnitude. At arcmin gamma BiHalofit is 2-10x the tree value,
so do NOT compare these at small gamma.

The install was already de-risked by fastnc_smoke_bihalofit.py; this script just
records the large-gamma BiHalofit Gamma^0..3 on the pinned geometry for the
orchestrator, at the pinned plumbing point (gamma=150', phi=60deg) plus a couple
nearby large-gamma points.

Output: outputs/stage1_fastnc_bihalofit_plumb.npz
"""

import os
import json
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.integrate import quad
from astropy.cosmology import wCDM

np.seterr(all='ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'outputs')
os.makedirs(OUTDIR, exist_ok=True)

# Reuse the configured-cosmology helpers from the tree-level driver to keep the
# cosmology / P(k) / growth EXACTLY identical between Step 1 and Step 2.
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "stage1_tree", os.path.join(HERE, "stage1_fastnc_tree.py"))
_tree = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_tree)

GRID = _tree.GRID
Z_S = _tree.Z_S
H = _tree.H
OMEGA_M = _tree.OMEGA_M
OMEGA_DE = _tree.OMEGA_DE
W0 = _tree.W0
N_S = _tree.N_S
PROJECTION = _tree.PROJECTION
arcmin_to_rad = _tree.arcmin_to_rad

import fastnc
from fastnc.bispectrum import BispectrumHalofit
from fastnc.fastnc import FastNaturalComponents


def build_bihalofit(Lmax=20, Lmax_diag=20):
    k_h, pk_h = _tree.load_pklin()
    s8 = _tree.sigma8_from_pklin(k_h, pk_h)
    cosmo = wCDM(H0=100 * H, Om0=OMEGA_M, Ode0=OMEGA_DE, w0=W0,
                 meta={'n': N_S, 'sigma8': s8})
    z_arr = np.linspace(0.0, 5.0, 200)
    D_arr = _tree.build_growth(z_arr)

    print(f"  cosmology: Om0={OMEGA_M}, h={H}, w0={W0}, n_s={N_S}, sigma8={s8:.6f}")
    bs = BispectrumHalofit(Lmax=Lmax, Lmax_diag=Lmax_diag)
    bs.set_cosmology(cosmo, ns=N_S, sigma8=s8)
    bs.set_pklin(k_h, pk_h.copy())
    bs.set_lgr(z_arr, D_arr)
    bs.set_source_distribution([np.array([Z_S])], [np.array([1.0])])
    bs.compute_kernel()
    bs.interpolate()
    bs.decompose()
    return bs


def eval_one_gamma(bs, gamma_arcmin, phi_rad, Lmax=20, Mmax=20):
    """log-uniform 2-point stencil; diagonal element for the requested gamma."""
    t1 = arcmin_to_rad(np.array([gamma_arcmin, gamma_arcmin * 3.0]))
    fnc = FastNaturalComponents(Lmax=Lmax, Mmax=Mmax, projection=PROJECTION,
                                t1=t1, phi=phi_rad, verbose=False)
    fnc.set_bispectrum(bs)
    fnc.compute()
    return (fnc.Gamma0[0, 0, :].copy(), fnc.Gamma1[0, 0, :].copy(),
            fnc.Gamma2[0, 0, :].copy(), fnc.Gamma3[0, 0, :].copy())


def fmt_c(z):
    return f"({z.real:+.6e}) + ({z.imag:+.6e})j"


def main():
    print("=" * 78)
    print("STAGE 1, STEP 1 (PLUMBING): BiHalofit large-gamma Gamma^0..3")
    print(f"fastnc {fastnc.__version__}; projection={PROJECTION}; z_s={Z_S}")
    print("=" * 78)

    bs = build_bihalofit(Lmax=20, Lmax_diag=20)

    # Pinned plumbing point + nearby large-gamma points.
    gamma_list = np.array([100.0, 150.0, 200.0, 300.0])   # arcmin (all large)
    phi_list_deg = np.array([60.0, 90.0])                  # the plumbing phi + one more
    phi_rad = np.deg2rad(phi_list_deg)

    n_g, n_p = gamma_list.size, phi_list_deg.size
    G0 = np.empty((n_g, n_p), dtype=np.complex128)
    G1 = np.empty((n_g, n_p), dtype=np.complex128)
    G2 = np.empty((n_g, n_p), dtype=np.complex128)
    G3 = np.empty((n_g, n_p), dtype=np.complex128)

    for i, g in enumerate(gamma_list):
        g0, g1, g2, g3 = eval_one_gamma(bs, g, phi_rad)
        G0[i], G1[i], G2[i], G3[i] = g0, g1, g2, g3

    print("\n--- BiHalofit Gamma^0..3 (cent), large-gamma plumbing grid ---")
    for i, g in enumerate(gamma_list):
        print(f"\n  gamma={g:.0f}':")
        for j, p in enumerate(phi_list_deg):
            tag = "  <-- PINNED PLUMBING POINT" if (abs(g - 150.0) < 1e-6
                                                    and abs(p - 60.0) < 1e-6) else ""
            print(f"    phi={p:5.1f}deg: G0={G0[i,j].real:+.6e}  "
                  f"G1={G1[i,j].real:+.6e}  G2={fmt_c(G2[i,j])}{tag}")

    # Highlight the pinned plumbing point
    ig = int(np.argmin(np.abs(gamma_list - 150.0)))
    ip = int(np.argmin(np.abs(phi_list_deg - 60.0)))
    print("\n[PINNED PLUMBING POINT] gamma=150', phi=60deg (BiHalofit, nonlinear):")
    print(f"    Gamma0 = {fmt_c(G0[ig, ip])}")
    print(f"    Gamma1 = {fmt_c(G1[ig, ip])}")
    print(f"    Gamma2 = {fmt_c(G2[ig, ip])}")
    print(f"    Gamma3 = {fmt_c(G3[ig, ip])}")

    convention_note = (
        "STEP-1 PLUMBING ONLY (not the physics validation). fastnc shipped "
        "BispectrumHalofit (BiHalofit, NONLINEAR matter bispectrum) on the pinned "
        "SAS isoceles geometry t1=t2=gamma, opening phi, projection='cent' "
        "(centroid, arXiv:2309.08601). Same cosmology/P(k)/growth/z_s=5 as the "
        "tree-level Step-2 run. Large-gamma (>=100') order-of-magnitude check vs "
        "the our-side tree-level Gamma^i (loose ~20%; tree~nonlinear only at large "
        "gamma). At arcmin gamma BiHalofit is 2-10x the tree value -- do NOT compare "
        "there. Lmax=Lmax_diag=Mmax=20. fastnc commit d94fd2a (v"
        f"{fastnc.__version__})."
    )

    outpath = os.path.join(OUTDIR, 'stage1_fastnc_bihalofit_plumb.npz')
    np.savez(
        outpath,
        gamma_arcmin=gamma_list,
        phi_deg=phi_list_deg,
        Gamma0=G0, Gamma1=G1, Gamma2=G2, Gamma3=G3,
        projection=PROJECTION,
        convention_note=convention_note,
        plumbing_point_gamma_arcmin=150.0,
        plumbing_point_phi_deg=60.0,
        plumbing_point_gamma0=G0[ig, ip],
        plumbing_point_gamma1=G1[ig, ip],
        plumbing_point_gamma2=G2[ig, ip],
        plumbing_point_gamma3=G3[ig, ip],
        bispectrum_model='BiHalofit (nonlinear)',
        z_s=Z_S,
        fastnc_version=fastnc.__version__,
        fastnc_commit='d94fd2a',
    )
    print(f"\nSaved BiHalofit plumbing grid -> {outpath}")
    print("STAGE 1 STEP 1: DONE.")


if __name__ == '__main__':
    main()
