"""make_appendix_curve_data.py -- FINE-grid cosmic-shear 3PCF curve data.

DATA ONLY (no figure). Produces a tidy npz with BOTH sides of the appendix
validation comparison (our SFT driving-field 3-cumulant vs the external fastnc
tree-level code), evaluated on a DENSE phi grid so a smooth curve can be drawn.

GRID (fine; everything else identical to stage1_grid.json):
  phi_deg      = 24 points evenly spaced in [10, 120] deg  (validated range)
  gamma_arcmin = [10.0, 50.0]                              (150' skipped to save time)
  z_s = 5, projection = 'cent'

Two compute environments are required (the two stage1 scripts use different
Python venvs and these cannot coexist in one interpreter):
  * OUR side  -- sft-wick env:
        /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
    Replicates (byte-identical) the compute setup of stage1_ours_v2_JACFIX.main():
    same cosmology, pk_file, _N_LAMBDA, kernel, common_low/high/mod_low/mod_high
    dicts, born3=ones, the fold, and the EXACT G0..G3 slot+centroid-phase
    assembly. The module-level helpers (_build_lambda_grid, _convergence_kernel,
    _sas_cosine_triples, centroid_phase, _DMOD_SPINS, _N_LAMBDA) are IMPORTED
    from stage1_ours_v2_JACFIX so no physics is duplicated; only the main()
    closure body (which is not importable) is replicated verbatim.
  * fastnc side -- fastnc venv:
        /Users/zzhang/projects/fastnc_venv/bin/python
    Reuses stage1_fastnc_tree.build_bispectrum + eval_gammas (the SAME tree
    computation, production setting Lmax=24,Lmax_diag=48,Mmax=24,epmu=1e-7).

USAGE
  Orchestrator (run under sft-wick env -- DEFAULT):
        sft-wick-python make_appendix_curve_data.py
    -> computes OUR side in-process, subprocess-calls the fastnc venv for the
       fnc side, then combines into outputs/appendix_fastnc_curve.npz.
  Sub-modes (used internally, but runnable standalone):
        <any python> make_appendix_curve_data.py --side ours --out <npz>
        fastnc-python make_appendix_curve_data.py --side fnc  --out <npz>

OUTPUT outputs/appendix_fastnc_curve.npz keys:
  phi_deg (n_phi,), gamma_arcmin (2,),
  ours_G0,_G1,_G2,_G3 / fnc_G0,_G1,_G2,_G3  complex (2, n_phi),
  ours_fi / fnc_fi  real (2, n_phi)  = sqrt(|G0|^2+|G1|^2+|G2|^2+|G3|^2),
  note  (str).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_OUT_DIR = _HERE / "outputs"
_FINAL_NPZ = _OUT_DIR / "appendix_fastnc_curve.npz"

_SFTWICK_PY = "/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python"
_FASTNC_PY = "/Users/zzhang/projects/fastnc_venv/bin/python"

# ---------------------------------------------------------------------------
# FINE GRID (validated range only; do NOT extrapolate past [10, 120] deg).
# ---------------------------------------------------------------------------
PHI_DEG_FINE = np.linspace(10.0, 120.0, 24)          # 24 pts
GAMMA_ARCMIN_FINE = np.array([10.0, 50.0])           # skip 150' to save time
Z_S = 5.0
PROJECTION = "cent"


def _frame_invariant(G0, G1, G2, G3):
    """Projection-independent 3PCF amplitude fi = sqrt(sum_i |Gamma^i|^2)."""
    return np.sqrt(np.abs(G0) ** 2 + np.abs(G1) ** 2
                   + np.abs(G2) ** 2 + np.abs(G3) ** 2)


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ===========================================================================
# OUR SIDE (sft-wick env). Replicates stage1_ours_v2_JACFIX.main() setup
# verbatim, reusing its importable module-level helpers so the physics and
# convention are byte-identical (we only re-create the main() closure body,
# which is not importable).
# ===========================================================================
def compute_ours(phi_deg: np.ndarray, gamma_arcmin: np.ndarray):
    M = _load_module(_HERE / "stage1_ours_v2_JACFIX.py", "stage1_ours_v2_JACFIX")

    # --- import the byte-identical helpers / constants ---------------------
    from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
    _build_lambda_grid = M._build_lambda_grid
    _convergence_kernel = M._convergence_kernel
    _sas_cosine_triples = M._sas_cosine_triples
    centroid_phase = M.centroid_phase
    _DMOD_SPINS = M._DMOD_SPINS
    _N_LAMBDA = M._N_LAMBDA

    phi_rad = np.radians(np.asarray(phi_deg, dtype=np.float64))

    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)

    z, chi_phys, lam_phys, lam_h, chi_s, lam_s = _build_lambda_grid(
        cosmo, Z_S, _N_LAMBDA)
    K = _convergence_kernel(chi_phys, z, chi_s)
    print(f"[ours] z_s={Z_S}, chi_s={chi_s:.2f} Mpc, lam_s={lam_s:.2f} Mpc, "
          f"_N_LAMBDA={_N_LAMBDA}")

    from canoes.sachs import (  # noqa: E402
        compute_kappa3_sigma3_high,
        compute_kappa3_zeta_table,
        kappa3_combine_low_high,
    )
    from canoes.sachs.kappa3 import (  # noqa: E402
        compute_kappa3_mod_zeta_equal_time,
        compute_kappa3_mod_sigma3_high,
    )

    # --- common_low/high/mod_low/mod_high: byte-identical to JACFIX ---------
    common_low = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128, spt_kind="tree_phi",
        units="physical", lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda", channels=("PPP",), parallelism="sequential",
        shell_block_size=4, triple_chunk_size=4096, n_workers=1,
        fuse_tree_phi_channels=True,
    )
    common_high = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000, ell_min=1e-3,
        n_ell=96, n_phi=64, channels=("PPP",), units="physical",
        lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda",
    )
    mod_low = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128, spt_kind="tree_phi",
        units="physical", lambda_convention="project", radial_measure="lambda",
    )
    mod_high = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000, ell_min=1e-3,
        n_ell=96, n_phi=64, units="physical", lambda_convention="project",
        radial_measure="lambda",
    )

    def build_zeta_channels(triples):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            low = compute_kappa3_zeta_table(triples, lam_h, **common_low)
            high = compute_kappa3_sigma3_high(triples, lam_h, **common_high)
            comb = kappa3_combine_low_high(low, high)
            zPPP = np.asarray(comb.zeta_PPP, dtype=np.float64)
            zD = {}
            for mu, sw in _DMOD_SPINS.items():
                dlo = compute_kappa3_mod_zeta_equal_time(
                    triples, lam_h, _spin_weights_dmod=sw, **mod_low)[1]
                dhi = compute_kappa3_mod_sigma3_high(
                    triples, lam_h, _spin_weights_dmod=sw, **mod_high)[1]
                zD[mu] = np.asarray(dlo + dhi, dtype=np.float64)
        return zPPP, zD[1], zD[2], zD[3]

    # JACFIX: born3 = ones (the (1+z)^-4 enters ONCE, from the vertex).
    born3 = np.ones_like(z)

    def fold(zeta):
        return np.trapezoid((K ** 3 * born3)[None, :] * zeta, lam_phys, axis=1)

    def natural_components(gamma_arcmin_val, phi_rad_arr):
        gamma_rad = np.radians(gamma_arcmin_val / 60.0)
        triples = _sas_cosine_triples(gamma_rad, phi_rad_arr)
        zPPP, zD1, zD2, zD3 = build_zeta_channels(triples)
        fPPP, fD1, fD2, fD3 = fold(zPPP), fold(zD1), fold(zD2), fold(zD3)
        sign = -1.0
        # SLOT RELABEL (Schneider-Lombardi / fastnc SAS convention), identical
        # to stage1_ours_v2_JACFIX.main(): APEX(leg3)->Gamma^1 (real),
        # BASE legs 1,2 -> Gamma^2,Gamma^3 so Gamma^2=conj(Gamma^3).
        G0 = sign * fPPP * centroid_phase(0, phi_rad_arr)
        G1 = sign * fD3 * centroid_phase(3, phi_rad_arr)   # APEX -> slot 1
        G2 = sign * fD1 * centroid_phase(1, phi_rad_arr)   # BASE leg 1 -> slot 2
        G3 = sign * fD2 * centroid_phase(2, phi_rad_arr)   # BASE leg 2 -> slot 3
        return (np.asarray(G0, np.complex128), np.asarray(G1, np.complex128),
                np.asarray(G2, np.complex128), np.asarray(G3, np.complex128))

    n_g, n_p = gamma_arcmin.size, phi_deg.size
    G0 = np.zeros((n_g, n_p), np.complex128)
    G1 = np.zeros((n_g, n_p), np.complex128)
    G2 = np.zeros((n_g, n_p), np.complex128)
    G3 = np.zeros((n_g, n_p), np.complex128)
    for ig, gam in enumerate(gamma_arcmin):
        t0 = time.time()
        g0, g1, g2, g3 = natural_components(gam, phi_rad)
        G0[ig], G1[ig], G2[ig], G3[ig] = g0, g1, g2, g3
        print(f"[ours]  gamma={gam:7.2f}' done in {time.time() - t0:6.1f}s "
              f"({phi_deg.size} phi)")
    fi = _frame_invariant(G0, G1, G2, G3)
    return dict(G0=G0, G1=G1, G2=G2, G3=G3, fi=fi)


# ===========================================================================
# fastnc SIDE (fastnc venv). Reuses stage1_fastnc_tree's tree computation.
# ===========================================================================
def compute_fnc(phi_deg: np.ndarray, gamma_arcmin: np.ndarray):
    F = _load_module(_HERE / "stage1_fastnc_tree.py", "stage1_fastnc_tree")
    assert F.PROJECTION == PROJECTION
    assert float(F.Z_S) == Z_S
    # Production setting -- identical to stage1_fastnc_tree.main() full grid.
    PROD = dict(Lmax=24, Lmax_diag=48, Mmax=24, epmu=1e-7)
    print(f"[fnc] building BispectrumSPT {PROD} (fastnc {F.fastnc.__version__})")
    t0 = time.time()
    bs = F.build_bispectrum(PROD["Lmax"], PROD["Lmax_diag"],
                            epmu=PROD["epmu"], verbose=True)
    print(f"[fnc] bispectrum built in {time.time() - t0:.1f}s")
    t0 = time.time()
    G0, G1, G2, G3 = F.eval_gammas(bs, gamma_arcmin, phi_deg,
                                   PROD["Lmax"], PROD["Mmax"], verbose=False)
    print(f"[fnc] eval_gammas ({gamma_arcmin.size} gamma x {phi_deg.size} phi) "
          f"in {time.time() - t0:.1f}s")
    fi = _frame_invariant(G0, G1, G2, G3)
    return dict(G0=G0, G1=G1, G2=G2, G3=G3, fi=fi)


# ===========================================================================
# Sub-mode driver: compute ONE side, dump to a plain npz.
# ===========================================================================
def _run_side(side: str, out_path: Path):
    if side == "ours":
        res = compute_ours(PHI_DEG_FINE, GAMMA_ARCMIN_FINE)
    elif side == "fnc":
        res = compute_fnc(PHI_DEG_FINE, GAMMA_ARCMIN_FINE)
    else:
        raise ValueError(side)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_path,
        phi_deg=PHI_DEG_FINE, gamma_arcmin=GAMMA_ARCMIN_FINE,
        G0=res["G0"], G1=res["G1"], G2=res["G2"], G3=res["G3"], fi=res["fi"],
    )
    assert np.all(np.isfinite(res["fi"])), f"{side}: non-finite fi"
    print(f"[{side}] saved -> {out_path}")


# ===========================================================================
# Orchestrator (sft-wick env): ours in-process + fnc via subprocess, combine.
# ===========================================================================
def _orchestrate():
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    ours_tmp = _OUT_DIR / "_appendix_curve_ours_tmp.npz"
    fnc_tmp = _OUT_DIR / "_appendix_curve_fnc_tmp.npz"

    # ---- OUR side, in-process (we are the sft-wick interpreter) -----------
    print("=" * 78)
    print("ORCHESTRATOR: computing OUR side in-process (sft-wick env)")
    print("=" * 78)
    t_ours0 = time.time()
    _run_side("ours", ours_tmp)
    t_ours = time.time() - t_ours0
    print(f"[orchestrator] OUR-side canoes fine sweep took {t_ours:.1f}s")

    # ---- fastnc side, subprocess into the fastnc venv ---------------------
    print("=" * 78)
    print("ORCHESTRATOR: computing fastnc side via subprocess (fastnc venv)")
    print("=" * 78)
    t_fnc0 = time.time()
    proc = subprocess.run(
        [_FASTNC_PY, str(Path(__file__).resolve()),
         "--side", "fnc", "--out", str(fnc_tmp)],
        cwd=str(_HERE),
    )
    if proc.returncode != 0:
        raise SystemExit(f"fastnc subprocess failed (rc={proc.returncode})")
    t_fnc = time.time() - t_fnc0
    print(f"[orchestrator] fastnc-side sweep took {t_fnc:.1f}s")

    # ---- combine ----------------------------------------------------------
    ours = np.load(ours_tmp)
    fnc = np.load(fnc_tmp)
    # sanity: both sides used the exact same grid
    assert np.allclose(ours["phi_deg"], PHI_DEG_FINE)
    assert np.allclose(fnc["phi_deg"], PHI_DEG_FINE)
    assert np.allclose(ours["gamma_arcmin"], GAMMA_ARCMIN_FINE)
    assert np.allclose(fnc["gamma_arcmin"], GAMMA_ARCMIN_FINE)

    note = (
        "Appendix validation curve data (DATA ONLY, no figure), built "
        "2026-06-10. FINE grid: phi_deg = 24 pts evenly spaced in [10,120] deg "
        "(validated range, no extrapolation); gamma_arcmin = [10.0, 50.0] "
        "(150' skipped to save time); z_s=5, projection='cent'. "
        "ours = SFT zeta_D (JACFIX slot+centroid fix: apex->Gamma^1, "
        "Gamma^2=conj(Gamma^3), centroid-projection phase); "
        "fnc = fastnc tree-level SPT natural components "
        "(Lmax=24,Lmax_diag=48,Mmax=24,epmu=1e-7). Both sides use the SAME "
        "cosmology, pk_file (PCAMBz0.txt), and conventions pinned in "
        "stage1_grid.json. fi = sqrt(|G0|^2+|G1|^2+|G2|^2+|G3|^2) is the "
        "frame-invariant (projection-independent) 3PCF amplitude. ours via "
        "stage1_ours_v2_JACFIX compute path (born3=ones, vertex carries "
        "(1+z)^-4); fnc via stage1_fastnc_tree.build_bispectrum + eval_gammas."
    )

    np.savez(
        _FINAL_NPZ,
        phi_deg=PHI_DEG_FINE,
        gamma_arcmin=GAMMA_ARCMIN_FINE,
        ours_G0=ours["G0"], ours_G1=ours["G1"],
        ours_G2=ours["G2"], ours_G3=ours["G3"], ours_fi=ours["fi"],
        fnc_G0=fnc["G0"], fnc_G1=fnc["G1"],
        fnc_G2=fnc["G2"], fnc_G3=fnc["G3"], fnc_fi=fnc["fi"],
        note=note,
    )
    print(f"\n[orchestrator] saved FINAL -> {_FINAL_NPZ}")

    _report(ours, fnc, t_ours)
    return t_ours


def _report(ours, fnc, t_ours):
    print("\n" + "=" * 78)
    print("SANITY / SUMMARY")
    print("=" * 78)
    phi = PHI_DEG_FINE
    gammas = GAMMA_ARCMIN_FINE

    all_finite = True
    for arr_name, arr in (("ours_fi", ours["fi"]), ("fnc_fi", fnc["fi"])):
        if not np.all(np.isfinite(arr)):
            all_finite = False
            print(f"  !! NON-FINITE in {arr_name}")
    for side, d in (("ours", ours), ("fnc", fnc)):
        for k in ("G0", "G1", "G2", "G3"):
            if not np.all(np.isfinite(d[k])):
                all_finite = False
                print(f"  !! NON-FINITE in {side}_{k}")
    print(f"  finiteness: {'ALL FINITE' if all_finite else 'FAILED'}")

    # gamma=10' full table
    ig10 = int(np.argmin(np.abs(gammas - 10.0)))
    print(f"\n  gamma = {gammas[ig10]:.1f}'  full phi table:")
    print(f"    {'phi[deg]':>9} {'ours_fi':>14} {'fnc_fi':>14} {'ratio':>10}")
    for j in range(phi.size):
        o = ours["fi"][ig10, j]
        f = fnc["fi"][ig10, j]
        r = o / f if f != 0 else np.nan
        print(f"    {phi[j]:9.3f} {o:14.6e} {f:14.6e} {r:10.4f}")

    # phi=60 sanity (~0.7)
    j60 = int(np.argmin(np.abs(phi - 60.0)))
    r60 = ours["fi"][ig10, j60] / fnc["fi"][ig10, j60]
    print(f"\n  SANITY: gamma=10', phi={phi[j60]:.2f} (nearest to 60) "
          f"ours/fnc fi ratio = {r60:.4f}  (expect ~0.7)")

    # gamma=50' ratio range
    ig50 = int(np.argmin(np.abs(gammas - 50.0)))
    r50 = ours["fi"][ig50] / fnc["fi"][ig50]
    print(f"\n  gamma = {gammas[ig50]:.1f}'  ratio range: "
          f"[{np.min(r50):.4f}, {np.max(r50):.4f}]  median={np.median(r50):.4f}")
    print(f"\n  OUR-side canoes fine sweep timing: {t_ours:.1f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", choices=["ours", "fnc"], default=None,
                    help="compute ONE side and dump to --out (sub-mode); "
                         "omit to run the full orchestrator")
    ap.add_argument("--out", default=None, help="output npz for --side mode")
    args = ap.parse_args()

    if args.side is not None:
        if args.out is None:
            raise SystemExit("--out is required with --side")
        _run_side(args.side, Path(args.out))
        return 0

    _orchestrate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
