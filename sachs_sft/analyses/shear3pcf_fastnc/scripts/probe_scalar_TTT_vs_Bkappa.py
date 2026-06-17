#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""DECISIVE: is the canoes/STAGE-0 convergence normalisation == the B_kappa-phase
reference normalisation, in the SCALAR (spin-0) sector?

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

We discovered the canoes/STAGE-0 W_can scalar weight and the B_kappa-phase
reference scalar weight differ by a NET ~15-52x (per-shell, L-independent radial
factor; +0.44 ell-slope from chi-weighting on b_delta).  STAGE-0 TTT closes to
0.85 against the STAGE-0 SPT reference (same convention as canoes).  The spin-2
verdict compares against the B_kappa-phase reference (different convention).

THIS PROBE folds the SAME canoes zeta_TTT (deployed, HIGH+LOW) on the 4 SAS
configs and compares it to the SCALAR (spin-0) convergence 3PCF computed from the
B_kappa-phase reference itself (its natural_components with ALL spins set to 0).

  - If canoes_TTT / Bkappa_scalar ~ 15-52x  (matching the spin-2 4-17x scale),
    the discrepancy is a GLOBAL convention/normalisation between the two
    references, NOT a spin-2 bug; the spin-2 verdict is CONFOUNDED.
  - If canoes_TTT / Bkappa_scalar ~ 1, the scalar genuinely agrees with B_kappa
    and the spin-2 4-17x is a REAL spin-2-specific bug.

Writes outputs/scalar_TTT_vs_Bkappa.npz.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"
sys.path.insert(0, str(ANALYSIS))
from stage1_bkappa_phase_reference import (  # noqa: E402
    BkappaLimber, BkappaInterp, triangle_geometry,
)

sys.path.insert(0, "/Users/zzhang/projects/canoes")
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.sachs import (  # noqa: E402
    compute_kappa3_sigma3_high, compute_kappa3_zeta_table, kappa3_combine_low_high,
)
from scipy.integrate import cumulative_trapezoid  # noqa: E402

_N_LAMBDA = 40
_CONFIGS = [(10.0, 60.0), (50.0, 60.0)]


def scalar_convergence_3pcf(bk, gamma_rad, phi_rad,
                            n_l=160, n_ang=96, l_min=1e-1, l_max=2.0e4):
    """Spin-0 convergence 3PCF Z(X1,X2,X3) = INT B_kappa e^{i pos} (NO spin phase).

    Mirrors natural_components() with all spins = 0 (drops e^{2i...}); this is the
    scalar convergence 3-point function at the triangle, the J_0 analogue.
    """
    X1, X2, X3, d2, d3, a1, a2, a3 = triangle_geometry(gamma_rad, phi_rad)
    lr = np.geomspace(l_min, l_max, n_l)
    dlnl = np.gradient(np.log(lr))
    wl = lr * lr * dlnl
    pa = np.linspace(0.0, 2.0 * np.pi, n_ang, endpoint=False)
    dpa = 2.0 * np.pi / n_ang
    L2 = lr[:, None]
    P2 = pa[None, :]
    l2x = L2 * np.cos(P2)
    l2y = L2 * np.sin(P2)
    w2 = wl[:, None] * dpa
    pos2 = l2x * d2[0] + l2y * d2[1]
    norm = 1.0 / (2.0 * np.pi) ** 4
    acc = 0.0 + 0.0j
    for i3 in range(n_l):
        l3 = lr[i3]
        w3r = wl[i3] * dpa
        for j3 in range(n_ang):
            p3 = pa[j3]
            l3x = l3 * np.cos(p3)
            l3y = l3 * np.sin(p3)
            l1x = -(l2x + l3x)
            l1y = -(l2y + l3y)
            L1 = np.sqrt(l1x ** 2 + l1y ** 2)
            Bk = bk(L1, L2, l3)
            pos = pos2 + (l3x * d3[0] + l3y * d3[1])
            acc += np.sum(Bk * np.exp(1j * pos) * w2 * w3r * norm)
    return complex(acc)


def main():
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    bk_exact = BkappaLimber(n_chi=160)
    bk = BkappaInterp(bk_exact, n_l=96, n_c=64)

    # canoes deployed zeta_TTT fold (HIGH+LOW), same as STAGE-0 / stage1
    zf = np.linspace(0, 5, 6000)
    chif = np.array([chi_of_z(z, Omega_m=cosmo.Omega_m, h=cosmo.h) for z in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(5.0 / _N_LAMBDA, 5.0, _N_LAMBDA)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    lam_h = lam_phys * cosmo.h
    chi_s = float(chif[-1])
    a = 1.0 / (1.0 + z)
    K = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s

    common_low = dict(pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128,
                      spt_kind="tree_phi", units="physical",
                      lambda_convention="project", radial_discretization="sample",
                      radial_measure="lambda", parallelism="sequential",
                      shell_block_size=4, triple_chunk_size=4096, n_workers=1,
                      fuse_tree_phi_channels=True)
    common_high = dict(pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000,
                       ell_min=1e-3, n_ell=96, n_phi=64, units="physical",
                       lambda_convention="project", radial_discretization="sample",
                       radial_measure="lambda")

    print("=" * 100)
    print("SCALAR convergence 3PCF: canoes zeta_TTT fold  vs  B_kappa-phase reference (spins=0)")
    print("=" * 100)
    print(f"{'(gam,phi)':>14s} {'|canoes_TTT|':>14s} {'|Bk_scalar|':>14s} "
          f"{'can/Bk':>10s} {'sign_can':>9s} {'sign_Bk':>9s}")
    rows = []
    for gam_arcmin, phi_deg in _CONFIGS:
        gamma_rad = np.radians(gam_arcmin / 60.0)
        phi_rad = np.radians(phi_deg)
        t3 = 2.0 * gamma_rad * np.sin(phi_rad / 2.0)
        cg = np.cos(gamma_rad)
        triples = np.array([[np.cos(t3), cg, cg]], dtype=np.float64)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            low = compute_kappa3_zeta_table(triples, lam_h, channels=("TTT",), **common_low)
            high = compute_kappa3_sigma3_high(triples, lam_h, channels=("TTT",), **common_high)
            comb = kappa3_combine_low_high(low, high)
        zTTT = np.asarray(comb.zeta_TTT, dtype=np.float64)[0]
        # canoes TTT carries sign (-1)^3 from gamma=-INT; the scalar fold:
        can_TTT = -1.0 * float(np.trapezoid((K ** 3) * zTTT, lam_phys))
        bk_scalar = scalar_convergence_3pcf(bk, gamma_rad, phi_rad)
        ratio = abs(can_TTT) / abs(bk_scalar) if bk_scalar != 0 else float("nan")
        rows.append((gam_arcmin, phi_deg, can_TTT, bk_scalar, ratio))
        print(f"{f'({gam_arcmin:.0f},{phi_deg:.0f})':>14s} {abs(can_TTT):14.4e} "
              f"{abs(bk_scalar):14.4e} {ratio:10.4f} "
              f"{np.sign(can_TTT):9.0f} {np.sign(bk_scalar.real):9.0f}")
    print("=" * 100)
    ratios = np.array([r[4] for r in rows])
    print(f"median can/Bk = {np.median(ratios):.3f}, range [{ratios.min():.3f}, {ratios.max():.3f}]")
    print("INTERPRETATION:")
    print("  ~15-52x  => GLOBAL convention mismatch (scalar too); spin-2 verdict CONFOUNDED.")
    print("  ~4-17x   => same scale as the spin-2 discrepancy => convention, not spin-2 bug.")
    print("  ~1       => scalar agrees with B_kappa; the spin-2 4-17x is a REAL spin-2 bug.")

    np.savez(OUT / "scalar_TTT_vs_Bkappa.npz",
             configs=np.array(_CONFIGS),
             canoes_TTT=np.array([r[2] for r in rows]),
             Bk_scalar=np.array([complex(r[3]) for r in rows]),
             ratio=ratios,
             note="canoes zeta_TTT fold vs B_kappa-phase scalar (spins=0) convergence 3PCF.")
    print(f"saved -> {OUT / 'scalar_TTT_vs_Bkappa.npz'}")


if __name__ == "__main__":
    main()
