"""Per-shell redshift bookkeeping of the collapsed FK vertex, tested
against the real dense tables.

Three measurements, feeding notes/theory_redshift_dependence.md:

A.  The explicit redshift chain. The factorized prediction

        zeta_pred(gamma; shell) = jac * (2pi)^-4 * R_TTT(z) * chi^-4 * h^4
                                  * D(z)^4 * (2pi)^2 S_cap(gamma; chi) H(chi)

    with every prefactor read from the source (input_ready_b_to_zeta_spec
    section 2: per Phi00 leg A(a)(1+z)^4, A(a) = -1.5 Om H0^2 (1+z),
    jac = (1+z)^-4, growth^4 applied to the z=0 tree B), is compared
    shell by shell against the dense tree table. The shells span
    z = 0.10 to 5.70, so a single mis-read power of (1+z) would tilt the
    ratio by a factor of ~6 across the table; a flat ratio certifies the
    whole chain.

B.  The nonlinear redshift excess: zeta_bihalofit / zeta_tree per shell
    (floor-matched tables, same grid, same cutoff), and its slope against
    the tree's D(z)^4, i.e. the effective growth exponent
    g_eff = d ln zeta_NL / d ln D.

C.  The per-shell cutoff undercount U(chi) = H(chi; converged) /
    H(chi; ell = 1000) for the tree spectrum. Because the deployed cutoff
    is fixed in ELL, the k-space cap k_cut = 1000/chi varies with the
    shell, so U depends on redshift; the fold-weighted average of U then
    depends on the source redshift (script redshift_fold_trends.py).

Run (canoes venv, from driver_field_emulators/code):
    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python redshift_pershell.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.special import j0

from b_model import load_pk_lin
from b_model.cosmology import FIDUCIAL, growth

HERE = Path(__file__).resolve().parent
PRODUCTS = HERE.parent / "products"
OUT_DIR = PRODUCTS / "redshift"
OUT_DIR.mkdir(exist_ok=True)

# Project closed-form background (D = a*chi, lambda/z/chi maps).
_SCRIPTS = (HERE.parent.parent / "SFT-lensing-paper-analyses" / "sachs_sft"
            / "scripts")
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
import D_callable as Dc  # noqa: E402

PK = load_pk_lin()
ELL_CUT = 60.0
ELL_MAX = 1000.0
ARCMIN = np.pi / (180.0 * 60.0)
GAMMAS_A = np.array([0.5, 2.0, 5.0, 10.0, 20.0])   # inside validity window

H0_HMPC = 100.0 / 299792.458   # h/Mpc, spec convention


def s_cap(gamma_rad: float, chi_hmpc: float, ell_max: float,
          n: int = 4000) -> float:
    """Capped soft factor INT_0^L dv v J0(v gamma) P(v/chi), P at z=0."""
    v = np.geomspace(1e-2, ell_max, n)
    p = PK(v / chi_hmpc)
    return float(np.trapezoid(v * j0(v * gamma_rad) * p, v))


def hard_factor(chi_hmpc: float, ell_max: float, n: int = 3000) -> float:
    """INT_60^L du u (17/7 - n_eff/2) P(u/chi), P at z=0."""
    u = np.geomspace(ELL_CUT, ell_max, n)
    k = u / chi_hmpc
    p = PK(k)
    neff = np.gradient(np.log(p), np.log(k))
    rbar = 17.0 / 7.0 - neff / 2.0
    return float(np.trapezoid(u * rbar * p, u))


def hard_factor_converged(chi_hmpc: float, k_hi: float = 10.0) -> float:
    """Same, upper limit at fixed k = k_hi h/Mpc (past the n_eff = -2
    turnover of the LINEAR spectrum, where u^2 P per log-u decays)."""
    return hard_factor(chi_hmpc, k_hi * chi_hmpc)


def raw_vertex_z0(gamma_rad: float, chi_hmpc: float, ell_max: float,
                  n_uv: int = 400, n_phi: int = 64) -> float:
    """Brute-force z = 0 collapsed TTT quadrature at shell chi (no
    factorization): INT du dv dphi u v W B_tree 2pi J0(v gamma).

    Same integrand as the builder's HIGH branch before its prefactors,
    so pref * D^4 * raw must reproduce the table up to quadrature
    resolution; any redshift TREND in that ratio would be a bookkeeping
    error, with the squeezed approximation entirely out of the loop.
    """
    from b_model.tree import b_tree_z0

    grid = np.geomspace(0.1, 2.2 * ell_max, n_uv)
    logw = np.gradient(np.log(grid))
    x, w_phi = np.polynomial.legendre.leggauss(n_phi)
    phi = np.pi * (x + 1.0)
    wp = np.pi * w_phi
    u = grid[:, None, None]
    v = grid[None, :, None]
    ph = phi[None, None, :]
    w = np.sqrt(u ** 2 + v ** 2 + 2.0 * u * v * np.cos(ph))
    mx = np.maximum(np.maximum(u, v), w)
    window = (mx > ELL_CUT) & (mx <= ell_max)
    bisp = b_tree_z0(w / chi_hmpc, np.broadcast_to(u, w.shape) / chi_hmpc,
                     np.broadcast_to(v, w.shape) / chi_hmpc, PK)
    integrand = (u * v) * (u * v) * window * bisp * \
        (2.0 * np.pi) * j0(v * gamma_rad)          # extra u,v: log measure
    return float(np.einsum("ijk,i,j,k->", integrand, logw, logw, wp))


def collapsed_rows(table, gammas_arcmin):
    """Nearest stored (1, c, c) rows to the requested gammas.

    The dense grid stores the production gamma values, so we snap to the
    closest stored row and return its ACTUAL gamma alongside the index.
    """
    cts = np.asarray(table["cosine_triples"], float)
    g_rows = np.degrees(np.arccos(np.clip(cts[:, 1], -1, 1))) * 60.0
    idx, actual = [], []
    for g in gammas_arcmin:
        j = int(np.argmin(np.abs(g_rows - g)))
        idx.append(j)
        actual.append(g_rows[j])
    return np.asarray(idx), np.asarray(actual)


def main() -> None:
    tree = np.load(PRODUCTS / "collapsed_dense_tree_cut1000.npz")
    bih = np.load(PRODUCTS / "collapsed_dense_bihalofit_cut1000.npz")
    lam = np.asarray(tree["lambda_shells_Mpc"], float)
    z_sh = np.asarray(tree["z_shells"], float)
    chi_mpc = np.asarray(tree["chi_shells_Mpc"], float)     # physical Mpc
    h = FIDUCIAL.h
    om = FIDUCIAL.Omega_m
    chi_sh = chi_mpc * h                                    # Mpc/h (builder)

    # Background consistency: stored z(lambda), chi(lambda) vs D_callable.
    z_ref = np.asarray([float(Dc.z_of_lambda(x)) for x in lam])
    chi_ref = np.asarray([float(Dc.chi_of_lambda(x)) for x in lam])
    dz = np.max(np.abs(z_sh - z_ref) / (1 + z_ref))
    dchi = np.max(np.abs(chi_mpc - chi_ref) / chi_ref)
    print(f"[bg] stored vs D_callable: max dz/(1+z) = {dz:.2e}, "
          f"max dchi/chi = {dchi:.2e}")

    rows, g_actual = collapsed_rows(tree, GAMMAS_A)
    print("[rows] requested gammas:", GAMMAS_A,
          "-> stored:", np.round(g_actual, 3))
    gammas = g_actual
    zeta_tree = np.asarray(tree["zeta_TTT"], float)[rows]      # (n_g, n_sh)
    zeta_bih = np.asarray(bih["zeta_TTT"], float)[rows]

    # ---- A: explicit redshift chain --------------------------------------
    D4 = growth(z_sh) ** 4
    a_leg = -1.5 * om * H0_HMPC ** 2 * (1.0 + z_sh)            # A(a), spec
    r_ttt = (a_leg * (1.0 + z_sh) ** 4) ** 3                   # three legs
    jac = (1.0 + z_sh) ** (-4)                                 # (dlam/dchi)^2
    pref = jac * (2 * np.pi) ** (-4) * r_ttt * chi_sh ** (-4) * h ** 4

    S = np.empty((gammas.size, lam.size))
    Hh = np.empty(lam.size)
    U = np.empty(lam.size)
    for j, chi in enumerate(chi_sh):
        Hh[j] = hard_factor(chi, ELL_MAX)
        U[j] = hard_factor_converged(chi) / Hh[j]
        for i, g in enumerate(gammas):
            S[i, j] = s_cap(g * ARCMIN, chi, ELL_MAX)
    pred = pref[None, :] * D4[None, :] * (2 * np.pi) ** 2 * S * Hh[None, :]

    ratio = zeta_tree / pred
    print("\n[A] table / factorized prediction (tree, cutoff 1000):")
    hdr = "  z_shell " + "".join(f"  g={g:>4.1f}'" for g in gammas)
    print(hdr)
    for j in range(lam.size):
        cells = "".join(f"  {ratio[i, j]:7.3f}" for i in range(gammas.size))
        print(f"  {z_sh[j]:7.3f} {cells}")
    med = np.median(ratio, axis=1)
    spread = np.max(ratio, axis=1) / np.min(ratio, axis=1)
    print("  per-gamma median over shells:",
          " ".join(f"{m:.3f}" for m in med))
    print("  per-gamma max/min over shells (z 0.10 -> 5.03):",
          " ".join(f"{s:.3f}" for s in spread))
    print("  [a single mis-read power of (1+z) would give max/min ~ 5.5]")

    # ---- A2: bookkeeping certification with the raw quadrature -----------
    # Same prefactor chain, but the z = 0 quadrature is now brute force:
    # no squeezed factorization anywhere. A flat ratio = 1 certifies every
    # explicit redshift factor; the drift seen in [A] is then attributed to
    # the factorization error E, which grows where the ell window maps to
    # steep-spectrum k at the near shells.
    print("\n[A2] table / (pref * D^4 * raw z=0 quadrature)  "
          "[no factorization]:")
    print(hdr)
    ratio2 = np.empty_like(ratio)
    for j, chi in enumerate(chi_sh):
        for i, g in enumerate(gammas):
            raw = raw_vertex_z0(g * ARCMIN, chi, ELL_MAX)
            ratio2[i, j] = zeta_tree[i, j] / (pref[j] * D4[j] * raw)
        cells = "".join(f"  {ratio2[i, j]:7.3f}" for i in range(gammas.size))
        print(f"  {z_sh[j]:7.3f} {cells}")
    print("  per-gamma max/min over shells:",
          " ".join(f"{np.max(ratio2[i])/np.min(ratio2[i]):.3f}"
                   for i in range(gammas.size)))

    # ---- B: nonlinear redshift excess ------------------------------------
    rnl = zeta_bih / zeta_tree
    print("\n[B] bihalofit / tree per shell (floor-matched):")
    print(hdr)
    for j in range(lam.size):
        cells = "".join(f"  {rnl[i, j]:7.3f}" for i in range(gammas.size))
        print(f"  {z_sh[j]:7.3f} {cells}")
    # effective growth exponent between adjacent shells at gamma = 0.5'
    lnD = np.log(growth(z_sh))
    lnZ = np.log(np.abs(zeta_bih[0]))
    g_eff = np.gradient(lnZ, lnD)
    print("  effective d ln zeta_NL / d ln D at 0.5' per shell "
          "(tree would give 4 plus the geometric drift):")
    lnZt = np.log(np.abs(zeta_tree[0]))
    g_eff_tree = np.gradient(lnZt, lnD)
    for j in range(lam.size):
        print(f"  z = {z_sh[j]:6.3f}   NL {g_eff[j]:7.2f}   "
              f"tree {g_eff_tree[j]:7.2f}   excess {g_eff[j]-g_eff_tree[j]:6.2f}")

    # ---- C: per-shell cutoff undercount ----------------------------------
    print("\n[C] tree cutoff undercount U = H(k<=10 h/Mpc) / H(ell<=1000):")
    for j in range(lam.size):
        print(f"  z = {z_sh[j]:6.3f}  chi = {chi_sh[j]:7.1f} Mpc/h  "
              f"k_cut = {1000.0/chi_sh[j]:6.3f} h/Mpc  U = {U[j]:6.2f}")

    np.savez(
        OUT_DIR / "pershell.npz",
        lam=lam, z=z_sh, chi_hmpc=chi_sh,
        gammas_arcmin=gammas,
        ratio_pred=ratio, ratio_raw=ratio2, rnl=rnl, U=U, H1000=Hh,
        D=growth(z_sh),
    )
    print(f"\n[saved] {OUT_DIR / 'pershell.npz'}")


if __name__ == "__main__":
    main()
