"""Direct numerical test of the collapsed-vertex factorization theorem.

Tests Proposition (3.1) of notes/theory_shape_amplitude_factorization.md
against the real linear P(k), by brute-force quadrature of the raw vertex
integral. First run in the review round of 2026-08-25; kept as a permanent
regression so the note's verification status stays reproducible.

    raw(gamma; L) = INT du dv dphi  u v  W  B_tree(w/chi, u/chi, v/chi)
                    x 2 pi J_0(v gamma)
    prediction:  raw = (2 pi)^2 S_cap(gamma; L) H(L) [1 + E]
    S_cap(gamma; L) = INT_0^L dv v J_0(v gamma) P(v/chi)
    H(L) = INT_60^L du u (17/7 - n(u/chi)/2) P(u/chi)

The soft factor carries the SAME window cap L as the hard one: the J_0
tail decays only as an oscillatory v^{-1} envelope, so S_cap approaches
the uncapped S only when gamma L >> 10. That correction to the note's
first draft came out of this test.

Run (canoes venv):
    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python verify_factorization_numeric.py
"""

from __future__ import annotations

import numpy as np
from scipy.special import j0

from b_model import load_pk_lin
from b_model.tree import b_tree_z0

CHI_H = 4445.0        # mid shell, Mpc/h
ELL_CUT = 60.0
ARCMIN = np.pi / (180.0 * 60.0)
GAMMAS = np.array([0.5, 1.0, 5.0, 10.0, 20.0, 40.0, 80.0, 150.0])

PK = load_pk_lin()


def raw_vertex(gamma_rad: float, ell_max: float,
               n_uv: int = 400, n_phi: int = 64) -> float:
    """Brute-force quadrature of the collapsed TTT vertex integrand."""
    grid = np.geomspace(0.1, 2.2 * ell_max, n_uv)
    logw = np.gradient(np.log(grid))
    x, w_phi = np.polynomial.legendre.leggauss(n_phi)
    phi = np.pi * (x + 1.0)
    wp = np.pi * w_phi
    u = grid[:, None, None]
    v = grid[None, :, None]
    ph = phi[None, None, :]
    w = np.sqrt(u**2 + v**2 + 2.0 * u * v * np.cos(ph))
    window = (np.maximum(np.maximum(u, v), w) > ELL_CUT) & \
             (np.maximum(np.maximum(u, v), w) <= ell_max)
    bisp = b_tree_z0(w / CHI_H, np.broadcast_to(u, w.shape) / CHI_H,
                     np.broadcast_to(v, w.shape) / CHI_H, PK)
    integrand = (u * v) * (u * v) * window * bisp * \
        (2.0 * np.pi) * j0(v * gamma_rad)          # extra u,v: log measure
    return float(np.einsum(
        "ijk,i,j,k->", integrand, logw, logw, wp / 1.0))


def soft_cap(gamma_rad: float, ell_max: float, n: int = 20000) -> float:
    v = np.geomspace(1e-2, ell_max, n)
    p = np.asarray([PK(np.array([vv / CHI_H]))[0] for vv in
                    v[:: max(1, n // 2000)]])
    vv = v[:: max(1, n // 2000)]
    integrand = vv * j0(vv * gamma_rad) * p
    return float(np.trapezoid(integrand, vv))


def hard_factor(ell_max: float, n: int = 2000) -> float:
    u = np.geomspace(ELL_CUT, ell_max, n)
    k = u / CHI_H
    p = np.asarray([PK(np.array([kk]))[0] for kk in k])
    lnp = np.log(p)
    neff = np.gradient(lnp, np.log(k))
    rbar = 17.0 / 7.0 - neff / 2.0
    return float(np.trapezoid(u * rbar * p, u))


def main() -> None:
    for ell_max in (1000.0, 4000.0):
        h = hard_factor(ell_max)
        print(f"\n=== ell_max = {ell_max:.0f}   H = {h:.4e} ===")
        print("  gamma[']   raw           (2pi)^2 S_cap H   ratio")
        for g_arcmin in GAMMAS:
            g = g_arcmin * ARCMIN
            raw = raw_vertex(g, ell_max)
            pred = (2.0 * np.pi) ** 2 * soft_cap(g, ell_max) * h
            print(f"  {g_arcmin:8.1f}  {raw:+.5e}   {pred:+.5e}   "
                  f"{raw / pred:7.3f}")

    # sign-change locations, ell_max = 1000
    ell_max = 1000.0
    gg = np.linspace(60.0, 110.0, 26)
    raws = np.array([raw_vertex(g * ARCMIN, ell_max) for g in gg])
    scap = np.array([soft_cap(g * ARCMIN, ell_max) for g in gg])

    def zero_of(y):
        s = np.where(np.diff(np.sign(y)))[0]
        if s.size == 0:
            return float("nan")
        i = s[0]
        return gg[i] - y[i] * (gg[i + 1] - gg[i]) / (y[i + 1] - y[i])

    print(f"\nSign change at ell_max=1000: raw at {zero_of(raws):.1f} arcmin, "
          f"S_cap at {zero_of(scap):.1f} arcmin")


if __name__ == "__main__":
    main()
