"""Is FK immune to the driving field's TWO-POINT structure, not merely to its scale?

The established fact is invariance under Sigma2 -> c*Sigma2 (a scale).  That is
weaker than what the input-matching question needs.  Here the exact expectation
engine is fed a Sigma2 with a DIFFERENT lambda profile and a DIFFERENT angular
(two-ray) structure, with zeta6 held FIXED, and the FK entry [0,3] compared.

Nothing in fk_expect_exact.py is modified: ``expectation`` already accepts
(V, A, Q, rho), so a local re-implementation of ``node_stats`` with a transform
applied to Sigma2 is enough.
"""
from __future__ import annotations

import sys
import numpy as np
sys.path.insert(0, "..")
import fk_expect_exact as ex

_DS, _BGM, _CORE = ex._DS, ex._BG, ex._CORE


def node_stats_mod(grid, cos_gamma, sigma_lambda, transform=None,
                   calibrate="smeared"):
    """node_stats with an arbitrary transform applied to Sigma2; zeta6 untouched."""
    N = grid.lam.size
    sig = float(sigma_lambda)
    V = np.empty((N, 6, 6)); Z = np.empty((N, 6, 6, 6))
    for k, lk in enumerate(grid.lam):
        S2 = _CORE._assemble_6x6(grid.builder, float(cos_gamma), float(lk))
        if transform is not None:
            S2 = transform(float(lk), S2, grid)
        V[k] = S2 / (2.0 * sig)
        Z[k] = _DS.zeta6(float(cos_gamma), float(lk))
    rho = float(np.exp(-grid.dlam / sig))
    A = np.empty_like(V); A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1.0 - rho * rho) * V[k]
    M = (V * (2.0 * sig)) if calibrate == "nominal" else ex.smeared_covariance(grid, A, rho)
    Q = np.array([_DS.solve_Q(M[k], Z[k]) for k in range(N)])
    return V, A, Q, rho, M, Z


def fk_exact(grid, cos_gamma, sigma, transform=None):
    V, A, Q, rho, M, Z = node_stats_mod(grid, cos_gamma, sigma, transform)
    T1, T2 = ex.expectation(grid, cos_gamma, sigma, V=V, A=A, Q=Q, rho=rho)
    S = (T1 + T2) + (T1 + T2).T
    return S[0, 3], T1, T2, M, Z


# ---- the transforms (zeta always held fixed) -------------------------------
def t_id(l, S, g): return S


def t_scale2(l, S, g): return 2.0 * S


def t_ramp(l, S, g):
    u = (l - g.lam[0]) / (g.lam[-1] - g.lam[0])
    return (0.4 + 3.2 * u) * S            # 8x swing in the lambda profile


def t_osc(l, S, g):
    u = (l - g.lam[0]) / (g.lam[-1] - g.lam[0])
    return (1.0 + 0.7 * np.sin(6.0 * np.pi * u)) * S


def t_nocross(l, S, g):
    out = S.copy(); out[0:3, 3:6] = 0.0; out[3:6, 0:3] = 0.0
    return out                            # rays uncorrelated at 2-point level


def t_diagonly(l, S, g):
    return np.diag(np.diag(S)).copy()     # kill every component correlation


def t_isotropic(l, S, g):
    """Replace Sigma2 by (tr S / 6) * I : same total power, no structure at all."""
    return (np.trace(S) / 6.0) * np.eye(6)


TRANSFORMS = [("identity", t_id), ("2 x Sigma2", t_scale2),
              ("lambda ramp 0.4->3.6", t_ramp), ("lambda oscillation +-70%", t_osc),
              ("cross-ray block -> 0", t_nocross), ("component-diagonal only", t_diagonly),
              ("isotropic  (tr/6) I", t_isotropic)]


if __name__ == "__main__":
    import time
    NL = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    GAM = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
    grid = ex.build_grid(n_lambda=NL)
    sigmas = [16.0, 8.0, 4.0, 2.0]
    base = {}
    print(f"# gamma = {GAM}'  n_lambda = {NL}   FK[0,3] exact expectation")
    print(f"# {'transform':<26}" + "".join(f"{'s=%g' % s:>14}" for s in sigmas)
          + "".join(f"{'ratio s=%g' % s:>12}" for s in sigmas))
    out = {}
    for name, fn in TRANSFORMS:
        vals = []
        for s in sigmas:
            t0 = time.time()
            v, T1, T2, M, Z = fk_exact(grid, cosg, s, fn)
            vals.append(v)
        out[name] = vals
        if name == "identity":
            base = dict(zip(sigmas, vals))
        rat = [vals[i] / base[sigmas[i]] for i in range(len(sigmas))]
        print(f"  {name:<26}" + "".join(f"{v:14.6e}" for v in vals)
              + "".join(f"{r:12.6f}" for r in rat))
    np.savez("fk_c_invariance.npz", names=np.array([n for n, _ in TRANSFORMS]),
             sigmas=np.array(sigmas),
             vals=np.array([out[n] for n, _ in TRANSFORMS]))
