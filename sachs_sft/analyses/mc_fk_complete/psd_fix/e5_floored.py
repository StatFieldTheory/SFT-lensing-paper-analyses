"""Does the eigen-floor make the Monte-Carlo realise a DIFFERENT covariance
from the one the exact expectation assumes?

``fk_complete_core`` draws its innovations with ``Lf[k] = _cholesky_psd(V[k])``,
so the covariance the simulation actually realises is ``Lf Lf^T`` -- the
eigen-FLOORED V.  ``node_stats`` builds the AR(1) variance ``A`` from the raw
``V``, negative eigenvalues and all.  Where V is indefinite the two disagree.
This re-runs the exact expectation with V replaced by ``Lf Lf^T`` to size that.
"""
from __future__ import annotations
import sys
import numpy as np
HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, HERE); sys.path.insert(0, HERE + "/psd_fix")
import fk_expect_exact as ex
import sigma2_repaired as R
from run_gate2 import analytic_fk
_CORE, _DS = ex._CORE, ex._DS

SIG, N = 8.0, 1000


def exact_with(grid, cosg, sig, V, floor):
    if floor:
        L = np.array([_CORE._cholesky_psd(V[k]) for k in range(V.shape[0])])
        V = np.einsum("kij,klj->kil", L, L, optimize=True)
    rho = float(np.exp(-grid.dlam / sig))
    A = np.empty_like(V); A[0] = V[0]
    for k in range(1, V.shape[0]):
        A[k] = rho * rho * A[k - 1] + (1 - rho * rho) * V[k]
    M = ex.smeared_covariance(grid, A, rho)
    Z = np.array([_DS.zeta6(cosg, float(l)) for l in grid.lam])
    Q = np.array([_DS.solve_Q(M[k], Z[k]) for k in range(V.shape[0])])
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    return float(((T1 + T2) + (T1 + T2).T)[0, 3])


for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    fold = float(analytic_fk([gam])[0])
    print(f"=== gamma={gam}'  sigma={SIG}  N={N}  fold={fold:.6e} ===")
    for tag, cls in (("baseline", None), ("R1Knots", R.R1Knots)):
        grid = ex.build_grid(n_lambda=N)
        if cls is not None:
            grid.builder = cls(background=grid.bg, apply_c0=False)
        V = np.array([_CORE._assemble_6x6(grid.builder, cosg, float(l))
                      for l in grid.lam]) / (2.0 * SIG)
        raw = exact_with(grid, cosg, SIG, V, False)
        flo = exact_with(grid, cosg, SIG, V, True)
        print(f"  {tag:<9} exact(raw V)     = {raw:.7e}  ({raw/fold:.5f} of fold)")
        print(f"  {'':<9} exact(floored V) = {flo:.7e}  ({flo/fold:.5f} of fold)"
              f"   floored/raw = {flo/raw:.7f}")
