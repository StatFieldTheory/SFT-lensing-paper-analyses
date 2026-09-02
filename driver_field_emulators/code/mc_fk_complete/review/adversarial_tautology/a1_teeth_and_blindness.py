"""ADVERSARIAL PROBE 1 -- is the gate-1 "placement identity" a tautology?

Two separate questions, both answered numerically on the real 6-component
two-ray system (no toy algebra):

  (A) TEETH.  Would ANY estimator pass?  Deliberately-wrong estimators are
      evaluated against the SAME discrete_truth and the deviation reported.
      If a wrong estimator passes at 1e-16 the gate is worthless.

  (B) BLINDNESS.  discrete_truth and expectation() share the kernels
      (response exponent, observable window, F vertex, Q, AR(1) statistics).
      Deform those CONSISTENTLY on both sides -- i.e. write down a different,
      physically WRONG stochastic system -- and ask whether the identity still
      holds.  If it does, the gate certifies the estimator's bookkeeping but
      says nothing about which diagram is being computed.

Nothing outside review/ is touched; the production modules are imported
read-only and the module globals that are monkeypatched (ex.F6, ex.FS6) are
restored.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                      # .../mc_fk_complete
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import fk_expect_exact as ex                                  # noqa: E402
from gate1_placement_share import discrete_truth              # noqa: E402


# --------------------------------------------------------------------------
# helpers: rebuild a Grid with a DIFFERENT (wrong) background / window
# --------------------------------------------------------------------------
def remake_grid(grid, D=None, w=None):
    """Rebuild every kernel of the Grid from a given D and node weight w.

    This reproduces build_grid() exactly, so the resulting Grid is a fully
    self-consistent discretisation of a (possibly WRONG) stochastic system:
    resp, Wd, Hd all follow from the same D and w.
    """
    N = grid.lam.size
    D = grid.D if D is None else np.asarray(D, float)
    if w is None:
        w = np.full(N, grid.dlam)
        w[0] = w[-1] = 0.5 * grid.dlam
    resp = np.empty(N)
    resp[0] = 1.0
    resp[1:] = (D[:-1] / D[1:]) ** 2
    tail = np.cumsum((w / D**2)[::-1])[::-1]
    Wd = D**2 * tail
    gj = grid.dlam * Wd[1:] / D[:-1] ** 4
    tail4 = np.zeros(N)
    tail4[:-1] = np.cumsum(gj[::-1])[::-1]
    Hd = D**4 * tail4
    return ex.Grid(grid.lam, grid.dlam, D, resp, Wd, Hd, grid.bg, grid.builder)


def kk(M):
    return float(M[0, 3])


def dev(new, truth):
    return float(np.abs(new - truth).max() / np.abs(truth).max())


# --------------------------------------------------------------------------
def main() -> int:
    N = 24
    sig = 60.0
    gam = 1.0
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=N, lam_min=406.0)
    V, A, Q, rho = ex.node_stats(grid, cosg, sig)
    print(f"grid N={N} dlam={grid.dlam:.2f} sigma={sig} gamma={gam}'  "
          f"sigma/dlam={sig/grid.dlam:.2f}")

    # ================= (A) TEETH ==========================================
    truth, Xi, R, U = discrete_truth(grid, A, Q, rho, ex.F6)
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)

    print("\n=== (A) TEETH: do WRONG estimators pass the same gate? ===")
    print(f"{'estimator':<46} {'kk/truth':>10} {'max rel dev':>12}  verdict")

    def report(name, M):
        d = dev(M, truth)
        print(f"{name:<46} {kk(M)/kk(truth):>10.5f} {d:>12.3e}  "
              f"{'PASS' if d < 1e-10 else 'FAIL'}")
        return d

    report("PRODUCTION  (T1+T2)+(T1+T2)^T", (T1 + T2) + (T1 + T2).T)
    report("T1 only  (simulate_fk_vr structure)", T1 + T1.T)
    report("T2 only", T2 + T2.T)
    report("T1+T2, no A<->B symmetrisation", T1 + T2)
    report("2*(T1+T2), double counted", 2 * ((T1 + T2) + (T1 + T2).T))
    report("(T1 + 2*T2) + transpose", (T1 + 2 * T2) + (T1 + 2 * T2).T)

    # wrong estimator variants that need re-running expectation()
    # (a) T2 built with F instead of FS  (no vertex-leg symmetrisation)
    FS_save = ex.FS6.copy()
    ex.FS6 = ex.F6.copy()
    _, T2_noFS = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    ex.FS6 = FS_save
    report("T2 with F instead of F+F^T on the vertex legs",
           (T1 + T2_noFS) + (T1 + T2_noFS).T)

    # (b) causality off-by-one: Phi[j] instead of Phi[j-1] in T1
    #     (i.e. F sees the state at its own node, not the previous one)
    T1_off = _t1_offbyone(grid, A, Q, rho)
    report("T1 with the F node reading its OWN state (off-by-one)",
           (T1_off + T2) + (T1_off + T2).T)

    # (c) drop the AR(1) cross-node correlation inside Phi (R -> diagonal)
    T1_diag, T2_diag = _expect_with_R(grid, A, Q, rho, diagonal=True)
    report("R[k,p] replaced by delta_kp A[k] (white-noise Phi)",
           (T1_diag + T2_diag) + (T1_diag + T2_diag).T)

    # ================= (B) BLINDNESS ======================================
    print("\n=== (B) BLINDNESS: change the PHYSICS consistently on both sides ===")
    print("Every row below is a DIFFERENT, deliberately wrong stochastic system.")
    print("The gate is re-run from scratch inside that wrong system.")
    print(f"{'wrong system':<46} {'kk (truth)':>12} {'max rel dev':>12}  verdict")

    def gate_in(name, g2, F6=None, Qx=None, Ax=None):
        F6u = ex.F6 if F6 is None else F6
        FS_bak, F_bak = ex.FS6.copy(), ex.F6.copy()
        if F6 is not None:
            ex.F6 = F6u
            ex.FS6 = F6u + F6u.transpose(0, 2, 1)
        Qu = Q if Qx is None else Qx
        Au = A if Ax is None else Ax
        tr, *_ = discrete_truth(g2, Au, Qu, rho, ex.F6)
        t1, t2 = ex.expectation(g2, cosg, sig, V=V, A=Au, Q=Qu, rho=rho)
        new = (t1 + t2) + (t1 + t2).T
        d = dev(new, tr)
        print(f"{name:<46} {kk(tr):>12.4e} {d:>12.3e}  "
              f"{'PASS' if d < 1e-10 else 'FAIL'}")
        ex.F6, ex.FS6 = F_bak, FS_bak
        return d

    gate_in("CONTROL: the true system", grid)

    # 1. wrong propagator power: (D_j/D_k)^1 instead of ^2, self-consistently
    gate_in("propagator (D_j/D_k)^1 instead of ^2",
            remake_grid(grid, D=grid.D ** 0.5))
    gate_in("propagator (D_j/D_k)^6 instead of ^2",
            remake_grid(grid, D=grid.D ** 3.0))
    # 2. no background at all: D = const  (flat space, no lensing kernel)
    gate_in("D = 1 (background switched off entirely)",
            remake_grid(grid, D=np.ones_like(grid.D)))
    # 3. nonsense observable window
    rng = np.random.default_rng(7)
    gate_in("random positive observable window w_k",
            remake_grid(grid, w=grid.dlam * (1.0 + rng.random(N))))
    # 4. arbitrary vertex F (NOT the paper's Sachs vertex, not even block diag)
    Frand = rng.standard_normal((6, 6, 6))
    gate_in("F replaced by a random (6,6,6) tensor", grid, F6=Frand)
    # 5. arbitrary Q (i.e. an arbitrary injected three-point cumulant)
    Qr = rng.standard_normal((N, 6, 6, 6))
    Qr = 0.5 * (Qr + Qr.transpose(0, 1, 3, 2))
    gate_in("Q replaced by a random symmetric tensor", grid, Qx=Qr)
    # 6. arbitrary node covariance A (i.e. wrong Sigma2 entirely)
    Ar = rng.standard_normal((N, 6, 6))
    Ar = np.einsum("kij,klj->kil", Ar, Ar)
    gate_in("Sigma2/A replaced by a random SPD sequence", grid, Ax=Ar)
    # 7. everything wrong at once
    gate_in("ALL of the above simultaneously",
            remake_grid(grid, D=grid.D ** 0.5,
                        w=grid.dlam * (1.0 + rng.random(N))),
            F6=Frand, Qx=Qr, Ax=Ar)
    return 0


# ---- wrong-estimator machinery -------------------------------------------
def _t1_offbyone(grid, A, Q, rho, block=64):
    """T1 with Phi[j] (state at the F node) instead of Phi[j-1]."""
    N = grid.lam.size
    dlam, Wd = grid.dlam, grid.Wd
    T1 = np.zeros((6, 6))
    cj1 = dlam * Wd
    for lo in range(0, N, block):
        ps = np.arange(lo, min(lo + block, N))
        Phi, R = ex._phi_block(grid, A, rho, ps)
        # shift back by one: Phi_off[j] = <s1[j] z_p>
        Phi_off = np.empty_like(Phi)
        Phi_off[:-1] = Phi[1:]
        Phi_off[-1] = Phi[-1]
        Ph = np.moveaxis(Phi_off, 1, 0).reshape(ps.size, N, 36)
        G = np.einsum("pjb,j,pjc->pbc", Ph, cj1, Ph, optimize=True)
        G = G.reshape(ps.size, 6, 6, 6, 6)
        T1 += np.einsum("p,Abc,pBuv,pbucv->AB",
                        dlam * Wd[ps], ex.F6, Q[ps], G, optimize=True)
    return T1


def _expect_with_R(grid, A, Q, rho, diagonal=False, block=64):
    """expectation() but with the AR(1) cross-node covariance replaced."""
    N = grid.lam.size
    dlam, Wd, D = grid.dlam, grid.Wd, grid.D
    T1 = np.zeros((6, 6)); T2 = np.zeros((6, 6))
    cj1 = dlam * Wd
    jm1 = np.maximum(np.arange(N) - 1, 0)
    Dj1 = D[jm1]
    jminus = np.arange(N) - 1
    idx = np.arange(N)
    for lo in range(0, N, block):
        ps = np.arange(lo, min(lo + block, N))
        lag = np.abs(idx[:, None] - ps[None, :])
        R = (rho ** lag)[:, :, None, None] * A[np.minimum(idx[:, None], ps[None, :])]
        if diagonal:
            R = np.where((lag == 0)[:, :, None, None], R, 0.0)
        phi = np.empty((N, ps.size, 6, 6))
        cur = np.zeros((ps.size, 6, 6))
        for k in range(N):
            cur = grid.resp[k] * cur + dlam * R[k]
            phi[k] = cur
        Phi = np.empty_like(phi); Phi[0] = 0.0; Phi[1:] = phi[:-1]
        Ph = np.moveaxis(Phi, 1, 0).reshape(ps.size, N, 36)
        G = np.einsum("pjb,j,pjc->pbc", Ph, cj1, Ph, optimize=True)
        G = G.reshape(ps.size, 6, 6, 6, 6)
        T1 += np.einsum("p,Abc,pBuv,pbucv->AB",
                        dlam * Wd[ps], ex.F6, Q[ps], G, optimize=True)
        U = np.where(jminus[:, None] >= ps[None, :],
                     dlam * (D[ps][None, :] / Dj1[:, None]) ** 2, 0.0)
        S = np.einsum("j,jp,jpbv->pbv", cj1, U, Phi, optimize=True)
        Psi = np.einsum("m,mpAu->pAu", dlam * Wd, R, optimize=True)
        T2 += np.einsum("Bbc,pcuv,pAu,pbv->AB", ex.FS6, Q[ps], Psi, S,
                        optimize=True)
    return T1, T2


if __name__ == "__main__":
    raise SystemExit(main())
