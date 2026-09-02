"""F-VERTEX-FREE calibration test, for every candidate calibration.

<kappa1_a kappa1_b kappa1_c> at O(Q) must equal its local (white-noise) value
sum_k dlam Wd_k^3 zeta_sym[k].  This involves NO F vertex, NO FK kernel and no
comparison with the paper's fold, so it cannot have been tuned to the FK answer.

If a calibration passes THIS and still gives a different FK, the FK sigma->0
limit is calibration-dependent and the attack succeeds.
If every calibration that passes this gives the same FK, it does not.

Psi[k] = sum_l dlam Wd_l R[k,l] is built by O(N) recursions (no N^2 R tensor).
"""
import sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/review/calib_attack")
import numpy as np, altcal
import fk_expect_exact as ex


def psi_exact(grid, A, rho):
    """Psi[k] = sum_l dlam Wd_l rho^|k-l| A[min(k,l)]   -- O(N)."""
    N = grid.lam.size
    wd = grid.dlam * grid.Wd
    P = np.zeros((N, 6, 6))
    for k in range(1, N):
        P[k] = rho * (P[k - 1] + wd[k - 1] * A[k - 1])
    S = np.zeros(N)
    for k in range(N - 2, -1, -1):
        S[k] = rho * (S[k + 1] + wd[k + 1])
    return P + wd[:, None, None] * A + S[:, None, None] * A


def sym(Z):
    return (Z + Z.transpose(0, 2, 1, 3) * 0)  # placeholder, unused


def full_sym(Z):
    """Full 6-fold symmetrisation of a (N,6,6,6) tensor on its last 3 axes."""
    return (Z + Z.transpose(0, 1, 3, 2) + Z.transpose(0, 2, 1, 3)
            + Z.transpose(0, 2, 3, 1) + Z.transpose(0, 3, 1, 2)
            + Z.transpose(0, 3, 2, 1)) / 6.0


def run(gam, N, sig, names):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=N)
    V, A, Z, rho = altcal.base_stats(grid, cosg, sig)
    Psi = psi_exact(grid, A, rho)
    wd = grid.dlam * grid.Wd
    Zs = full_sym(Z)
    local = np.einsum("k,kabc->abc", grid.dlam * grid.Wd ** 3, Zs, optimize=True)
    out = {}
    for nm in names:
        M = altcal.calib_matrices(grid, V, A, rho, sig, nm)
        Q = np.array([ex._DS.solve_Q(M[k], Z[k]) for k in range(N)])
        t1 = np.einsum("k,kaij,kib,kjc->abc", wd, Q, Psi, Psi, optimize=True)
        tot = t1 + t1.transpose(1, 0, 2) + t1.transpose(2, 1, 0)
        out[nm] = (tot, local)
    return out


if __name__ == "__main__":
    gam = float(sys.argv[1]); pairs = sys.argv[2]; names = sys.argv[3].split(",")
    chans = [(0, 0, 0), (0, 0, 3), (0, 3, 3), (0, 1, 1), (3, 3, 3)]
    for item in pairs.split(","):
        Ns, sigs = item.split(":"); N, sig = int(Ns), float(sigs)
        out = run(gam, N, sig, names)
        print(f"\ngamma={gam}'  N={N}  sigma={sig}"
              f"   <k1^3> / local  (F-vertex-free)")
        hdr = "  ".join(f"{str(c):>10}" for c in chans)
        print(f"{'calib':>12}  {hdr}")
        for nm in names:
            tot, loc = out[nm]
            row = "  ".join(f"{tot[c]/loc[c]:>10.4f}" for c in chans)
            print(f"{nm:>12}  {row}", flush=True)
