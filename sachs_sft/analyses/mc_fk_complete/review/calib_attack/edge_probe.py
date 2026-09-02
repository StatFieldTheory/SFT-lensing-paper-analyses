"""How much does the grid-edge behaviour of B matter?

NOTES admits B/Sigma2 falls to 0.35-0.56 at the first/last nodes.  Two tests:

(1) STATIONARY CONTINUATION.  Replace B by the value it would have if the
    AR(1) chain had a burn-in (an infinite stationary past) and continued past
    the source, i.e. add the two missing exponential tails.  This is a
    genuinely different, arguably more "natural" calibration; it removes the
    edge deficit entirely.  How far does FK move?

(2) EDGE EXCISION.  Zero Q on the first/last n nodes.  This bounds the TOTAL
    contribution those nodes make, edge-calibration error included.
"""
import sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/review/calib_attack")
import numpy as np, altcal, fold
import fk_expect_exact as ex

gam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0152
pairs = sys.argv[2] if len(sys.argv) > 2 else "1000:8,1999:4"
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
fv, _ = fold.fold_at(gam)
print(f"gamma={gam}'  fold={fv:.6e}")
for item in pairs.split(","):
    Ns, sigs = item.split(":"); N, sig = int(Ns), float(sigs)
    grid = ex.build_grid(n_lambda=N)
    V, A, Z, rho = altcal.base_stats(grid, cosg, sig)
    S2 = V * (2 * sig)
    B = altcal.calib_matrices(grid, V, A, rho, sig, "smeared")
    idx = np.arange(N)
    tailL = grid.dlam * (rho ** (idx + 1)) / (1.0 - rho)
    tailR = grid.dlam * (rho ** (N - idx)) / (1.0 - rho)
    Bsc = B + tailL[:, None, None] * A[0] + tailR[:, None, None] * A[N - 1]
    out = {}
    for nm, M in (("smeared", B), ("edgefix_statcont", Bsc)):
        Q = np.array([ex._DS.solve_Q(M[k], Z[k]) for k in range(N)])
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho,
                                block=128 if N <= 2000 else 64)
        out[nm] = float(((T1 + T2) + (T1 + T2).T)[0, 3])
    # edge excision, on the author's own B
    Q0 = np.array([ex._DS.solve_Q(B[k], Z[k]) for k in range(N)])
    ned = int(round(3 * sig / grid.dlam))
    for n in (1, ned, 2 * ned):
        Q = Q0.copy(); Q[:n] = 0.0; Q[N - n:] = 0.0
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho,
                                block=128 if N <= 2000 else 64)
        out[f"excise{n}"] = float(((T1 + T2) + (T1 + T2).T)[0, 3])
    base = out["smeared"]
    print(f"  N={N} sigma={sig} dlam={grid.dlam:.4f}  (3sigma = {ned} nodes)")
    r0 = np.array([np.abs(B[k]).max()/np.abs(S2[k]).max() for k in range(N)])
    r1 = np.array([np.abs(Bsc[k]).max()/np.abs(S2[k]).max() for k in range(N)])
    print(f"      B/S2   at k=0,1,2,N-2,N-1: "
          f"{r0[0]:.3f} {r0[1]:.3f} {r0[2]:.3f} ... {r0[-2]:.3f} {r0[-1]:.3f}")
    print(f"      Bsc/S2 at k=0,1,2,N-2,N-1: "
          f"{r1[0]:.3f} {r1[1]:.3f} {r1[2]:.3f} ... {r1[-2]:.3f} {r1[-1]:.3f}")
    for nm, v in out.items():
        print(f"      {nm:>18}  FK={v:+.6e}  /fold={v/fv:+.6f}  "
              f"/smeared={v/base:+.6f}", flush=True)
