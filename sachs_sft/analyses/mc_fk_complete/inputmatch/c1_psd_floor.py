"""C1. Does the PSD floor in _cholesky_psd make the MC field differ from the
field the exact Wick expectation assumes?

The MC draws innovations with L[k] = _cholesky_psd(V[k]); when V[k] has a
negative eigenvalue the floor makes L L^T = V_floored != V.  The exact
expectation propagates V as given.  So the two are NOT the same estimator.
Quantify: (a) how many nodes floor and how much kernel weight they carry,
(b) the exact expectation re-evaluated with the FLOORED variance recursion.
"""
from __future__ import annotations
import sys, time
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import fk_expect_exact as ex
_CORE = ex._CORE

GAM = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
SIG = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
NL = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
grid = ex.build_grid(n_lambda=NL)
V, A, Q, rho = ex.node_stats(grid, cosg, SIG, calibrate="smeared")
N = grid.lam.size

# --- (a) which nodes floor -------------------------------------------------
floored, negmag = [], []
Vf = np.empty_like(V)
for k in range(N):
    L = _CORE._cholesky_psd(V[k])
    Vf[k] = L @ L.T
    w = np.linalg.eigvalsh(V[k])
    if w[0] < 0:
        floored.append(k); negmag.append(w[0] / abs(w).max())
floored = np.array(floored, int)
kern = grid.dlam * grid.Wd * grid.Hd
kern_frac = kern[floored].sum() / kern.sum() if floored.size else 0.0
dV = np.abs(Vf - V).max(axis=(1, 2)) / np.abs(V).max(axis=(1, 2))
print(f"gamma={GAM}' sigma={SIG} N={N}")
print(f"  nodes with a negative eigenvalue in V: {floored.size}/{N}"
      f"   kernel-weight fraction {kern_frac:.4%}")
if floored.size:
    print(f"  worst relative negative eigenvalue {min(negmag):.3e}")
print(f"  max_k ||Vf-V||/||V||  = {dV.max():.3e}   (weighted mean "
      f"{np.average(dV, weights=kern/kern.sum()):.3e})")

# --- (b) exact expectation with the floored variance recursion -------------
Af = np.empty_like(Vf)
Af[0] = Vf[0]
for k in range(1, N):
    Af[k] = rho * rho * Af[k - 1] + (1.0 - rho * rho) * Vf[k]
t0 = time.time()
T1, T2 = ex.expectation(grid, cosg, SIG, V=V, A=A, Q=Q, rho=rho)
e_true = float(((T1 + T2) + (T1 + T2).T)[0, 3])
T1f, T2f = ex.expectation(grid, cosg, SIG, V=Vf, A=Af, Q=Q, rho=rho)
e_floor = float(((T1f + T2f) + (T1f + T2f).T)[0, 3])
print(f"  exact (A from V, what run_sweep prints) = {e_true:.7e}")
print(f"  exact (A from FLOORED V, what the MC actually simulates) = {e_floor:.7e}")
print(f"  ratio floored/true = {e_floor/e_true:.6f}   "
      f"[{time.time()-t0:.0f}s]")
np.savez(f"inputmatch/_c1_g{GAM}_s{SIG}_N{NL}.npz", floored=floored,
         kern_frac=kern_frac, e_true=e_true, e_floor=e_floor, dVmax=dV.max())
