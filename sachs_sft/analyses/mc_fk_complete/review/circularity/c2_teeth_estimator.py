"""TEETH, injected on the ESTIMATOR side (not on the reference kernel).

``verify_teeth.py`` in the parent folder perturbs the discrete LOCAL REFERENCE
kernel.  That is a proxy for the analytic side, and it leaves the object the
paper actually claims to have validated -- the stochastic-dynamics estimator --
untouched.  Here every defect is injected into the estimator itself (the exact
Wick expectation of the Monte-Carlo, which the MC reproduces to 0.2%), and the
comparison target is the UNMODIFIED paper fold.

Reported: exact/fold at sigma_lambda = 8 and 4, and the NOTES extrapolation
2 r(4) - r(8) to sigma_lambda -> 0.
"""
import _env, time, numpy as np, fk_expect_exact as ex, inj_exact as ij
from run_gate2 import analytic_fk

N = 1000
GAMMAS = (1.0, 5.0, 17.3)
SIGMAS = (8.0, 4.0)
DEFECTS = ["none", "acausal_U", "acausal_Phi", "acausal_both", "U_pow1",
           "F_no_leg_sym", "F_coef_half", "F_legswap", "no_transpose",
           "zeta_x1.3", "drop_dlam_T2"]

grid = ex.build_grid(n_lambda=N)
grid_p1 = ij.build_grid(n_lambda=N, resp_pow=1.0)
grid_nw = ij.build_grid(n_lambda=N, end_weights=False)

res = {}
for gam in GAMMAS:
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    ana = float(analytic_fk([gam])[0])
    for sig in SIGMAS:
        t0 = time.time()
        V, A, Q, rho = ex.node_stats(grid, cosg, sig)
        for d in DEFECTS:
            T1, T2 = ij.expectation(grid, cosg, sig, V, A, Q, rho, defect=d)
            res[(gam, sig, d)] = ij.fk_kk(T1, T2, d) / ana
        # grid-level defects: same nodes (V,A,Q do not depend on resp/weights)
        T1, T2 = ij.expectation(grid_p1, cosg, sig, V, A, Q, rho)
        res[(gam, sig, "resp_pow1")] = ij.fk_kk(T1, T2) / ana
        T1, T2 = ij.expectation(grid_nw, cosg, sig, V, A, Q, rho)
        res[(gam, sig, "no_end_weights")] = ij.fk_kk(T1, T2) / ana
        print(f"  [{gam}' sigma={sig}] {time.time()-t0:.0f}s", flush=True)

ALL = DEFECTS + ["resp_pow1", "no_end_weights"]
LBL = {"none": "BASELINE (no defect)",
       "acausal_U": "causality: F node may precede the drive (T2 leg)",
       "acausal_Phi": "causality: F sees the NEW state (s1[j] not s1[j-1])",
       "acausal_both": "causality: both of the above",
       "U_pow1": "one propagator on the F leg (D^1 not D^2 in U)",
       "F_no_leg_sym": "F+F^T -> F on the F-vertex legs (T2)",
       "F_coef_half": "F_101 = F_202 : -2 -> -1  (lost factor 2)",
       "F_legswap": "transpose the two F input legs",
       "no_transpose": "omit the A<->B symmetrisation",
       "zeta_x1.3": "scale zeta by 1.3 on the ESTIMATOR SIDE ONLY",
       "drop_dlam_T2": "drop one dlam measure in T2",
       "resp_pow1": "response (D_/D)^2 -> ^1 (whole dynamics)",
       "no_end_weights": "drop the trapezoid end weights"}

print("\nexact expectation / paper fold, sigma_lambda -> 0 (2 r(4) - r(8))\n")
hdr = f"{'defect':<50}" + "".join(f"{g:>10.1f}'" for g in GAMMAS)
print(hdr); print("-" * len(hdr))
out = {}
for d in ALL:
    row = []
    for gam in GAMMAS:
        r8, r4 = res[(gam, 8.0, d)], res[(gam, 4.0, d)]
        row.append(2 * r4 - r8)
    out[d] = row
    print(f"{LBL[d]:<50}" + "".join(f"{v:>11.4f}" for v in row))
print("\n(raw, no extrapolation)  sigma=8 / sigma=4")
for d in ALL:
    s = "".join(f"  {res[(g,8.0,d)]:.4f}/{res[(g,4.0,d)]:.4f}" for g in GAMMAS)
    print(f"{LBL[d]:<50}{s}")
np.savez("_c2_teeth.npz", keys=np.array([f"{k[0]}|{k[1]}|{k[2]}" for k in res]),
         vals=np.array(list(res.values())))
