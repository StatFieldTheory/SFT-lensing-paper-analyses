"""Is the sigma->0 = 1 result partly TAUTOLOGICAL?

Under the smeared calibration solve_Q is exact, so the injected cumulant is
identically the SYMMETRISED tabulated zeta at every sigma.  Split

   FK_exact(sigma)/fold  =  [FK_exact(sigma)/FK_localref]  x  [FK_localref/fold]

The first factor -> 1 as sigma->0 and is where the assembly is actually tested.
The second is sigma-independent: it is the discrete local fold of the injected
cumulant against the paper's analytic FK, i.e. the part the extrapolation does
NOT test.  Report it, and how much of it is the leg-symmetrisation of solve_Q.
"""
import sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/review/calib_attack")
import numpy as np, altcal, fold
import fk_expect_exact as ex
import k1cubed

N = 1999
grid = ex.build_grid(n_lambda=N)
print(f"N={N} dlam={grid.dlam:.4f}")
print(f"{'gamma':>10} {'fold':>13} {'ref_sym':>13} {'ref_raw':>13} "
      f"{'sym/raw':>8} {'sym/fold':>9} {'raw/fold':>9}")
for gam in (0.5005, 1.0152, 2.0621, 5.3042, 17.2755, 27.7051):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    Z = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    Zs = k1cubed.full_sym(Z)
    fv, _ = fold.fold_at(gam)
    r_sym = ex.fk_reference(grid, Zs)[0, 3]
    r_raw = ex.fk_reference(grid, Z)[0, 3]
    print(f"{gam:>10.4f} {fv:>13.6e} {r_sym:>13.6e} {r_raw:>13.6e} "
          f"{r_sym/r_raw:>8.4f} {r_sym/fv:>9.4f} {r_raw/fv:>9.4f}", flush=True)
