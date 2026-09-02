"""Cross-table replication: is the 1.000 one lucky number, or a kernel statement?

For each (vertex table, sft-wick fold) pair already on disk, bind the estimator
to EXACTLY the callable module the fold was run with, and compare.  The tables
differ in the bispectrum model (tree vs bihalofit), the ell cutoff (960 ->
15360), the radial refinement (r1/r2/r4) and the leg-permutation logic (legacy
vs perm-aware), so each pair weights the kernel differently and carries a
different amplitude.  A kernel error that happened to cancel against one table
cannot cancel against all of them.

Each row is: exact expectation of the estimator (sigma_lambda -> 0 by the NOTES
extrapolation) divided by that table's own fold, at the fold's own gamma nodes.
"""
import _env, importlib.util, sys, time
import numpy as np
from pathlib import Path
import fk_expect_exact as ex
from _fold import fold_nodes

PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/products")
PAIRS = [
    ("permclosed_cut15360 (PRODUCTION, perm-aware)",
     "_variant_table_permclosed_cut15360_permfix/callable__variant_table_permclosed_cut15360_permfix.py",
     "table_permclosed_cut15360_permfix_xi.npz"),
    ("permclosed_cut1000  (perm-aware)",
     "_variant_table_permclosed_cut1000_permfix/callable__variant_table_permclosed_cut1000_permfix.py",
     "table_permclosed_cut1000_permfix_xi.npz"),
    ("permclosed_cut1000_r1 (perm-aware, coarse radial)",
     "_variant_table_permclosed_cut1000_r1_permfix/callable__variant_table_permclosed_cut1000_r1_permfix.py",
     "table_permclosed_cut1000_r1_permfix_xi.npz"),
    ("permclosed_cut1000_r2 (perm-aware)",
     "_variant_table_permclosed_cut1000_r2_permfix/callable__variant_table_permclosed_cut1000_r2_permfix.py",
     "table_permclosed_cut1000_r2_permfix_xi.npz"),
    ("permclosed_cut1000  (LEGACY sorted-leg callable)",
     "_variant_table_permclosed_cut1000_legacy/callable__variant_table_permclosed_cut1000_legacy.py",
     "table_permclosed_cut1000_legacy_xi.npz"),
    ("tree_cut1000_r4     (legacy)",
     "_variant_table_tree_cut1000_r4/callable__variant_table_tree_cut1000_r4.py",
     "table_tree_cut1000_r4_xi.npz"),
    ("tree_cut15360_r4    (legacy)",
     "_variant_table_tree_cut15360_r4/callable__variant_table_tree_cut15360_r4.py",
     "table_tree_cut15360_r4_xi.npz"),
    ("bihalofit_cut960_r4 (legacy, NONLINEAR Bk)",
     "_variant_table_bihalofit_cut960_r4/callable__variant_table_bihalofit_cut960_r4.py",
     "table_bihalofit_cut960_r4_xi.npz"),
    ("bihalofit_cut15360_r4 (legacy, NONLINEAR Bk)",
     "_variant_table_bihalofit_cut15360_r4/callable__variant_table_bihalofit_cut15360_r4.py",
     "table_bihalofit_cut15360_r4_xi.npz"),
]

N = 1000
SIGMAS = (8.0, 4.0)
grid = ex.build_grid(n_lambda=N)
g0, _ = fold_nodes(PROD / PAIRS[0][2])
gsel = [float(g0[np.argmin(np.abs(g0 - t))]) for t in (1.0, 5.0, 17.3)]
print("fold gamma nodes:", [f"{g:.4f}'" for g in gsel])
print(f"\n{'table / fold pair':<46}" + "".join(f"{g:>9.2f}'" for g in gsel)
      + f"{'fold@1arcmin':>15}")
print("-" * 96)
out = []
for i, (label, cal, foldf) in enumerate(PAIRS):
    spec = importlib.util.spec_from_file_location(f"_vk3_{i}", PROD / cal)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[f"_vk3_{i}"] = mod
    spec.loader.exec_module(mod)
    ex._DS._k3 = mod
    gg, vv = fold_nodes(PROD / foldf)
    row = []
    for gs in gsel:
        j = int(np.argmin(np.abs(gg - gs)))
        fold = float(vv[j])
        cosg = float(np.cos(np.deg2rad(gg[j] / 60.0)))
        rr = {}
        for sig in SIGMAS:
            V, A, Q, rho = ex.node_stats(grid, cosg, sig)
            T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
            rr[sig] = float(((T1 + T2) + (T1 + T2).T)[0, 3]) / fold
        row.append(2 * rr[4.0] - rr[8.0])
    j0 = int(np.argmin(np.abs(gg - gsel[0])))
    print(f"{label:<46}" + "".join(f"{v:>10.4f}" for v in row)
          + f"{vv[j0]:>15.4e}", flush=True)
    out.append((label, row, float(vv[j0])))
a = np.array([r for _, rr, _ in out for r in rr])
print(f"\nall {a.size} (table, gamma) cells: median {np.median(a):.4f}  "
      f"min {a.min():.4f}  max {a.max():.4f}")
amps = np.array([v for _, _, v in out])
print(f"fold amplitude at 1' spans {amps.max()/amps.min():.2f}x across tables")
