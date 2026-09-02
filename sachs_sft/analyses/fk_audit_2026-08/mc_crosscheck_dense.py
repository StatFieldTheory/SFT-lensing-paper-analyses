"""Independent cross-check of the dense-grid FK curve.

`analyses/mc_sachs_2pt/fk_analytic.py` re-assembles the FK two-point
diagram with its own quadrature: its own causal simplex, its own response
product and its own F-K index contraction. It reads the SAME zeta table as
the production fold, so what it validates is the ASSEMBLY, not the input
statistics. That is exactly the check the dense-grid rebuild needs: the
table changed, and we want to know that sft-wick still folds it correctly.

This script points both sides at a chosen table and compares.

Usage (from ``sachs_sft/analyses/fk_audit_2026-08``)::

    SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
    $SFT mc_crosscheck_dense.py ../products/collapsed_dense_tree_cut1000.npz \
                                ../products/collapsed_dense_tree_cut1000_xi.npz
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

SS = (Path(__file__).resolve().parents[2])
MC = SS / "analyses" / "mc_sachs_2pt"
K3 = SS / "callables" / "kappa3_vertex" / "equal_time_limber"
for p in (str(MC), str(K3), str(SS / "callables" / "C_propagator" / "corr_op")):
    if p not in sys.path:
        sys.path.insert(0, p)


def load_xi(path: Path):
    """(gamma_arcmin, cos_gamma, xi_kk, t_final) from a sweep result."""
    d = np.load(path, allow_pickle=True)  # locally produced artifact
    m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
    xs, ys = d["x"][m], d["y"][m]
    vs = np.asarray(d["value"], float)[m]
    t_final = float(np.unique(d["t_final"])[0])
    gam, cos = [], []
    for xi, yi in zip(xs, ys):
        xi = np.asarray(xi, float) / np.linalg.norm(np.asarray(xi, float))
        yi = np.asarray(yi, float) / np.linalg.norm(np.asarray(yi, float))
        c = float(np.clip(xi @ yi, -1.0, 1.0))
        cos.append(c)
        gam.append(math.degrees(math.acos(c)) * 60.0)
    order = np.argsort(gam)
    return (np.asarray(gam)[order], np.asarray(cos)[order], vs[order], t_final)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("table", type=Path, help="kappa3 NPZ under test")
    ap.add_argument("xi", type=Path, help="the FK sweep result folded from it")
    args = ap.parse_args()

    # Redirect the kappa3 callable at the table under test BEFORE the MC
    # modules import it, so both sides read the same statistics.
    import equal_time_limber_kappa3_callable as k3

    k3.TABLE_PATH = args.table.resolve()
    if hasattr(k3, "_CACHE"):
        k3._CACHE = None
    for name in ("_cache", "_TABLE_CACHE"):
        if hasattr(k3, name):
            setattr(k3, name, None)

    import fk_analytic  # noqa: E402  (after the redirect)

    gam, cos, xi_prod, t_final = load_xi(args.xi)
    print(f"table : {args.table.name}")
    print(f"fold  : {args.xi.name}   (t_final = {t_final:.4f})")
    print()
    print(f"{'gamma[arcmin]':>13} {'independent':>14} {'production':>14} {'ratio':>9}")
    ratios = []
    for i, (g, c) in enumerate(zip(gam, cos)):
        if g > 700.0:  # beyond this both sides are at the cancellation floor
            continue
        indep = float(fk_analytic.fk_analytic(c, t_final)[0, 0])
        prod = xi_prod[i]
        r = indep / prod if prod != 0 else np.nan
        ratios.append(r)
        if i % 4 == 0 or i < 3:
            print(f"{g:13.3f} {indep:14.5e} {prod:14.5e} {r:9.4f}")
    ratios = np.asarray(ratios)
    good = np.isfinite(ratios) & (ratios > 0)
    print()
    print(f"  median ratio = {np.median(ratios[good]):.4f}   "
          f"spread = {np.std(ratios[good]):.4f}   "
          f"range [{ratios[good].min():.4f}, {ratios[good].max():.4f}]")
    print("  (a flat ratio near 1 means the production fold reproduces an")
    print("   independent assembly of the same diagram from the same table)")


if __name__ == "__main__":
    main()
