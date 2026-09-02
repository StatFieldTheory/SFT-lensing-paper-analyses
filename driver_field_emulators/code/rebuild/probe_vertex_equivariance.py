"""Is the three-point vertex equivariant under simultaneous leg/direction swaps?

The cumulant <X_a(n_i) X_b(n_j) X_c(n_k)> is symmetric under permuting the
(leg, direction) pairs together.  The Monte-Carlo's ``driver_stats.zeta6``
relies on that: it evaluates the vertex once per ORDERED ray triple and drops
the value into the (6,6,6) tensor at the legs' component indices, asserting in
its docstring that the result is "symmetric by construction".

``solve_Q`` then fully symmetrises whatever it is handed, so any failure of
equivariance is silently averaged away inside the Monte-Carlo while the
analytic fold keeps it.  The two sides would then be using different vertices.

This probe evaluates the same isoceles triangle in all three leg orderings and
reports how far apart they are.  The collapsed family (1, cos g, cos g) is the
one the FK two-point function queries, and its three orderings are

    (n1,n1,n2) -> (1, cos g, cos g)
    (n1,n2,n1) -> (cos g, cos g, 1)
    (n2,n1,n1) -> (cos g, 1, cos g)

which are the same triangle written three ways.  If the table stores only one
of the three, the other two are reached by nearest-neighbour interpolation and
need not agree.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"
_FIX = _REPO / "driver_field_emulators" / "code" / "callable_fixed"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--callable", choices=("fixed", "deployed"), default="fixed")
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.5, 1.0, 2.6, 6.7, 17.3, 44.4, 114.3, 293.9, 957.2])
    ap.add_argument("--lams", type=float, nargs="+", default=[800.0, 1500.0, 2100.0])
    a = ap.parse_args()

    sys.path.insert(0, str(_MC))
    sys.path.insert(0, str(_FIX))
    import driver_stats as ds
    if a.callable == "fixed":
        import perm_aware_kappa3_callable as k3
        k3.TABLE_PATH = a.table.resolve()
        k3._CACHE = None
    else:
        sys.path.insert(0, str(
            _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables"
            / "kappa3_vertex" / "equal_time_limber"))
        import equal_time_limber_kappa3_callable as k3
        k3.TABLE_PATH = a.table.resolve()
    ds._k3 = k3

    n1 = ds._N1
    print(f"[probe] callable={a.callable}, table={a.table.name}\n")
    # Split by how many legs carry spin.  The all-scalar entry (Phi00 on every
    # leg) has no basis convention at all, so if IT is non-equivariant the fault
    # is the table lookup; if only the spin entries move, the fault is the
    # spin-basis convention under leg permutation.
    spinclass = {0: [(0, 0, 0)], 1: [], 2: [], 3: []}
    for a3 in range(3):
        for b3 in range(3):
            for c3 in range(3):
                n = sum(1 for x in (a3, b3, c3) if x > 0)
                if n:
                    spinclass[n].append((a3, b3, c3))
    print(f"{'gamma':>8} {'lam':>7} {'asym all':>10} | "
          f"{'0 spin':>10} {'1 spin':>10} {'2 spin':>10} {'3 spin':>10}")
    for gamma in a.gammas:
        cg = float(np.cos(np.deg2rad(gamma / 60.0)))
        n2 = ds._n2_of_cos(cg)
        for lam in a.lams:
            A = ds.zeta_tensor(n1, n1, n2, lam)     # legs (a,b,c) at (n1,n1,n2)
            B = ds.zeta_tensor(n1, n2, n1, lam)     # legs (a,c,b) at (n1,n1,n2)
            C = ds.zeta_tensor(n2, n1, n1, lam)     # legs (c,a,b) at (n1,n1,n2)
            # equivariance: A[a,b,c] == B[a,c,b] == C[c,a,b]
            Bp = np.transpose(B, (0, 2, 1))
            Cp = np.transpose(C, (1, 2, 0))
            nA = np.linalg.norm(A)
            dB = np.linalg.norm(A - Bp) / max(nA, 1e-300)
            dC = np.linalg.norm(A - Cp) / max(nA, 1e-300)
            # the piece solve_Q would discard: A minus its full symmetrisation
            S = sum(np.transpose(A, p) for p in
                    [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1),
                     (2, 1, 0)]) / 6.0
            asym = np.linalg.norm(A - S) / max(nA, 1e-300)
            per = {}
            for n, cells in spinclass.items():
                idx = tuple(np.array(cells).T)
                num = np.linalg.norm((A - S)[idx])
                den = np.linalg.norm(A[idx])
                per[n] = num / max(den, 1e-300)
            print(f"{gamma:8.1f} {lam:7.0f} {asym:10.3e} | "
                  + " ".join(f"{per[n]:10.3e}" for n in (0, 1, 2, 3)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
