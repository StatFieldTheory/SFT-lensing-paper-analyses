"""ADVERSARIAL PROBE 2 -- where (if anywhere) is the NON-tautological link?

Probe 1 showed the gate-1 identity survives ANY consistent change of the
physics.  So the question the attack really poses is: which comparison in the
chain actually discriminates the paper's FK diagram from a wrong one?

Two candidate references are tested for discriminating power:

  INTERNAL   fk_reference(grid, zeta_injected)  -- the discrete white-noise
             (SPEC section 3) fold built from the SAME grid kernels.
  EXTERNAL   the paper's converged permutation-closed fold,
             products/table_permclosed_cut15360_permfix_xi.npz, produced by the
             sft-wick workflow, which knows nothing about this discretisation.

For each of five stochastic systems (one true, four deliberately wrong but
internally consistent -- all of which PASS gate 1 at 4e-16) we compute the exact
expectation at two sigma_lambda on a fixed sigma/dlam ladder, extrapolate
linearly to sigma_lambda = 0, and report the ratio to BOTH references.

If the internal reference gives 1.000 for the wrong systems too, it is as blind
as gate 1 and the entire non-tautological content sits in the external fold.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import fk_expect_exact as ex                       # noqa: E402
import _bootstrap                                  # noqa: E402
from a1_teeth_and_blindness import remake_grid     # noqa: E402


def analytic_fk(gammas):
    d = np.load(_bootstrap.FOLD, allow_pickle=True)
    m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)), -1, 1)))) * 60)
    g = np.asarray(g); v = np.asarray(d["value"], float)[m]
    o = np.argsort(g)
    return np.interp(gammas, g[o], v[o])


def variants(grid):
    rng = np.random.default_rng(7)
    N = grid.lam.size
    w_rand = grid.dlam * (1.0 + rng.random(N))
    return [
        ("TRUE system (production)", grid),
        ("WRONG: propagator (D_j/D_k)^1", remake_grid(grid, D=grid.D ** 0.5)),
        ("WRONG: propagator (D_j/D_k)^6", remake_grid(grid, D=grid.D ** 3.0)),
        ("WRONG: D = 1 (no background)",
         remake_grid(grid, D=np.ones_like(grid.D))),
        ("WRONG: random observable window", remake_grid(grid, w=w_rand)),
    ]


def main() -> int:
    ladder = [(1000, 8.0), (1999, 4.0)]        # fixed sigma/dlam = 4.19
    gammas = [1.0, 5.0, 17.3]
    res = {}
    for gam in gammas:
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        for N, sig in ladder:
            base = ex.build_grid(n_lambda=N)
            V, A, Q, rho = ex.node_stats(base, cosg, sig, calibrate="smeared")
            M = ex.calib_matrix(base, V, A, rho, sig, "smeared")
            zi = ex.zeta_injected(M, Q)
            for name, g in variants(base):
                T1, T2 = ex.expectation(g, cosg, sig, V=V, A=A, Q=Q, rho=rho)
                exact = ((T1 + T2) + (T1 + T2).T)[0, 3]
                internal = ex.fk_reference(g, zi)[0, 3]
                res[(gam, sig, name)] = (exact, internal)
            print(f"  done gamma={gam} sigma={sig} N={N}", flush=True)

    print("\n" + "=" * 96)
    print("Exact expectation, extrapolated to sigma_lambda -> 0 by 2*r(4) - r(8),")
    print("against the INTERNAL discrete reference and the EXTERNAL paper fold.")
    print("Every system below passes the gate-1 placement identity at ~4e-16.")
    print("=" * 96)
    hdr = f"{'system':<34}"
    for gam in gammas:
        hdr += f" | {'g=' + str(gam):>19}"
    print(hdr)
    print(f"{'':<34}" + " | " + " | ".join(
        [f"{'int':>8} {'EXTERNAL':>10}"] * len(gammas)))
    names = [n for n, _ in variants(ex.build_grid(n_lambda=16))]
    for name in names:
        row = f"{name:<34}"
        for gam in gammas:
            ana = float(analytic_fk([gam])[0])
            e8, i8 = res[(gam, 8.0, name)]
            e4, i4 = res[(gam, 4.0, name)]
            r_int = 2 * (e4 / i4) - (e8 / i8)
            r_ext = 2 * (e4 / ana) - (e8 / ana)
            row += f" | {r_int:>8.4f} {r_ext:>10.4f}"
        print(row)
    print("\n'int' = exact / (that system's OWN discrete white-noise fold)")
    print("'EXTERNAL' = exact / (paper fold, sft-wick, same zeta table)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
