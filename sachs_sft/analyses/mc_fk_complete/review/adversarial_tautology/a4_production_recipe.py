"""ADVERSARIAL PROBE 4 -- discrimination under the recipe the PAPER actually uses.

Probes 2-3 used a fixed sigma/dlam ladder.  The published number
(NOTES gate 2/4, analyse_final.py) uses FIXED n_lambda = 1000 with
sigma_lambda = 8 and 4 and the linear extrapolation 2*r(4) - r(8).  This probe
re-runs the discrimination test in exactly that recipe, so the numbers quoted
are the resolving power of the paper's own check.

Each row is a self-consistent stochastic system that passes the gate-1
placement identity at ~4e-16 (probe 1) and converges to its OWN SPEC-section-3
white-noise reference at 1.000 (probe 2).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import fk_expect_exact as ex                        # noqa: E402
from a1_teeth_and_blindness import remake_grid      # noqa: E402
from a2_external_link import analytic_fk            # noqa: E402


def build_rows(base):
    rng = np.random.default_rng(7)
    N = base.lam.size
    rows = [("TRUE system (production)", base)]
    for p in (0.5, 0.98, 0.99, 1.01, 1.02, 1.05, 1.10, 3.0):
        rows.append((f"prop exponent {2*p:.2f}  (true = 2.00)",
                     remake_grid(base, D=base.D ** p)))
    rows.append(("D = 1 (background off)",
                 remake_grid(base, D=np.ones_like(base.D))))
    rows.append(("random observable window",
                 remake_grid(base, w=base.dlam * (1.0 + rng.random(N)))))
    rows.append(("trapezoid end weights dropped",
                 remake_grid(base, w=np.full(N, base.dlam))))
    return rows


def main() -> int:
    N = 1000
    gammas = [1.0, 5.0, 17.3]
    res, refs = {}, {}
    for gam in gammas:
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        for sig in (8.0, 4.0):
            base = ex.build_grid(n_lambda=N)
            V, A, Q, rho = ex.node_stats(base, cosg, sig, calibrate="smeared")
            M = ex.calib_matrix(base, V, A, rho, sig, "smeared")
            zi = ex.zeta_injected(M, Q)
            for name, g in build_rows(base):
                T1, T2 = ex.expectation(g, cosg, sig, V=V, A=A, Q=Q, rho=rho)
                res[(gam, sig, name)] = ((T1 + T2) + (T1 + T2).T)[0, 3]
                refs[(gam, sig, name)] = ex.fk_reference(g, zi)[0, 3]
            print(f"  done gamma={gam} sigma={sig}", flush=True)

    names = [n for n, _ in build_rows(ex.build_grid(n_lambda=16))]
    print("\n" + "=" * 92)
    print("PRODUCTION RECIPE: N=1000 fixed, sigma=8 and 4, extrapolated 2*r(4)-r(8)")
    print("  'own ref'  = / that system's own SPEC-3 discrete white-noise fold")
    print("  'PAPER'    = / the sft-wick converged fold (the only external number)")
    print("=" * 92)
    head = f"{'system':<36}"
    for gam in gammas:
        head += f" |{'g=' + str(gam) + ' own':>12}{'PAPER':>9}"
    print(head)
    for name in names:
        row = f"{name:<36}"
        for gam in gammas:
            ana = float(analytic_fk([gam])[0])
            e8, e4 = res[(gam, 8.0, name)], res[(gam, 4.0, name)]
            i8, i4 = refs[(gam, 8.0, name)], refs[(gam, 4.0, name)]
            row += f" |{2*(e4/i4) - (e8/i8):>12.4f}{2*(e4/ana) - (e8/ana):>9.4f}"
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
