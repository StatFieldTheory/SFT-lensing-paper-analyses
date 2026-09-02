"""ADVERSARIAL PROBE 3 -- resolving power of the ONLY non-blind comparison.

Probe 2 showed the external paper fold is the only reference that separates the
true stochastic system from wrong ones.  Caricatures (D = 1) are separated by
factors of 2-3.  The question that matters for the paper's claim is how small a
kernel error survives: graded errors in the propagator exponent, plus the
2-per-mille trapezoid convention.

Each row is again a fully self-consistent system that passes gate 1 at ~4e-16.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import fk_expect_exact as ex                          # noqa: E402
from a1_teeth_and_blindness import remake_grid        # noqa: E402
from a2_external_link import analytic_fk              # noqa: E402


def main() -> int:
    ladder = [(1000, 8.0), (1999, 4.0)]
    gammas = [1.0, 5.0, 17.3]
    powers = [0.90, 0.95, 0.98, 0.99, 1.00, 1.01, 1.02, 1.05, 1.10]
    res = {}
    for gam in gammas:
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        for N, sig in ladder:
            base = ex.build_grid(n_lambda=N)
            V, A, Q, rho = ex.node_stats(base, cosg, sig, calibrate="smeared")
            rows = [(f"prop exponent {2*p:.2f}", remake_grid(base, D=base.D ** p))
                    for p in powers]
            wflat = np.full(N, base.dlam)
            rows.append(("trapezoid end weights dropped",
                         remake_grid(base, w=wflat)))
            for name, g in rows:
                T1, T2 = ex.expectation(g, cosg, sig, V=V, A=A, Q=Q, rho=rho)
                res[(gam, sig, name)] = ((T1 + T2) + (T1 + T2).T)[0, 3]
            print(f"  done gamma={gam} sigma={sig}", flush=True)

    names = [f"prop exponent {2*p:.2f}" for p in powers] \
        + ["trapezoid end weights dropped"]
    print("\n" + "=" * 78)
    print("exact expectation extrapolated to sigma->0, / PAPER FOLD")
    print("(the true system is 'prop exponent 2.00'; all rows pass gate 1)")
    print("=" * 78)
    print(f"{'system':<34} " + " ".join(f"{'g=' + str(g):>11}" for g in gammas))
    for name in names:
        row = f"{name:<34} "
        for gam in gammas:
            ana = float(analytic_fk([gam])[0])
            r = 2 * (res[(gam, 4.0, name)] / ana) - (res[(gam, 8.0, name)] / ana)
            row += f"{r:>11.4f} "
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
