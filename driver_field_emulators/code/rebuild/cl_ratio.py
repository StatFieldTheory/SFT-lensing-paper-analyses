"""FK relative to Order-0 in harmonic space, deployed against corrected.

The draft's low-multipole claim is a statement about C_ell, so it has to be
checked in C_ell. This uses the paper's OWN curved-sky transform, imported
from the figure generator rather than reimplemented, so the comparison
inherits its Wigner-d kernels, apodisation and DC subtraction exactly.

Usage::

    python cl_ratio.py --deployed <xi.npz> --corrected <xi.npz> --o0 <xi.npz>
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
_A3 = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "analysis3"

#: observable -> (component combos, Wigner-d indices (m, n))
SPECS = {
    "kappa-kappa": ([((0, 0), +1.0)], (0, 0)),
    "xi_plus": ([((1, 1), +1.0), ((2, 2), +1.0)], (2, 2)),
    "xi_minus": ([((1, 1), +1.0), ((2, 2), -1.0)], (2, -2)),
    "kappa_gamma_t": ([((0, 1), -1.0)], (2, 0)),
}


def load_generator():
    if str(_A3) not in sys.path:
        sys.path.insert(0, str(_A3))
    spec = importlib.util.spec_from_file_location(
        "cl_gen", _A3 / "plot_analysis3_cl_decomposition.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--deployed", type=Path, required=True)
    ap.add_argument("--corrected", type=Path, required=True)
    ap.add_argument("--o0", type=Path, required=True)
    ap.add_argument("--observable", default="kappa-kappa", choices=sorted(SPECS))
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    gen = load_generator()
    combos, (m, n) = SPECS[args.observable]

    gamma, o0_groups = gen.load_sweep_order(args.o0, 0)
    _, dep_groups = gen.load_sweep_order(args.deployed, 2)
    _, cor_groups = gen.load_sweep_order(args.corrected, 2)
    o0 = gen._combine(o0_groups, combos)
    dep = gen._combine(dep_groups, combos)
    cor = gen._combine(cor_groups, combos)

    setup = gen.build_curved_matrix(gamma, gen.ELL, m, n)
    cl_o0 = gen.forward_curved(o0, setup)
    cl_dep = gen.forward_curved(dep, setup)
    cl_cor = gen.forward_curved(cor, setup)

    print(f"\n{args.observable}: FK relative to Order-0, curved-sky transform")
    print(f"{'ell':>6} {'C_O0':>13} {'FK dep':>13} {'FK cor':>13} "
          f"{'dep/O0':>10} {'cor/O0':>10}")
    for i, ell in enumerate(gen.ELL):
        base = cl_o0[i]
        if base == 0:
            continue
        print(f"{int(ell):6d} {base:+13.4e} {cl_dep[i]:+13.4e} {cl_cor[i]:+13.4e} "
              f"{cl_dep[i] / base:+10.4f} {cl_cor[i] / base:+10.4f}")

    if args.out:
        np.savez(args.out, ell=gen.ELL, cl_o0=cl_o0, cl_fk_deployed=cl_dep,
                 cl_fk_corrected=cl_cor, observable=args.observable)
        print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
