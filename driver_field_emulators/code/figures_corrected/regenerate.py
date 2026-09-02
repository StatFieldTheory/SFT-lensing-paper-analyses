"""Regenerate the vertex-dependent paper figures against a corrected FK result.

The paper's generators pin their inputs and outputs with module-level
constants and write into the repository. Rather than copy and edit each
one, which would fork six scripts and let them drift, this driver imports
each generator as a module, rebinds the two constants that matter, and
calls its `main()`. The plotting code, styling and observable definitions
are therefore the paper's own, unmodified, and the only difference is which
FK file is read and where the figure is written.

Nothing here writes into `figures/`, `sections/` or any deployed run
directory. Output goes to `driver_field_emulators/figures/corrected/`.

Usage::

    python regenerate.py --fk-npz ../../products/<corrected>_xi.npz \
                         --figure nlo cl eb
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[2]
_ANALYSIS3 = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft"
              / "analyses" / "analysis3")
_MC = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft"
       / "analyses" / "mc_sachs_2pt")
_OUT = _REPO / "driver_field_emulators" / "figures" / "corrected"

#: figure key -> (directory, module file, output-stem attribute)
GENERATORS = {
    "nlo": (_ANALYSIS3, "plot_analysis3_nlo_decomposition.py", "OUT_STEM"),
    "cl": (_ANALYSIS3, "plot_analysis3_cl_decomposition.py", "OUT_STEM"),
    "eb": (_ANALYSIS3, "plot_cl_EB_polarization.py", "OUT_STEM"),
    "mc": (_MC, "fig_xi_channels.py", None),
}


def _load(directory: Path, filename: str):
    """Import a generator from its own directory so its sibling imports work."""
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
    spec = importlib.util.spec_from_file_location(
        Path(filename).stem, directory / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _rebind_fk(module, fk_npz: Path) -> int:
    """Point every FK_NPZ this module can see at `fk_npz`. Returns the count.

    `plot_cl_EB_polarization` reads the constant off the cl module rather
    than defining its own, so the rebinding has to follow module references,
    not just the module's own namespace.
    """
    count = 0
    if hasattr(module, "FK_NPZ"):
        module.FK_NPZ = fk_npz
        count += 1
    for value in vars(module).values():
        if hasattr(value, "FK_NPZ") and getattr(value, "__name__", "").startswith("plot_"):
            value.FK_NPZ = fk_npz
            count += 1
    return count


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fk-npz", type=Path, required=True)
    ap.add_argument("--figure", nargs="+", default=["nlo", "cl", "eb"],
                    choices=sorted(GENERATORS))
    ap.add_argument("--suffix", default="corrected",
                    help="appended to each output stem")
    ap.add_argument("--out-dir", type=Path, default=_OUT)
    args = ap.parse_args()

    fk = args.fk_npz.resolve()
    if not fk.exists():
        raise SystemExit(f"FK npz not found: {fk}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for key in args.figure:
        directory, filename, stem_attr = GENERATORS[key]
        cwd = Path.cwd()
        try:
            module = _load(directory, filename)
            rebound = _rebind_fk(module, fk)
            if rebound == 0:
                raise SystemExit(f"{filename}: no FK_NPZ constant to rebind")
            if stem_attr is not None:
                old = getattr(module, stem_attr)
                setattr(module, stem_attr,
                        args.out_dir / f"{Path(old).name}_{args.suffix}")
            print(f"[fig] {key}: {filename}, {rebound} FK reference(s) rebound")
            code = module.main()
            if code:
                raise SystemExit(f"{filename} returned {code}")
        finally:
            import os
            os.chdir(cwd)
    print(f"[fig] wrote into {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
