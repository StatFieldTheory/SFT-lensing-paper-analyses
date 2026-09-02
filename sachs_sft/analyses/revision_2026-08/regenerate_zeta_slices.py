"""Regenerate the zeta-slices figure at the converged multipole cutoff.

Every other number and figure in the revision uses ell_max = 15360; the
zeta-slices figure alone is still drawn from an ell <= 1000 evaluation
(``ell_band_decomp.py`` ships _HIGH_WINDOWS ending at 1000). That matters
because the cutoff changes the figure, not just its normalisation: at shell 8
the normalised zeta_TTT at gamma = 5' moves from -0.87 to -0.27 and the peak
|zeta| grows by a factor 10.6 between the two cutoffs, so the printed panels
show a cumulant falling about three times more slowly in angle than the one
the paper's numbers come from, with annotated peak values an order of
magnitude low.

The exact Wigner-3j LOW ladder stays at ell <= 60 and is untouched; only
flat-sky Born-Limber HIGH windows are appended, which the smoke run times at
about 1.5 s per window per three gammas.

Neither production script is edited: both are imported and their module-level
constants rebound, the same pattern as the other regeneration drivers here.

Usage::

    <sft-wick python> regenerate_zeta_slices.py            # build the npz
    <PyCCL python>    regenerate_zeta_slices.py --plot-only # draw from it
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[4]
_ETL = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables"
        / "kappa3_vertex" / "equal_time_limber")
_OUT_NPZ = _ETL / "outputs" / "zeta_bands_cut15360.npz"
_FIGDIR = _HERE / "outputs"

# LOW stays exact Wigner-3j at ell <= 60; HIGH is flat-sky Born-Limber, so the
# extension to the converged cutoff is four more disjoint windows.
HIGH_WINDOWS = [(60, 125), (125, 250), (250, 500), (500, 1000),
                (1000, 2000), (2000, 4000), (4000, 8000), (8000, 15360)]


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    sys.path.insert(0, str(path.parent))
    spec.loader.exec_module(module)
    return module


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-gamma", type=int, default=28)
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--gamma-max", type=float, default=None,
                    help="clip the angular axis. At ell_max = 15360 the single "
                         "wide Limber window and the sum of the eight disjoint "
                         "windows agree to 3-6%% out to about 85 arcmin and "
                         "diverge by factors of several beyond 117, where |zeta| "
                         "is already below 2%% of its peak: past that the panels "
                         "would show numerically unresolved cancellation as "
                         "banded structure.")
    a = ap.parse_args()

    _FIGDIR.mkdir(parents=True, exist_ok=True)

    if not a.plot_only:
        band = _load(_ETL / "ell_band_decomp.py", "ell_band_decomp")
        band._HIGH_WINDOWS = HIGH_WINDOWS
        print(f"[zeta] HIGH windows -> {HIGH_WINDOWS[-1]} (converged cutoff)",
              flush=True)
        sys.argv = ["ell_band_decomp.py", "--n-gamma", str(a.n_gamma),
                    "--out", str(_OUT_NPZ)]
        rc = band.main()
        if rc:
            return rc
        print(f"[zeta] bands -> {_OUT_NPZ}", flush=True)
        return 0

    npz = _OUT_NPZ
    if a.gamma_max is not None:
        import numpy as np
        d = dict(np.load(_OUT_NPZ, allow_pickle=True))
        g = np.asarray(d["gamma_arcmin"], float)
        keep = g <= a.gamma_max
        for k, v in list(d.items()):
            v = np.asarray(v)
            if v.ndim >= 1 and v.shape[0] == g.size:
                d[k] = v[keep]
        npz = _OUT_NPZ.with_name(
            _OUT_NPZ.stem + f"_gmax{int(a.gamma_max)}.npz")
        np.savez(npz, **d)
        print(f"[zeta] clipped to gamma <= {a.gamma_max}' "
              f"({keep.sum()} of {g.size} points) -> {npz.name}", flush=True)

    plot = _load(_ETL / "plot_zeta_figures.py", "plot_zeta_figures")
    plot._NPZ = npz
    plot._FIGDIR = _FIGDIR
    print(f"[zeta] plotting {npz.name} -> {_FIGDIR}", flush=True)
    return plot.main() or 0


if __name__ == "__main__":
    raise SystemExit(main())
