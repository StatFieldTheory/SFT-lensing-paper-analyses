"""BiHalofit backend, reusing the fastnc implementation.

We call the BiHalofit fitting formula of Takahashi et al. 2020
(arXiv:1911.07886) through the ``Halofit`` class of fastnc
(https://github.com/git-sunao/fastnc, MIT, Sunao Sugiyama), the same code
already used for the appendix shear-3PCF cross-check of this project.

fastnc is installed with ``--no-deps`` (its full dependency set pulls in
astropy, pandas and mpi4py, none of which the BiHalofit path needs), so
``import fastnc`` fails at package import time. ``halofit.py`` itself only
needs numpy and scipy, so we load that single module directly from the
installed package directory and bypass ``fastnc/__init__.py``.

Units: k in h/Mpc, P in (Mpc/h)^3, B in (Mpc/h)^6, matter, equal time.
"""

from __future__ import annotations

import importlib.util
import sysconfig
from functools import lru_cache
from pathlib import Path
from types import ModuleType

import numpy as np
from numpy.typing import NDArray

from .cosmology import FIDUCIAL, Fiducial, growth, load_pk_table

#: Elements per BiHalofit call. The coefficient table is a ~20-field
#: structured array, so a full (96, 96, 64) grid in one call would allocate
#: over 100 MB; chunking keeps the footprint flat.
CHUNK = 200_000

#: z grid on which fastnc builds its r_sigma / neff / sigma8(z) splines.
#: Must cover every shell of the kappa3 table (z_max = 5.7).
Z_SPLINE_MAX = 6.5


def _module_search_paths() -> list[Path]:
    candidates = [Path(sysconfig.get_paths()["purelib"]) / "fastnc" / "halofit.py"]
    try:  # editable / source checkouts
        import fastnc  # noqa: F401  (import may fail on heavy deps)
    except Exception:
        pass
    return candidates


@lru_cache(maxsize=1)
def load_fastnc_halofit() -> ModuleType:
    """Import fastnc's ``halofit`` module without executing the package init."""
    for path in _module_search_paths():
        if path.exists():
            spec = importlib.util.spec_from_file_location("_fastnc_halofit", path)
            if spec is None or spec.loader is None:  # pragma: no cover
                continue
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise ImportError(
        "fastnc halofit.py not found. Install it with:\n"
        "    python -m pip install --no-deps fastnc"
    )


@lru_cache(maxsize=4)
def _halofit_instance(cosmo: Fiducial, pk_table: str) -> object:
    """Build (and cache) a fastnc ``Halofit`` for one cosmology and P(k)."""
    module = load_fastnc_halofit()
    k_h, p_m = load_pk_table(pk_table)
    z_grid = np.linspace(0.0, Z_SPLINE_MAX, 66)
    cosmo_dict = {
        "Om0": cosmo.Omega_m,
        "Ode0": cosmo.Omega_de,
        "ns": cosmo.n_s,
        "w0": -1.0,
        "wa": 0.0,
        "fnu0": 0.0,
        "sigma8": cosmo.sigma8,
    }
    # fastnc 1.2.3's Halofit.__init__ calls set_pklin before set_cosmology,
    # and set_pklin's sigma8 renormalisation reads self.cosmo, so the
    # all-arguments constructor raises AttributeError. Build empty (the
    # normalisation short-circuits on k is None) and set in dependency order.
    # Halofit renormalises pklin in place, so pass copies and never let the
    # caller's table be mutated.
    halofit = module.Halofit()
    halofit.set_cosmology(cosmo_dict)
    halofit.set_lgr(z_grid, np.asarray(growth(z_grid, cosmo), dtype=np.float64))
    halofit.set_pklin(k_h.copy(), p_m.copy())
    return halofit


def b_bihalofit(
    k1: NDArray[np.float64],
    k2: NDArray[np.float64],
    k3: NDArray[np.float64],
    z: float,
    cosmo: Fiducial = FIDUCIAL,
    pk_table: str = "PCAMBz0.txt",
) -> NDArray[np.float64]:
    """Nonlinear matter bispectrum from BiHalofit, in (Mpc/h)^6.

    The triangles are assumed closed (``all_physical=True``): the kappa3
    quadrature builds the third leg as ``|l2 + l3|``, so the triangle
    inequality holds by construction and fastnc's own float-tolerance
    rejection would only mask degenerate-but-valid cells as NaN.
    """
    halofit = _halofit_instance(cosmo, pk_table)
    shape = np.broadcast(k1, k2, k3).shape
    flat = [np.broadcast_to(np.asarray(k, dtype=np.float64), shape).ravel()
            for k in (k1, k2, k3)]
    out = np.empty(flat[0].size, dtype=np.float64)
    for start in range(0, flat[0].size, CHUNK):
        stop = min(start + CHUNK, flat[0].size)
        sl = slice(start, stop)
        z_chunk = np.full(stop - start, float(z), dtype=np.float64)
        out[sl] = halofit.get_bihalofit(
            flat[0][sl].copy(), flat[1][sl].copy(), flat[2][sl].copy(),
            z_chunk, all_physical=True, which=["Bh3", "Bh1"],
        )
    return out.reshape(shape)


def baryon_ratio_bihalofit(
    k1: NDArray[np.float64],
    k2: NDArray[np.float64],
    k3: NDArray[np.float64],
    z: float,
) -> NDArray[np.float64]:
    """Takahashi et al. baryon-correction ratio R_b for the bispectrum.

    The IllustrisTNG-calibrated companion fit of BiHalofit. It is a cheap
    alternative to the baccoemu boost, defined only for z < 5 (fastnc
    returns 1 above that).

    Availability: released fastnc 1.2.3 on PyPI ships halofit.py without
    this helper; it exists in the upstream working tree. If it is missing,
    use the baccoemu boost instead (the planned baryon route) rather than
    silently returning unity.
    """
    module = load_fastnc_halofit()
    if not hasattr(module, "get_Rb_bihalofit"):
        raise NotImplementedError(
            "the installed fastnc halofit.py has no get_Rb_bihalofit; use the "
            "baccoemu boost for baryons, or install a fastnc build that "
            "provides it"
        )
    shape = np.broadcast(k1, k2, k3).shape
    flat = [np.broadcast_to(np.asarray(k, dtype=np.float64), shape).ravel()
            for k in (k1, k2, k3)]
    out = np.empty(flat[0].size, dtype=np.float64)
    for start in range(0, flat[0].size, CHUNK):
        stop = min(start + CHUNK, flat[0].size)
        sl = slice(start, stop)
        z_chunk = np.full(stop - start, float(z), dtype=np.float64)
        out[sl] = module.get_Rb_bihalofit(
            flat[0][sl].copy(), flat[1][sl].copy(), flat[2][sl].copy(), z_chunk
        )
    return out.reshape(shape)
