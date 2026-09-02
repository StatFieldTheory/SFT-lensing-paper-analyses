"""``B_model`` dispatcher: the callable the kappa3 vertex rebuild consumes.

Public entry point::

    from b_model import make_b_delta_fn
    b_delta_fn = make_b_delta_fn("bihalofit")
    B = b_delta_fn(k1, k2, k3, z)     # (Mpc/h)^6

Contract (fixed by ``compute_kappa3_sigma3_high``; see
``notes/input_ready_b_to_zeta_spec.md`` section 2):

* ``k1, k2, k3`` are float64 arrays of a common shape, in h/Mpc, forming
  closed triangles; ``z`` is a scalar (one shell per call)
* the return has that same shape, in (Mpc/h)^6
* the value is the MATTER bispectrum: the vertex applies the per-leg
  Poisson factor A(a) itself, and the caller must bypass the hardwired
  ``growth**4`` because the z dependence now lives inside this function
* the function is symmetric under permutations of its three k arguments
  and never mutates its inputs

Models
------
``tree``            tree-level SPT with D(z)^4 growth (reproduces the
                    deployed table when wired into the build)
``bihalofit``       BiHalofit gravity-only, blended to tree below the
                    low-k numerical seam
``bihalofit_baryon``BiHalofit times the Takahashi baryon ratio R_b
"""

from __future__ import annotations

from typing import Callable, Literal

import numpy as np
from numpy.typing import NDArray

from .cosmology import FIDUCIAL, PK_TABLE_3PT, Fiducial, growth, load_pk_lin
from .halofit_backend import b_bihalofit, baryon_ratio_bihalofit
from .tree import b_tree_z0

Model = Literal["tree", "bihalofit", "bihalofit_baryon"]

#: Low-k seam of the tree/BiHalofit blend, in h/Mpc, applied to the LARGEST
#: leg of each triangle: tree below, BiHalofit above.
#:
#: BiHalofit does NOT reduce to tree level in the linear limit: measured
#: against the exact tree assembly it sits 2-4% high for every k <= 1e-3
#: (1.024 at k = 1e-4, 1.043 at k = 1e-3, equilateral z = 0), a known
#: limitation of the fitting form, well inside its quoted ~10% accuracy but
#: plainly wrong where tree is exact. The seam is therefore placed at the
#: crossover where the two models agree (ratio 1.001 at k = 6.3e-3, z = 0),
#: not deeper in the linear regime where the mismatch is largest.
#:
#: In production this blend never activates: the HIGH branch only keeps
#: triangles whose largest leg exceeds ell = 60, i.e. k_max > 60 / chi_h >
#: 0.038 h/Mpc even at the farthest shell, so every production cell has
#: weight 1. The blend guards out-of-domain and diagnostic calls only.
#: The residual model mismatch across the ell_cut = 60 dispatch boundary
#: (tree in LOW, BiHalofit in HIGH, differing by 1-3% at k ~ 0.04 h/Mpc)
#: is a separate, documented systematic measured by the seam boundary test.
K_BLEND_LO = 5.0e-3
K_BLEND_HI = 3.0e-2

#: BiHalofit calibration ceiling (Takahashi et al. 2020 quote ~10% to
#: k = 3 h/Mpc and ~15% to k = 10 h/Mpc). Exceeding it is allowed but
#: recorded, since the vertex reaches k ~ 3.8 h/Mpc at the nearest shell.
K_CALIBRATED_MAX = 10.0


def _smoothstep(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """C^1 Hermite ramp on [0, 1]."""
    t = np.clip(x, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _blend_weight(k_max: NDArray[np.float64]) -> NDArray[np.float64]:
    """Weight of the nonlinear model: 0 below the seam, 1 above it."""
    lo, hi = np.log(K_BLEND_LO), np.log(K_BLEND_HI)
    with np.errstate(divide="ignore"):
        log_k = np.log(np.maximum(k_max, 1e-300))
    return _smoothstep((log_k - lo) / (hi - lo))


def make_b_delta_fn(
    model: Model = "bihalofit",
    cosmo: Fiducial = FIDUCIAL,
    pk_table: str = PK_TABLE_3PT,
    strict: bool = True,
    growth_fn: Callable[[float], float] | None = None,
) -> Callable[..., NDArray[np.float64]]:
    """Build the ``b_delta_fn(k1, k2, k3, z)`` callable for one model.

    Parameters
    ----------
    model
        Which bispectrum model to serve (see module docstring).
    cosmo
        Fiducial cosmology; must match the P(k) table's normalisation.
    pk_table
        CAMB table name; defaults to the one the 3-point build reads.
    strict
        Raise on non-finite output instead of returning it. Keep this on in
        production: a silent NaN would propagate into the zeta table.
    growth_fn
        Optional ``D(z)`` override for the tree branch, normalised D(0) = 1.
        The default (``cosmology.growth``) solves the same ODE canoes does
        but interpolates it on its own grid, which leaves a ~3e-7 relative
        difference against the vertex's ``lightcone.M[0] / (1 + z)``. Pass
        the lightcone's own growth to make the tree branch reproduce the
        internal assembly to machine precision, as the regression gate does.
    """
    if model not in ("tree", "bihalofit", "bihalofit_baryon"):
        raise ValueError(f"unknown B_model {model!r}")
    pk_lin = load_pk_lin(pk_table, n_s=cosmo.n_s)
    growth_of_z = growth_fn if growth_fn is not None else (lambda z: growth(z, cosmo))

    def b_delta_fn(
        k1: NDArray[np.float64],
        k2: NDArray[np.float64],
        k3: NDArray[np.float64],
        z: float,
    ) -> NDArray[np.float64]:
        k1 = np.asarray(k1, dtype=np.float64)
        k2 = np.asarray(k2, dtype=np.float64)
        k3 = np.asarray(k3, dtype=np.float64)
        if np.any(k1 <= 0) or np.any(k2 <= 0) or np.any(k3 <= 0):
            raise ValueError("all leg moduli must be positive")
        d4 = float(growth_of_z(z)) ** 4
        tree = d4 * b_tree_z0(k1, k2, k3, pk_lin)
        if model == "tree":
            out = tree
        else:
            nonlinear = b_bihalofit(k1, k2, k3, z, cosmo=cosmo, pk_table=pk_table)
            if model == "bihalofit_baryon":
                nonlinear = nonlinear * baryon_ratio_bihalofit(k1, k2, k3, z)
            weight = _blend_weight(np.maximum(np.maximum(k1, k2), k3))
            out = weight * nonlinear + (1.0 - weight) * tree
        if strict and not np.all(np.isfinite(out)):
            n_bad = int(np.sum(~np.isfinite(out)))
            raise FloatingPointError(
                f"B_model({model!r}) produced {n_bad} non-finite values at z={z}"
            )
        return out

    b_delta_fn.model = model  # type: ignore[attr-defined]
    b_delta_fn.cosmo = cosmo  # type: ignore[attr-defined]
    b_delta_fn.pk_table = pk_table  # type: ignore[attr-defined]
    return b_delta_fn
