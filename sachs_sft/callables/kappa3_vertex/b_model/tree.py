"""Tree-level SPT matter bispectrum, assembled exactly as canoes does.

This is the reference branch of ``B_model``: it must reproduce the internal
tree assembly of ``compute_kappa3_sigma3_high`` term by term, so that the
rebuilt pipeline with ``b_delta_fn = B_model(model="tree")`` returns the
deployed table (and hence the +4.29e-6 FK kk baseline) unchanged.

Conventions
-----------
* wavenumbers in h/Mpc, P in (Mpc/h)^3, B in (Mpc/h)^6, matter (no Poisson
  factor: the vertex applies A(a) per leg internally)
* closed triangles, k1 + k2 + k3 = 0 as vectors, so the cosine between the
  legs (i, j) follows from the third modulus:
  ``mu_ij = (k_l^2 - k_i^2 - k_j^2) / (2 k_i k_j)``
  This is algebraically identical to the ``-(u + v cos phi) / w`` form the
  canoes quadrature uses for its (u, v, phi) parameterisation.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

ArrayLike = NDArray[np.float64]


def f2_kernel(k_a: ArrayLike, k_b: ArrayLike, mu: ArrayLike) -> ArrayLike:
    """Symmetrised second-order SPT mode-coupling kernel F2^(s)."""
    return 5.0 / 7.0 + 0.5 * mu * (k_a / k_b + k_b / k_a) + (2.0 / 7.0) * mu**2


def closure_cosines(
    k1: ArrayLike, k2: ArrayLike, k3: ArrayLike
) -> tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Return ``(mu_12, mu_13, mu_23)`` for a closed triangle."""
    mu_12 = (k3**2 - k1**2 - k2**2) / (2.0 * k1 * k2)
    mu_13 = (k2**2 - k1**2 - k3**2) / (2.0 * k1 * k3)
    mu_23 = (k1**2 - k2**2 - k3**2) / (2.0 * k2 * k3)
    return mu_12, mu_13, mu_23


def b_tree_z0(
    k1: ArrayLike, k2: ArrayLike, k3: ArrayLike, pk_lin
) -> ArrayLike:
    """Tree-level matter bispectrum at z = 0.

    Parameters
    ----------
    k1, k2, k3
        Leg moduli in h/Mpc, broadcastable, forming closed triangles.
    pk_lin
        Callable returning the linear matter P(k, z=0) in (Mpc/h)^3.
    """
    p1, p2, p3 = pk_lin(k1), pk_lin(k2), pk_lin(k3)
    mu_12, mu_13, mu_23 = closure_cosines(k1, k2, k3)
    f12 = f2_kernel(k1, k2, mu_12)
    f13 = f2_kernel(k1, k3, mu_13)
    f23 = f2_kernel(k2, k3, mu_23)
    return 2.0 * (f12 * p1 * p2 + f13 * p1 * p3 + f23 * p2 * p3)


def b_tree(
    k1: ArrayLike, k2: ArrayLike, k3: ArrayLike, z: float, pk_lin, growth_fn
) -> ArrayLike:
    """Tree-level matter bispectrum at redshift ``z``.

    Uses the separable growth scaling ``B(z) = D(z)^4 B(0)`` that the paper
    adopts (Eq. appendix bdelta growth): three linear legs carry one D each
    and the second-order leg carries D^2, so the product is D^4 on the
    equal-time diagonal. Exact only within the EdS-kernel approximation.
    """
    return float(growth_fn(z)) ** 4 * b_tree_z0(k1, k2, k3, pk_lin)
