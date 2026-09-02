"""B_model: matter-bispectrum input for the SFT lensing kappa3 vertex.

Serves ``b_delta_fn(k1, k2, k3, z)`` in the exact units and conventions the
canoes ``equal_time_limber`` build consumes: k in h/Mpc, B in (Mpc/h)^6,
matter (no Poisson factor), equal time, one scalar redshift per call.

See ``../../notes/input_ready_b_to_zeta_spec.md`` for the full chain from
this callable to the configuration-space zeta table and on to sft-wick.
"""

from .cosmology import (
    FIDUCIAL,
    PK_TABLE_2PT,
    PK_TABLE_3PT,
    SIGMA8_2PT,
    SIGMA8_3PT,
    Fiducial,
    growth,
    load_pk_lin,
    load_pk_table,
    sigma8_of_table,
)
from .core import K_BLEND_HI, K_BLEND_LO, make_b_delta_fn
from .halofit_backend import b_bihalofit, baryon_ratio_bihalofit
from .tree import b_tree, b_tree_z0, closure_cosines, f2_kernel

__all__ = [
    "FIDUCIAL",
    "Fiducial",
    "K_BLEND_HI",
    "K_BLEND_LO",
    "PK_TABLE_2PT",
    "PK_TABLE_3PT",
    "SIGMA8_2PT",
    "SIGMA8_3PT",
    "b_bihalofit",
    "b_tree",
    "b_tree_z0",
    "baryon_ratio_bihalofit",
    "closure_cosines",
    "f2_kernel",
    "growth",
    "load_pk_lin",
    "load_pk_table",
    "make_b_delta_fn",
    "sigma8_of_table",
]
