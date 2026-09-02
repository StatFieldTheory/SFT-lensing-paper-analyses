"""Integration tests: B_model wired into the canoes kappa3 HIGH branch.

These are the gates that protect the rebuild against "new plumbing" errors
before any nonlinear number is trusted:

* the tree branch routed through ``b_delta_fn`` must reproduce the internal
  tree assembly (the deployed table's physics) to numerical noise
* the nonlinear branch must change the vertex in the expected direction and
  by a bounded amount
* the metadata must record which bispectrum produced the table

Run from ``sachs_sft/callables/kappa3_vertex``::

    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python -m pytest b_model/tests -q
"""

from __future__ import annotations

import numpy as np
import pytest

from canoes.sachs.kappa3 import (
    compute_kappa3_mod_sigma3_high,
    compute_kappa3_sigma3_high,
)

from b_model import FIDUCIAL, load_pk_lin, make_b_delta_fn

CHANNELS = ("TTT", "TTP", "TPP", "PPP")

#: Tree-vs-internal agreement is limited by the growth-factor interpolation
#: route, not by the bispectrum: canoes interpolates D on its lightcone
#: z-grid while b_model interpolates the growth ODE on its a-grid, giving
#: dD/D ~ 3e-7 and hence dB/B ~ 4 dD/D ~ 1e-6. That is four orders below the
#: 1e-3 boundary-validation gate. Pass ``growth_fn`` to remove it entirely.
TREE_REGRESSION_RTOL = 5.0e-6

#: Collapsed family (1, c, c) plus open triangles, the shapes the FK fold
#: and the general table respectively sample.
TRIPLES = np.array(
    [
        [1.0, 0.999999, 0.999999],  # gamma ~ 5 arcmin, corner of the FK slice
        [1.0, 0.96, 0.96],          # gamma ~ 16 deg
        [1.0, 0.12, 0.12],          # gamma ~ 83 deg, widest production point
        [0.8, 0.7, 0.6],            # open triangle
        [1.0, 1.0, 1.0],            # exact corner
    ]
)

#: Extreme shells of the deployed grid, in lambda_project [Mpc/h].
LAMBDA_SHELLS = np.array([396.63, 1200.0, 2326.61]) * FIDUCIAL.h


def _kwargs():
    return dict(
        pk_matter_today=load_pk_lin(),
        cosmo=FIDUCIAL,
        ell_cut=60,
        ell_high_max=1000,
        n_ell=32,
        n_phi=32,
    )


def test_tree_via_callable_reproduces_internal_assembly():
    """The regression gate: new plumbing must not change the physics."""
    kwargs = _kwargs()
    reference = compute_kappa3_sigma3_high(TRIPLES, LAMBDA_SHELLS, **kwargs)
    routed = compute_kappa3_sigma3_high(
        TRIPLES, LAMBDA_SHELLS, **kwargs, b_delta_fn=make_b_delta_fn("tree")
    )
    for channel in CHANNELS:
        expected = getattr(reference, f"zeta_{channel}")
        got = getattr(routed, f"zeta_{channel}")
        np.testing.assert_allclose(got, expected, rtol=TREE_REGRESSION_RTOL)


def test_tree_via_callable_reproduces_modulus_channels():
    kwargs = _kwargs()
    reference = compute_kappa3_mod_sigma3_high(TRIPLES, LAMBDA_SHELLS, **kwargs)
    routed = compute_kappa3_mod_sigma3_high(
        TRIPLES, LAMBDA_SHELLS, **kwargs, b_delta_fn=make_b_delta_fn("tree")
    )
    for expected, got in zip(reference, routed):
        np.testing.assert_allclose(got, expected, rtol=TREE_REGRESSION_RTOL)


def test_metadata_records_the_bispectrum_source():
    kwargs = _kwargs()
    internal = compute_kappa3_sigma3_high(TRIPLES, LAMBDA_SHELLS, **kwargs)
    assert internal.cosmo_meta["b_delta_source"] == "tree_spt_D4"
    routed = compute_kappa3_sigma3_high(
        TRIPLES, LAMBDA_SHELLS, **kwargs, b_delta_fn=make_b_delta_fn("bihalofit")
    )
    assert routed.cosmo_meta["b_delta_source"] == "b_delta_fn:bihalofit"


def test_nonlinear_vertex_is_enhanced_and_bounded():
    """BiHalofit must amplify the vertex, strongly near and weakly far.

    zeta is a signed band-cancelling integral (the ell bands carry opposite
    signs, see the archived ell-band table), so the nonlinear-to-tree ratio
    is NOT bounded below by 1 and can even change sign when the boosted
    high-ell band overturns the low-ell one. The physical expectations that
    do hold are on magnitudes: strong amplification at the nearest shell,
    where the hard legs are most nonlinear, and near-identity at the
    farthest shell, where growth suppression leaves the vertex tree-like.
    """
    kwargs = _kwargs()
    tree = compute_kappa3_sigma3_high(TRIPLES, LAMBDA_SHELLS, **kwargs)
    nonlinear = compute_kappa3_sigma3_high(
        TRIPLES, LAMBDA_SHELLS, **kwargs, b_delta_fn=make_b_delta_fn("bihalofit")
    )
    assert np.all(np.isfinite(nonlinear.zeta_TTT))

    amplification = np.abs(nonlinear.zeta_TTT) / np.abs(tree.zeta_TTT)
    assert np.all(amplification[:, 0] > 2.0), f"near shell not amplified: {amplification[:, 0]}"
    assert np.all(amplification[:, 0] < 200.0), f"runaway: {amplification[:, 0]}"
    np.testing.assert_allclose(amplification[:, -1], 1.0, atol=0.02)
    assert np.all(amplification[:, 0] > amplification[:, -1])


def test_callable_shape_and_finiteness_are_enforced():
    """A malformed b_delta_fn must fail loud, never silently corrupt a table."""
    kwargs = _kwargs()
    with pytest.raises(ValueError, match="shape"):
        compute_kappa3_sigma3_high(
            TRIPLES, LAMBDA_SHELLS, **kwargs,
            b_delta_fn=lambda k1, k2, k3, z: np.ones(3),
        )
    with pytest.raises(ValueError, match="non-finite"):
        compute_kappa3_sigma3_high(
            TRIPLES, LAMBDA_SHELLS, **kwargs,
            b_delta_fn=lambda k1, k2, k3, z: np.full_like(k1, np.nan),
        )
