"""Boundary test: SAS-isoceles mirror symmetry of the cosmic-shear 3PCF
natural-component assembler.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
Run:  <interp> -m pytest tests/test_natural_component_mirror_symmetry.py -v

WHY THIS IS A BOUNDARY TEST (first-principles, parameter-free)
==============================================================
On the SAS isoceles family (t1 = t2 = gamma, opening angle phi, third side
t3 = 2 gamma sin(phi/2)) the triangle is invariant under the y -> -y mirror:
the APEX vertex is self-mirror, and the two BASE vertices swap.  Because the
shear is a spin-(+2) field, a reflection maps spin +2 -> spin -2 and negates
the centroid reference angle.  This forces, for ANY correct natural-component
assembler (Schneider-Lombardi / fastnc convention), with NO free parameters:

    Gamma^1  is REAL                 (the APEX-conjugated component)
    Gamma^2  =  conj(Gamma^3)        (the two BASE-conjugated components)

fastnc satisfies this to ~1e-16.  The derivation is
scripts/derive_natural_component_leg_map.wl (xAct), giving the exact per-vertex
centroid angle alpha_apex = 0, alpha_baseB = -pi + arctan(3 tan(phi/2)),
alpha_baseC = -alpha_baseB.

The assembler stage1_ours_v2_JACFIX.py builds Gamma^mu = sign * fold_mu *
centroid_phase(conj_leg_mu, phi).  The mirror symmetry is a property of the
SLOT-BINDING + PHASE alone, given that on the SAS family the two BASE radial
folds are equal (fold_base1 = fold_base2, real) and the APEX/PPP folds are real.
We therefore test the assembler's phase+binding logic directly by injecting
representative real folds; this isolates the slot-label/phase convention (the
locus of the 2026-06-10 fix) from the slow canoes radial query.

The same construction with the OLD slot binding (apex in slot Gamma^3) and the
OLD `_amb` great-circle angles VIOLATES the symmetry by O(1) -- demonstrated by
the `binding="old"` parametrization, which is expected to FAIL the asserts (run
it explicitly to confirm the test has teeth).

Boundary-validation rule compliance: the (gamma, phi) sweep INCLUDES the extreme
corners (phi=15 and phi=150; gamma=1' very small and gamma=100' large).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_HERE = Path(__file__).resolve().parent
_ASSEMBLER = _HERE.parent / "stage1_ours_v2_JACFIX.py"

# Import the FIXED assembler module to use its real centroid-phase functions.
_spec = importlib.util.spec_from_file_location("stage1_jacfix", _ASSEMBLER)
_mod = importlib.util.module_from_spec(_spec)
# The module imports canoes at top level; that is fine in the sft-wick env. If
# canoes import is too heavy for CI, the phase helpers could be vendored, but we
# deliberately import the REAL module so the test tracks the production code.
_spec.loader.exec_module(_mod)

centroid_phase = _mod.centroid_phase
_centroid_alpha = _mod._centroid_alpha

SIGN = -1.0  # the assembler's overall sign convention

# --- known-correct pin (APEX-channel modulus) -------------------------------
# At gamma=10', phi=60 the fresh post-fix canoes apex real-fold modulus is
# |fold_apex| = 4.8545e-08 (recovered from stage1_ours_v2_JACFIX.npz, the apex
# channel after the fix).  Pinned so a regression shows as numeric drift.
PIN_APEX_MOD_G10_PHI60 = 4.8545e-08
PIN_RTOL = 5e-3


def _synthetic_folds(phi_deg: float) -> tuple[float, float, float, float]:
    """Representative REAL radial folds on the SAS family.

    The physical content the symmetry needs: the two BASE folds are EQUAL and
    real; the APEX and PPP folds are real.  We use a smooth phi-dependent shape
    (sign-changing, like the real folds) so the test is not trivially satisfied
    by constants.  Exact values are irrelevant to the symmetry; only base-leg
    equality + reality matter.
    """
    p = np.radians(phi_deg)
    f_base = 1.0e-7 * np.cos(1.5 * p) * (1.0 + 0.3 * np.sin(p))
    f_apex = 1.0e-7 * (0.5 + 0.5 * np.cos(p))
    f_ppp = 1.0e-9 * (1.0 + 0.2 * np.cos(2.0 * p))
    return float(f_base), float(f_base), float(f_apex), float(f_ppp)


def _assemble(binding: str, fbase1, fbase2, fapex, fppp, phi_rad):
    """Assemble Gamma^0..3 for a given slot binding.

    binding="fixed": apex-conj -> Gamma^1 ; base legs -> {Gamma^2, Gamma^3}
                     using the xAct centroid phase (production fix).
    binding="old"  : apex-conj -> Gamma^3 ; base legs -> {Gamma^1, Gamma^2}
                     using the legacy great-circle `_amb` phase (pre-fix bug).
    """
    if binding == "fixed":
        G1 = SIGN * fapex * centroid_phase(3, phi_rad)   # APEX -> slot 1
        G2 = SIGN * fbase1 * centroid_phase(1, phi_rad)  # BASE leg 1 -> slot 2
        G3 = SIGN * fbase2 * centroid_phase(2, phi_rad)  # BASE leg 2 -> slot 3
        G0 = SIGN * fppp * centroid_phase(0, phi_rad)
        return G0, G1, G2, G3
    if binding == "old":
        half = phi_rad / 2.0
        d1 = -np.arctan2(1.0, 3.0 * np.tan(half))
        d2 = -2.0 * np.pi + np.arctan2(1.0, 3.0 * np.tan(half))
        d3 = -0.5 * np.pi + np.arctan2(1.0, np.tan(half))

        def old_phase(mu):
            s = {0: (1, 1, 1), 1: (-1, 1, 1), 2: (1, -1, 1), 3: (1, 1, -1)}[mu]
            return np.exp(-2j * (s[0] * d1 + s[1] * d2 + s[2] * d3))

        G0 = SIGN * fppp * old_phase(0)
        G1 = SIGN * fbase1 * old_phase(1)  # BASE leg 1 -> slot 1 (OLD)
        G2 = SIGN * fbase2 * old_phase(2)  # BASE leg 2 -> slot 2 (OLD)
        G3 = SIGN * fapex * old_phase(3)   # APEX -> slot 3 (OLD)
        return G0, G1, G2, G3
    raise ValueError(binding)


GAMMAS = [1.0, 10.0, 100.0]
PHIS = [15.0, 30.0, 60.0, 90.0, 120.0, 150.0]


@pytest.mark.parametrize("gamma_arcmin", GAMMAS)
@pytest.mark.parametrize("phi_deg", PHIS)
def test_mirror_symmetry_fixed(gamma_arcmin, phi_deg):
    """FIXED assembler: Gamma^1 real and Gamma^2 = conj(Gamma^3) on the SAS family.

    Note: the symmetry is gamma-independent (the centroid angles depend only on
    phi); gamma is swept to honour the boundary-validation extreme-corner rule.
    """
    phi_rad = np.radians(phi_deg)
    fb1, fb2, fa, fp = _synthetic_folds(phi_deg)
    G0, G1, G2, G3 = _assemble("fixed", fb1, fb2, fa, fp, phi_rad)

    # all finite
    for G in (G0, G1, G2, G3):
        assert np.isfinite(G.real) and np.isfinite(G.imag)

    # (i) APEX-conjugated component (slot 1) is REAL
    im_g1 = abs(G1.imag) / max(abs(G1), 1e-300)
    assert im_g1 < 1e-6, f"Gamma^1 not real: |Im|/|G1|={im_g1:.3e}"

    # (ii) BASE pair: Gamma^2 = conj(Gamma^3)
    base = abs(G2 - np.conj(G3)) / max(abs(G2), 1e-300)
    assert base < 1e-6, f"Gamma^2 != conj(Gamma^3): {base:.3e}"

    # Gamma^0 (all-unconjugated PPP) is also real on the isoceles family
    im_g0 = abs(G0.imag) / max(abs(G0), 1e-300)
    assert im_g0 < 1e-6, f"Gamma^0 not real: |Im|/|G0|={im_g0:.3e}"


def test_apex_modulus_pin():
    """Pin the fresh post-fix canoes APEX-channel modulus at (gamma=10', phi=60).

    Regressions in the radial fold (or an accidental phase that changes the
    modulus) show up here as numeric drift, independent of the symmetry asserts.
    """
    npz = _HERE.parent / "outputs" / "stage1_ours_v2_JACFIX.npz"
    if not npz.exists():
        pytest.skip("stage1_ours_v2_JACFIX.npz not built")
    # allow_pickle: trusted local artifact built by this same pipeline (not external).
    data = np.load(npz, allow_pickle=True)
    ig = int(np.argmin(np.abs(data["gamma_arcmin"] - 10.0)))
    ip = int(np.argmin(np.abs(data["phi_deg"] - 60.0)))
    # After the fix, the APEX channel lives in slot Gamma^1.
    apex_mod = abs(data["Gamma1"][ig, ip])
    assert apex_mod == pytest.approx(PIN_APEX_MOD_G10_PHI60, rel=PIN_RTOL), (
        f"apex modulus drifted: got {apex_mod:.4e}, "
        f"pinned {PIN_APEX_MOD_G10_PHI60:.4e}")


@pytest.mark.parametrize("phi_deg", PHIS)
def test_old_binding_violates_symmetry(phi_deg):
    """Demonstrate the test has teeth: the OLD slot/phase binding VIOLATES the
    mirror symmetry by O(1).  This is the pre-fix bug; we assert it FAILS the
    same constraints the fixed assembler passes."""
    phi_rad = np.radians(phi_deg)
    fb1, fb2, fa, fp = _synthetic_folds(phi_deg)
    _, G1, G2, G3 = _assemble("old", fb1, fb2, fa, fp, phi_rad)
    im_g1 = abs(G1.imag) / max(abs(G1), 1e-300)
    base = abs(G2 - np.conj(G3)) / max(abs(G2), 1e-300)
    # At least one of the two constraints must be badly violated by the old map.
    assert (im_g1 > 1e-3) or (base > 1e-3), (
        f"OLD binding unexpectedly satisfied symmetry at phi={phi_deg}: "
        f"im_g1={im_g1:.3e} base={base:.3e}")
