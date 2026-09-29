"""Record a finite-angle pure-E round trip using the accepted Order-0 product.

The stored xi-plus and xi-minus were generated from the same tapered cl_ww
spectrum, with no B input. This control measures recovery on the existing
40-angle grid. It is not an error bound for a different input such as FK.
Run with the PyCCL Python environment. The output must be a fresh JSON path.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS3 = HERE.parent / "analysis3"
sys.path.insert(0, str(ANALYSIS3))
import plot_analysis3_cl_decomposition as transform  # noqa: E402

DEFAULT_SOURCE = (
    HERE.parent / "analysis1/r1_aligned"
    / "source_folded_full_boundary_n4096_k2048.npz"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--out", type=Path, default=HERE / "pure_e_transform_check.json")
    args = parser.parse_args()
    if args.out.suffix != ".json":
        parser.error("--out must have the .json suffix")
    if args.out.exists():
        parser.error(f"output exists, choose a fresh path: {args.out}")

    with np.load(args.source, allow_pickle=False) as data:
        ell = np.asarray(data["ell"], dtype=int)
        gamma = np.asarray(data["gamma_arcmin"], dtype=float)
        cl_ee = np.asarray(data["cl_ww"], dtype=float)
        taper = np.asarray(data["ell_taper"], dtype=float)
        xi_plus = np.asarray(data["xi_plus"], dtype=float)
        xi_minus = np.asarray(data["xi_minus"], dtype=float)
    if not np.array_equal(ell, np.arange(2, 5001)):
        raise ValueError("expected the accepted integer ell = 2..5000 band")
    if gamma.shape != (40,) or not np.all(np.diff(gamma) > 0):
        raise ValueError("expected the increasing 40-angle production grid")
    if xi_plus.shape != gamma.shape or xi_minus.shape != gamma.shape:
        raise ValueError("stored spin correlation shapes differ from the angle grid")
    if cl_ee.shape != ell.shape or taper.shape != ell.shape:
        raise ValueError("stored spectrum or taper shape differs from ell")
    if not all(np.all(np.isfinite(arr)) for arr in (gamma, cl_ee, taper, xi_plus, xi_minus)):
        raise ValueError("source contains non-finite values")

    taper_start = float(ell[0] + 0.8 * (ell[-1] - ell[0]))
    expected_taper = np.ones(ell.size)
    tapered = ell > taper_start
    expected_taper[tapered] = 0.5 * (
        1.0 + np.cos(np.pi * (ell[tapered] - taper_start) / (ell[-1] - taper_start))
    )
    if not np.array_equal(taper, expected_taper):
        raise ValueError("source taper differs from the documented cosine window")

    production_ell = transform.ELL.astype(int)
    true_ee = (cl_ee * taper)[production_ell - ell[0]]
    if np.any(true_ee <= 0):
        raise ValueError("this relative-residual control requires positive input EE")
    band_prefactor = production_ell * (production_ell + 1) / (2.0 * np.pi)
    high = production_ell >= 50
    peak = float(np.max(band_prefactor[high] * true_ee[high]))
    recoveries = {}
    for n_fine in (20000, 40000):
        plus_setup = transform.build_curved_matrix(
            gamma, production_ell, 2, 2, n_fine=n_fine, apodise=False,
        )
        minus_setup = transform.build_curved_matrix(
            gamma, production_ell, 2, -2, n_fine=n_fine, apodise=False,
        )
        plus = transform.forward_curved(xi_plus, plus_setup, dc_subtract=False)
        minus = transform.forward_curved(xi_minus, minus_setup, dc_subtract=False)
        recovered_ee = 0.5 * (plus + minus)
        recovered_bb = 0.5 * (plus - minus)
        recoveries[str(n_fine)] = {
            "recovered_ee": recovered_ee.tolist(),
            "recovered_bb": recovered_bb.tolist(),
            "ee_minus_input": (recovered_ee - true_ee).tolist(),
            "recovered_ee_over_input_ee": (recovered_ee / true_ee).tolist(),
            "recovered_bb_over_input_ee": (recovered_bb / true_ee).tolist(),
            "max_abs_bb_bandpower_over_input_ee_peak_ell_ge_50": float(
                np.max(np.abs(band_prefactor[high] * recovered_bb[high])) / peak
            ),
        }

    coarse, fine = recoveries["20000"], recoveries["40000"]
    difference = {}
    for key in ("recovered_ee", "recovered_bb"):
        delta = np.asarray(fine[key]) - np.asarray(coarse[key])
        difference[key + "_fine_minus_coarse"] = delta.tolist()
        difference[key + "_max_abs_bandpower_change_over_input_ee_peak_ell_ge_50"] = float(
            np.max(np.abs(band_prefactor[high] * delta[high])) / peak
        )
    result = {
        "schema_version": 1,
        "control": "Stored pure-E Order-0 xi-plus and xi-minus on the production angle grid",
        "scope": (
            "Residuals include finite angular coverage and interpolation of this input. "
            "The integration-grid comparison does not test missing angular samples. "
            "These results are not a universal FK error bound or a subtraction template."
        ),
        "source": str(args.source.resolve()),
        "source_sha256": sha256(args.source),
        "script_sha256": sha256(Path(__file__)),
        "transform_sha256": sha256(Path(transform.__file__)),
        "source_spectrum": "cl_ww * ell_taper",
        "source_ell_band": [int(ell[0]), int(ell[-1])],
        "source_taper": {
            "kind": "cosine over the final 20 percent of the ell band",
            "starts_at_ell": taper_start,
            "ends_at_ell": int(ell[-1]),
            "array_sha256": hashlib.sha256(taper.tobytes()).hexdigest(),
        },
        "gamma_arcmin": gamma.tolist(),
        "ell": production_ell.tolist(),
        "true_ee": true_ee.tolist(),
        "true_bb": np.zeros(production_ell.size).tolist(),
        "transform_settings": {"dc_subtract": False, "apodise": False},
        "recoveries": recoveries,
        "integration_grid_comparison": difference,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Wrote {args.out}")
    for n_fine, recovery in recoveries.items():
        print(n_fine, recovery["max_abs_bb_bandpower_over_input_ee_peak_ell_ge_50"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
