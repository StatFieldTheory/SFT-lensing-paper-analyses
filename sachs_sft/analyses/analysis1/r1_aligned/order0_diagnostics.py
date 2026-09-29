"""Bounded diagnostics of the fixed C table's affine interpolation.

Alternative interpolants are diagnostics, not approved replacement models.
No covariance table or production output is modified or rebuilt.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline, interp1d

HERE = Path(__file__).resolve().parent
SACHS = HERE.parents[2]
C_TABLE = SACHS / "callables/C_propagator/corr_op/corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz"
CANOES = Path("/Users/zzhang/projects/angular_statistics/canoes")


def integration_weights(nodes: np.ndarray, lower: float, upper: float, method: str) -> np.ndarray:
    """Integrate cardinal source interpolation functions numerically/exactly."""
    identity = np.eye(nodes.size)
    if method == "cubic":
        return CubicSpline(nodes, identity, axis=0).integrate(lower, upper)
    x = nodes
    y = identity
    if method == "linear_zero":
        x = np.r_[0.0, nodes]
        y = np.vstack([np.zeros(nodes.size), identity])
    support = np.r_[lower, x[(x > lower) & (x < upper)], upper]
    basis = interp1d(x, y, axis=0, fill_value="extrapolate")
    return np.trapezoid(basis(support), support, axis=0)


def check_high_ell_windows() -> None:
    """Factor the existing high-ell Weyl source window without rebuilding C."""
    sys.path.insert(0, str(CANOES / "src"))
    from canoes.cosmo.lightcone import build_lightcone, chi_of_lambda, z_of_chi
    from canoes.sachs.sft_input.corr_op.table import _build_shared_chi_grid
    from canoes.sachs.sft_input.corr_op.build import _load_pk_phi
    from canoes.sachs.sft_input.corr_op._limber import _id_basis_windows
    from canoes.nuell.transfer._sachs_decomp import get_sachs_decomposition
    from canoes.nuell.transfer.sachs_driver import COEFS_PSI0

    with np.load(C_TABLE) as table:
        nodes = np.asarray(table["lambda_grid"])
        ell_build = np.asarray(table["ell_grid"])
        xx = np.asarray(table["cl_xx"])
        md = json.loads(table["cl_table_metadata_json"].item())
        config = json.loads(table["config_json"].item())
    source = json.loads((HERE / "outputs/order0_table_integral/metadata.json").read_text())
    final = source["lambda_source_mpc"]
    chi_source = source["chi_source_mpc"]
    cosmo = SimpleNamespace(**config["cosmology"])
    source_chi = np.asarray(md["resolved_source_chi_grid"])
    chi = _build_shared_chi_grid(chi_min=50.0, lambda_grid=source_chi, n_chi_per_pair=256)
    lightcone, opz = build_lightcone(cosmo, chi, z_max=1200, distance_unit="Mpc")
    opz_sources = np.interp(source_chi, chi, opz)
    selection = np.unique([np.argmin(np.abs(ell_build - target)) for target in (1000, 2000, 4000)])
    ell = ell_build[selection]
    g, pref = _id_basis_windows(
        decomp=get_sachs_decomposition(COEFS_PSI0), chi=chi,
        lam=source_chi, opz_chi=opz, opz_lam=opz_sources,
        m_id=lightcone.M[0], lightcone=lightcone,
        l2_arr=ell.astype(float) * (ell + 1),
    )
    # Recover the inner-distance-only factor from a source that covers the grid.
    source_distance = source_chi[-1] / opz_sources[-1]
    inner_factor = g[0, -1] * source_distance**2
    source_floor_chi = float(chi_of_lambda(cosmo, np.array([nodes[0]]), distance_unit="Mpc", z_max=1200)[0])
    dense_chi = np.geomspace(50.0, chi_source, 100001)
    dense_z = z_of_chi(cosmo, dense_chi, distance_unit="Mpc", z_max=1200)
    dense_a = 1.0 / (1.0 + dense_z)
    dense_distance = dense_a * dense_chi
    integrand = dense_a**2 / dense_distance**2
    cumulative = cumulative_trapezoid(integrand, dense_chi, initial=0)
    lower = np.maximum(chi, source_floor_chi)
    tail = np.interp(lower, dense_chi, cumulative[-1] - cumulative, left=0, right=0)
    continuum_window = inner_factor * tail
    full_source_tail = np.interp(chi, dense_chi, cumulative[-1] - cumulative, left=0, right=0)
    weights = {
        name: integration_weights(nodes, float(nodes[0]), final, name)
        for name in ("linear", "cubic")
    }
    windows = {name: weight @ g[0] for name, weight in weights.items()}
    windows["continuum_source"] = continuum_window
    windows["continuum_source_no_artificial_floor"] = inner_factor * full_source_tail
    power = _load_pk_phi(
        CANOES / "examples/data/PCAMB_pyccl_stf_fid_z0.txt",
        omega_m=cosmo.Omega_m, h_param=cosmo.h, n_s=0.97, distance_unit="Mpc",
    )
    with np.load(HERE / "ccl_padded_n4096.npz") as ccl:
        reference = np.asarray(ccl["cl_ee"])
    results = {}
    for slot, (index, multipole) in enumerate(zip(selection, ell, strict=True)):
        kernel = power(multipole / chi) / chi**2
        row = {}
        for name, window in windows.items():
            value = float(pref[0, slot]**2 * np.trapezoid(kernel * window**2, chi))
            row[name] = {"cl_ee": value, "ratio_to_ccl": value / reference[multipole - 2]}
            if name in weights:
                table_value = float(weights[name] @ xx[index] @ weights[name])
                row[name]["factorization_ratio_to_stored_table"] = value / table_value
        for name in ("continuum_source", "continuum_source_no_artificial_floor"):
            window = windows[name]
            refined_kernel = power((multipole + 0.5) / chi) / chi**2
            value = float(pref[0, slot]**2 * np.trapezoid(refined_kernel * window**2, chi))
            row[name + "_ell_plus_half"] = {"cl_ee": value, "ratio_to_ccl": value / reference[multipole - 2]}
        results[str(multipole)] = row
    output = {"scope": "High-ell Weyl ID branch only, identical stored source and inner-radial grids. Source-continuum window numerically integrates the actual response factor before forming its covariance.",
              "source_floor_chi_mpc": source_floor_chi, "chi_source_mpc": chi_source,
              "results": results}
    path = HERE / "outputs/high_ell_window_diagnostic.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


def main() -> None:
    if "--high-ell-window" in sys.argv:
        check_high_ell_windows()
        return
    sys.path.insert(0, str(CANOES / "src"))
    from canoes.sachs.sft_input.corr_op.table import _interpolate_cl_to_dense_ell
    from ccl_reference import transform

    out = HERE / "outputs" / "interpolation_diagnostics_v2"
    out.mkdir(exist_ok=False)
    source = json.loads((HERE / "outputs/order0_table_integral/metadata.json").read_text())
    upper = float(source["lambda_source_mpc"])
    with np.load(C_TABLE) as table:
        nodes = np.asarray(table["lambda_grid"])
        ell_build = np.asarray(table["ell_grid"])
        ell = np.asarray(table["ell_assembly_grid"])
        saved = {name: np.asarray(table[name]) for name in ("cl_pp", "cl_px", "cl_xx")}
    weights = {
        "linear": integration_weights(nodes, float(nodes[0]), upper, "linear"),
        "cubic": integration_weights(nodes, float(nodes[0]), upper, "cubic"),
        "linear_zero": integration_weights(nodes, 0.0, upper, "linear_zero"),
        "linear_extrapolate": integration_weights(nodes, 0.0, upper, "linear"),
    }
    spectra = {method: {} for method in weights}
    for channel, cube in saved.items():
        dense = _interpolate_cl_to_dense_ell(cube, ell_build, ell)
        for method, weight in weights.items():
            spectra[method][channel] = np.einsum("i,lij,j->l", weight, dense, weight, optimize=True)
    with np.load(HERE / "outputs/order0_table_integral/order0_observables.npz") as reference:
        theta = np.asarray(reference["theta_arcmin"])
        baseline = {key: np.asarray(reference[key]) for key in ("kk", "xip", "xim", "kgt")}
    ccl_path = HERE / "ccl_padded_n4096.npz"
    with np.load(ccl_path) as ccl:
        ccl_observables = {key: np.asarray(ccl[key]) for key in baseline}
        ccl_spectra = {"cl_pp": np.asarray(ccl["cl_kk"]), "cl_px": -np.asarray(ccl["cl_ke"]), "cl_xx": np.asarray(ccl["cl_ee"])}
    summary = {
        "source_lambda_mpc": upper,
        "interpretation": "Alternative interpolation and extrapolation checks of a fixed covariance table, not physical predictions or convergence claims.",
        "ccl_reference": str(ccl_path),
        "methods": {},
        "actual_low_affine_policy": "Zero when either source affine coordinate is below 406 Mpc. Earlier assumed extrapolation was incorrect.",
        "batch_geometry_policy": "The current canoes batch implementation rounds cos(theta) to 8 decimals before angular-channel construction, unlike the scalar callable.",
    }
    names = {"kk": "xi_kappa", "xip": "xi_plus", "xim": "xi_minus", "kgt": "xi_kappa_gamma"}
    for method, cl in spectra.items():
        transformed, _ = transform(ell, cl["cl_pp"], cl["cl_xx"], -cl["cl_px"], theta, 0.2)
        observable = {key: transformed[value] for key, value in names.items()}
        if method == "linear":
            rounded_theta = np.rad2deg(np.arccos(np.round(np.cos(np.deg2rad(theta / 60)), 8))) * 60
            rounded, _ = transform(ell, cl["cl_pp"], cl["cl_xx"], -cl["cl_px"], rounded_theta, 0.2)
            rounded_observable = {key: rounded[value] for key, value in names.items()}
            np.savez_compressed(out / "linear_batch_rounded.npz", theta_arcmin=theta,
                                effective_theta_arcmin=rounded_theta, ell=ell,
                                **rounded_observable, **cl)
        result = {}
        for key in baseline:
            peak = float(np.max(np.abs(ccl_observables[key])))
            result[key] = {
                "small_angle": float(observable[key][0]),
                "change_from_linear_over_ccl_peak": float(np.max(np.abs(observable[key] - baseline[key])) / peak),
                "max_residual_over_ccl_peak": float(np.max(np.abs(observable[key] - ccl_observables[key])) / peak),
                "small_angle_ratio_to_ccl": float(observable[key][0] / ccl_observables[key][0]),
            }
        result["spectral_ratio_to_ccl"] = {
            str(multipole): {channel: float(cl[channel][multipole - 2] / ccl_spectra[channel][multipole - 2]) for channel in cl}
            for multipole in (100, 500, 800, 1000, 2000, 4000)
        }
        summary["methods"][method] = result
        np.savez_compressed(out / f"{method}.npz", theta_arcmin=theta, ell=ell, **observable, **cl, affine_weights=weights[method])
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
