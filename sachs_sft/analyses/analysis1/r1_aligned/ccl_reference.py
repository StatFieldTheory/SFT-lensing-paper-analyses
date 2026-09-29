"""Independent CCL projection with the shared inputs of the R1 comparison.

The matter spectrum, background, and growth are prescribed inputs. CCL alone
performs the radial angular-spectrum integration. The convergence and shear
tracers differ through CCL's documented angular prefactors, obtained from its
API rather than approximated as equal. No nonlinear optical vertices enter.

Official API references:
https://ccl.readthedocs.io/en/latest/api/pyccl.cosmology.html
https://ccl.readthedocs.io/en/latest/api/pyccl.tracers.html
https://ccl.readthedocs.io/en/latest/api/pyccl.cells.html
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyccl as ccl
from scipy.interpolate import CubicSpline
from scipy.special import eval_legendre

HERE = Path(__file__).resolve().parent
SACHS_ROOT = HERE.parents[2]
DEFAULT_CANOES = Path("/Users/zzhang/projects/angular_statistics/canoes")
OMEGA_M = 0.3160919980475834
H = 0.6711
OMEGA_B = 0.02207 / H**2


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(args: argparse.Namespace):
    """Load shared physical inputs and construct independent CCL tracers."""
    sys.path.insert(0, str(args.canoes_root / "src"))
    from canoes.cosmo.growth import growth_D
    from canoes.cosmo.lightcone import _zchi_growth_tables

    raw = np.loadtxt(args.power_table)
    if raw.ndim != 2 or raw.shape[1] != 2:
        raise ValueError("The input power table must have two columns.")
    k = raw[:, 0] * H
    p0 = raw[:, 1] / H**3
    if np.any(np.diff(k) <= 0) or np.any(p0 <= 0):
        raise ValueError("The input power table must be sorted and positive.")

    z, chi, hubble, growth = _zchi_growth_tables(
        SimpleNamespace(Omega_m=OMEGA_M, h=H),
        z_max=args.background_zmax,
        n_z=args.background_nz,
        distance_unit="Mpc",
    )
    a = (1.0 / (1.0 + z))[::-1].copy()
    d = growth[::-1].copy()
    growth_rate = CubicSpline(np.log(a), np.log(d)).derivative()(np.log(a))
    a_pk = np.geomspace(a[0], 1.0, args.power_na)
    d_pk = np.asarray(growth_D(1.0 / a_pk - 1.0, Omega_m=OMEGA_M))
    pk = {
        "a": a_pk,
        "k": k,
        "delta_matter:delta_matter": d_pk[:, None] ** 2 * p0[None, :],
    }
    cosmo = ccl.CosmologyCalculator(
        Omega_c=OMEGA_M - OMEGA_B,
        Omega_b=OMEGA_B,
        h=H,
        n_s=0.97,
        sigma8=0.81,
        Omega_g=0.0,
        Neff=0.0,
        m_nu=0.0,
        background={"a": a, "chi": chi[::-1].copy(), "h_over_h0": hubble[::-1].copy()},
        growth={"a": a, "growth_factor": d, "growth_rate": growth_rate},
        pk_linear=pk,
        pk_nonlin=pk,
    )
    full_chi, full_kernel = ccl.get_kappa_kernel(
        cosmo, z_source=args.z_s, n_samples=args.kernel_samples
    )
    if not 0 < args.chi_min < full_chi[-1]:
        raise ValueError("chi_min must lie strictly inside the source distance.")
    # FFTLog needs a wide, zero-padded radial range. Starting its logarithmic
    # range at the nonzero physical cut creates severe low-ell edge aliasing.
    # Represent the sharp physical cut by a tiny resolved transition, retaining
    # the original smooth kernel above it and exact zeros below the transition.
    edge_left = args.chi_min - args.kernel_edge_width
    if edge_left <= args.fkem_chi_min:
        raise ValueError("The physical cut must exceed the padded FFTLog minimum.")
    edge_nodes = args.chi_min + args.kernel_edge_width * np.array([0, 1, 2, 4, 8, 16, 32, 64])
    radial_chi = np.unique(np.r_[0.0, np.geomspace(args.fkem_chi_min, edge_left, 256),
                                edge_nodes, full_chi[full_chi > args.chi_min]])
    radial_kernel = np.interp(radial_chi, full_chi, full_kernel)
    radial_kernel[radial_chi <= edge_left] = 0.0
    pad_max = getattr(args, "radial_zero_pad_max", None)
    if pad_max is not None:
        if not full_chi[-1] < pad_max < chi[-1]:
            raise ValueError("Upper zero padding must exceed the source and remain inside the supplied background.")
        padding = np.geomspace(full_chi[-1] * (1.0 + 1e-9), pad_max, 256)
        radial_chi = np.r_[radial_chi, padding]
        radial_kernel = np.r_[radial_kernel, np.zeros(padding.size)]
    convergence = ccl.Tracer()
    shear = ccl.Tracer()
    for tracer, derivative in ((convergence, 1), (shear, 2)):
        tracer.add_tracer(
            cosmo,
            kernel=(radial_chi, radial_kernel),
            der_bessel=-1,
            der_angles=derivative,
        )
    metadata = {
        "cosmology": {"Omega_m": OMEGA_M, "Omega_b": OMEGA_B, "h": H,
                      "n_s": 0.97, "sigma8_label_of_input_table": 0.81,
                      "Omega_g": 0.0, "Neff": 0.0, "m_nu": 0.0},
        "background": "canoes matter-plus-Lambda input arrays in physical Mpc",
        "growth": "canoes.cosmo.growth.growth_D, normalized at z=0",
        "power": "same tabulated z=0 linear CAMB power on both sides, scaled by D(a)^2",
        "power_table": str(args.power_table.resolve()),
        "power_table_sha256": sha256(args.power_table),
        "power_table_units": {"input_k": "h/Mpc", "input_P": "(Mpc/h)^3",
                              "ccl_k": "1/Mpc", "ccl_P": "Mpc^3"},
        "power_k_range_inverse_mpc": [float(k[0]), float(k[-1])],
        "power_extrapolation": "CCL Pk2D defaults, low-k order 1, high-k order 2",
        "z_s": args.z_s,
        "chi_source_mpc": float(full_chi[-1]),
        "radial_support_mpc": [float(edge_left), float(full_chi[-1])],
        "radial_tabulation_range_mpc": [float(radial_chi[0]), float(radial_chi[-1])],
        "upper_zero_padding_mpc": pad_max,
        "physical_radial_cut_mpc": args.chi_min,
        "radial_cut_transition_width_mpc": args.kernel_edge_width,
        "kernel": "CCL get_kappa_kernel with shared physical chi_min and explicit zero padding outside physical support",
        "tracer_der_bessel": -1,
        "convergence_der_angles": 1,
        "shear_der_angles": 2,
        "source_hashes": {
            str(Path(__file__).resolve()): sha256(Path(__file__)),
            str(args.canoes_root / "src/canoes/cosmo/lightcone.py"):
                sha256(args.canoes_root / "src/canoes/cosmo/lightcone.py"),
            str(args.canoes_root / "src/canoes/cosmo/growth.py"):
                sha256(args.canoes_root / "src/canoes/cosmo/growth.py"),
        },
    }
    return cosmo, convergence, shear, radial_chi, radial_kernel, metadata


def transform(ells, cl_kk, cl_ee, cl_ke, gamma_arcmin, taper_fraction):
    """Apply the same finite full-sky angular band and taper to all spectra."""
    sys.path.insert(0, str(SACHS_ROOT / "scripts"))
    from spin_rotation import _wigner_d_array

    taper = np.ones_like(ells, dtype=float)
    start = float(ells[0] + (1.0 - taper_fraction) * (ells[-1] - ells[0]))
    if taper_fraction > 0:
        mask = ells > start
        taper[mask] = 0.5 * (1.0 + np.cos(np.pi * (ells[mask] - start) / (ells[-1] - start)))
    weight = (2.0 * ells + 1.0) * taper / (4.0 * np.pi)
    output = {key: np.empty(gamma_arcmin.size) for key in
              ("xi_kappa", "xi_plus", "xi_minus", "xi_kappa_gamma")}
    for index, angle in enumerate(np.deg2rad(gamma_arcmin / 60.0)):
        cosine = np.cos(angle)
        output["xi_kappa"][index] = np.sum(weight * cl_kk * eval_legendre(ells, cosine))
        output["xi_plus"][index] = np.sum(weight * cl_ee * _wigner_d_array(ells, 2, 2, cosine))
        output["xi_minus"][index] = np.sum(weight * cl_ee * _wigner_d_array(ells, 2, -2, cosine))
        output["xi_kappa_gamma"][index] = np.sum(weight * cl_ke * _wigner_d_array(ells, 0, 2, cosine))
    return output, taper


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--canoes-root", type=Path, default=DEFAULT_CANOES)
    parser.add_argument("--power-table", type=Path)
    parser.add_argument("--z-s", type=float, default=5.0)
    parser.add_argument("--chi-min", type=float, default=50.0)
    parser.add_argument("--ell-max", type=int, default=5000)
    parser.add_argument("--l-limber", type=int, default=800)
    parser.add_argument("--fkem-nchi", type=int, default=1024)
    parser.add_argument("--fkem-chi-min", type=float, default=1e-6)
    parser.add_argument("--kernel-edge-width", type=float, default=1e-3)
    parser.add_argument("--kernel-samples", type=int, default=4096)
    parser.add_argument("--radial-zero-pad-max", type=float,
                        help="Extend the identically zero kernel beyond the source to lower FKEM's finite k minimum")
    parser.add_argument("--background-zmax", type=float, default=1200.0)
    parser.add_argument("--background-nz", type=int, default=5001)
    parser.add_argument("--power-na", type=int, default=256)
    parser.add_argument("--taper-fraction", type=float, default=0.2)
    parser.add_argument("--angle-file", type=Path, help="NPZ containing the shared theta_arcmin grid.")
    parser.add_argument("--probe", action="store_true", help="Only a few multipoles, no real-space sums.")
    args = parser.parse_args()
    if args.power_table is None:
        args.power_table = args.canoes_root / "examples/data/PCAMB_pyccl_stf_fid_z0.txt"
    if args.out.suffix != ".npz":
        parser.error("--out must end in .npz")
    if args.out.exists() or args.out.with_suffix(".json").exists():
        parser.error("Output exists. Choose a new name to preserve earlier artifacts.")
    if args.z_s <= 0 or args.z_s > args.background_zmax or args.ell_max < 3:
        parser.error("Require 0 < z_s <= background_zmax and ell_max >= 3")
    if not 0 <= args.taper_fraction <= 1:
        parser.error("Require 0 <= taper_fraction <= 1")
    start = time.monotonic()
    print(f"CCL {ccl.__version__}: preparing prescribed input arrays", flush=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cosmo, convergence, shear, radial_chi, radial_kernel, metadata = prepare(args)
        if args.probe:
            ells = np.unique([2, 10, 100, min(args.l_limber, args.ell_max), args.ell_max])
            ells = ells[(ells >= 2) & (ells <= args.ell_max)]
        else:
            ells = np.arange(2, args.ell_max + 1)
        print(f"Computing {ells.size} multipoles with FKEM through ell={args.l_limber}", flush=True)
        cl_kk = np.asarray(ccl.angular_cl(
            cosmo, convergence, convergence, ells.astype(float),
            l_limber=args.l_limber,
            non_limber_integration_method="FKEM",
            fkem_Nchi=args.fkem_nchi,
            fkem_chi_min=args.fkem_chi_min,
        ))
        f_kappa = convergence.get_f_ell(ells)[0]
        f_shear = shear.get_f_ell(ells)[0]
        factor = f_shear / f_kappa
        cl_ee = cl_kk * factor**2
        cl_ke = cl_kk * factor
        k_ranges = []
        for multipole in np.unique([2, 10, 100, min(args.l_limber, args.ell_max)]):
            if multipole not in ells or multipole > args.l_limber:
                continue
            cached_k, _ = convergence._get_fkem_fft(
                convergence._trc[0], args.fkem_nchi, args.fkem_chi_min,
                float(radial_chi[-1]), int(multipole), cosmo,
            )
            if cached_k is not None:
                k_ranges.append({"ell": int(multipole), "k_min_inverse_mpc": float(cached_k[0]),
                                 "k_max_inverse_mpc": float(cached_k[-1]), "n_k": int(cached_k.size)})
        metadata["fkem_output_k_ranges"] = k_ranges
        metadata["finite_k_range_caveat"] = (
            "FKEM integrates only its returned finite k grid. Upper radial zero padding "
            "reduces, but does not eliminate, the omitted low-k tail. This numerical "
            "effect must not be attributed to the Ricci operator."
        )
        payload = {
            "ells": ells, "Cl_kappa": cl_kk, "Cl_EE": cl_ee, "Cl_kappa_E": cl_ke,
            "ell": ells, "cl_kk": cl_kk, "cl_ee": cl_ee, "cl_ke": cl_ke,
            "f_kappa": f_kappa, "f_shear": f_shear,
            "radial_chi_mpc": radial_chi, "radial_kernel": radial_kernel,
        }
        if args.probe:
            # Check the API-prefactor reuse against separately constructed tracers.
            direct_ee = ccl.angular_cl(cosmo, shear, shear, ells.astype(float),
                                      l_limber=args.l_limber, fkem_Nchi=args.fkem_nchi,
                                      fkem_chi_min=args.fkem_chi_min)
            direct_ke = ccl.angular_cl(cosmo, convergence, shear, ells.astype(float),
                                      l_limber=args.l_limber, fkem_Nchi=args.fkem_nchi,
                                      fkem_chi_min=args.fkem_chi_min)
            np.testing.assert_allclose(cl_ee, direct_ee, rtol=2e-9, atol=0)
            np.testing.assert_allclose(cl_ke, direct_ke, rtol=2e-9, atol=0)
            metadata["direct_tracer_prefactor_check"] = "passed, rtol=2e-9"
        else:
            if args.angle_file is None:
                import yaml

                config_path = SACHS_ROOT / "sftwick_outputs/2PCF/C_corr_op_O0/config_L2.yaml"
                config = yaml.safe_load(config_path.read_text())
                positions = config["sweep"]["positions_grid"]
                first = np.asarray(positions["x"][0], dtype=float)
                second = np.asarray(positions["y"], dtype=float)
                first /= np.linalg.norm(first)
                second /= np.linalg.norm(second, axis=1)[:, None]
                gamma_arcmin = np.sort(np.rad2deg(np.arccos(np.clip(second @ first, -1, 1))) * 60)
                metadata["angle_config_sha256"] = sha256(config_path)
                metadata["angle_config"] = str(config_path)
            else:
                with np.load(args.angle_file, allow_pickle=False) as angles:
                    gamma_arcmin = np.asarray(angles["theta_arcmin"], dtype=float)
                metadata["angle_file_sha256"] = sha256(args.angle_file)
            if gamma_arcmin.ndim != 1 or np.any(np.diff(gamma_arcmin) <= 0):
                raise ValueError("The angular grid must be one-dimensional and ascending.")
            correlations, taper = transform(ells, cl_kk, cl_ee, cl_ke, gamma_arcmin, args.taper_fraction)
            payload.update(correlations, gamma_arcmin=gamma_arcmin, theta_arcmin=gamma_arcmin,
                           kk=correlations["xi_kappa"], xip=correlations["xi_plus"],
                           xim=correlations["xi_minus"], kgt=correlations["xi_kappa_gamma"],
                           ell_taper=taper)
        if not all(np.all(np.isfinite(value)) for value in payload.values()):
            raise RuntimeError("Non-finite output detected")
        metadata.update(
            pyccl_version=ccl.__version__,
            settings={key: str(value) if isinstance(value, Path) else value
                      for key, value in vars(args).items()},
            integration="FKEM for ell <= l_limber and Limber for ell > l_limber",
            elapsed_seconds=time.monotonic() - start,
            warnings=[{"category": item.category.__name__, "message": str(item.message)} for item in caught],
            angular_spectra="CCL convergence spectrum with exact CCL tracer angular factors for shear/cross",
            real_space_transform="direct Legendre/Wigner-d sums, no flat-sky approximation",
            optical_order="linear projection, no Sachs F or driving-cumulant K insertions",
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, **payload, metadata_json=json.dumps(metadata, sort_keys=True))
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"out": str(args.out), "elapsed_seconds": metadata["elapsed_seconds"],
                      "warnings": metadata["warnings"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
