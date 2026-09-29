"""Project the registered source-folded Sachs operators with independent CCL.

The canoes coefficient decomposition defines the operator. CCL evaluates its
radial projection, including first and second Bessel derivatives, without the
canoes covariance table or its t-form angular-spectrum assembler. No amplitude
calibration is applied. Outputs retain the negative Ricci-Weyl cross sign.
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

import ccl_reference as standard
import source_folded_reference as folded

HERE = Path(__file__).resolve().parent


def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def angular_mode(term, cosmo, chi_s):
    """Identify a registered prefactor by comparison with actual CCL values."""
    probes = np.array([2, 3, 10, 37, 100, 500, 800, 1300], dtype=float)
    target = np.asarray(term.ell_prefactor_fn(probes * (probes + 1)))
    for mode in range(3):
        tracer = ccl.Tracer()
        tracer.add_tracer(cosmo, kernel=(np.array([0., chi_s]), np.ones(2)),
                          der_angles=mode)
        candidate = np.asarray(tracer.get_f_ell(probes)[0])
        for sign in (1., -1.):
            # CCL uses an approximation to its spin-2 angular factor. Retain
            # its documented implementation and record the sub-2-ppm mismatch.
            if np.allclose(target, sign * candidate, rtol=2e-6, atol=0):
                error = float(np.max(np.abs(sign * candidate / target - 1)))
                return mode, sign, error
    raise ValueError(f"Unsupported angular prefactor in {term.label}")


def build_tracers(args, cosmo, chi_s, chi, lc, opz, vectors):
    from canoes.nuell.cls._compute_cl_windowed_sachs import _effective_window_on_chi
    from canoes.nuell.transfer._sachs_decomp import get_sachs_decomposition
    from canoes.nuell.transfer.basis import BASIS_TO_AB, DERIV_BASIS
    from canoes.sachs.sft_input.corr_op.build import C_KM_S

    # Zero padding protects FKEM's logarithmic radial transform at the sharp
    # physical chi_min cut. The nonzero windows are independently sampled from
    # the registered coefficients. No source-table interpolation enters.
    edge_left = args.chi_min - args.kernel_edge_width
    edge_nodes = args.chi_min + args.kernel_edge_width * np.array([0, 1, 2, 4, 8, 16, 32, 64])
    radial_chi = np.unique(np.r_[0., np.geomspace(args.fkem_chi_min, edge_left, 256),
                                edge_nodes, chi])
    valid = radial_chi >= args.chi_min
    q = radial_chi[valid]
    # Use the same lightcone interpolation convention as canoes' shifted-window
    # helper. The dense chi grid provides H, Hprime, M[a] and redshift values.
    q_opz = np.interp(q, chi, opz)
    q_a = ccl.scale_factor_of_chi(cosmo, q)
    q_growth = np.asarray(ccl.growth_factor(cosmo, q_a))
    poisson_amplitude = 1.5 * standard.OMEGA_M * (100 * standard.H / C_KM_S) ** 2
    lk = np.linspace(np.log(1e-8), np.log(1e4), 256)
    tracers, records, kernels = {}, {}, {}
    for name, vector in vectors.items():
        tracer_terms = []
        records[name] = []
        for basis_index, terms in enumerate(get_sachs_decomposition(vector)):
            time_order, radial_order = BASIS_TO_AB[basis_index]
            for term in terms:
                angle_order, sign, angle_error = angular_mode(term, cosmo, chi_s)
                effective = _effective_window_on_chi(
                    q, q_opz, chi_s, float(opz[-1]),
                    np.asarray(term.chi_kernel_fn(q, lc)),
                    np.interp(q, chi, lc.M[time_order]), chi_pow=0,
                    apply_density_correction=True,
                )
                kernel = np.zeros_like(radial_chi)
                # FKEM multiplies each radial window by matter growth. M[a] in
                # effective already supplies the desired potential evolution.
                kernel[valid] = sign * effective / q_growth
                # This is precisely the _load_pk_phi Poisson amplitude plus
                # the k power associated with the registered radial derivative.
                bessel_order = radial_order
                transfer_power = radial_order - 2
                if radial_order == 0 and args.id_bessel == "lensing":
                    # CCL explicitly provides j_l(x)/x^2 for potential-like
                    # lensing terms. Move chi^2 into the kernel and remove k^-2
                    # from the transfer, preserving the same radial integrand.
                    bessel_order = -1
                    transfer_power = 0
                    kernel[valid] *= q**2
                log_transfer = np.log(poisson_amplitude) + transfer_power * lk
                tracer = ccl.Tracer()
                tracer.add_tracer(
                    cosmo, kernel=(radial_chi, kernel), transfer_k=(lk, log_transfer),
                    der_bessel=bessel_order, der_angles=angle_order, is_logt=True,
                    extrap_order_lok=1, extrap_order_hik=1,
                )
                tracer_terms.append(tracer)
                records[name].append({
                    "label": term.label, "basis": DERIV_BASIS[basis_index],
                    "time_derivative": time_order, "radial_derivative": radial_order,
                    "ccl_der_bessel": bessel_order, "ccl_der_angles": angle_order,
                    "angular_sign_in_kernel": sign, "transfer_k_power": transfer_power,
                    "angular_prefactor_relative_error_at_probes": angle_error,
                    "poisson_amplitude": poisson_amplitude,
                })
                kernels[term.label] = kernel
        tracers[name] = tracer_terms
    return tracers, radial_chi, records, kernels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--ells", default="2,10,100,500,800")
    parser.add_argument("--n-chi", type=int, default=4096)
    parser.add_argument("--source-nchi", type=int, default=100001)
    parser.add_argument("--fkem-nchi", type=int, default=4096)
    parser.add_argument("--chi-min", type=float, default=50.)
    parser.add_argument("--kernel-edge-width", type=float, default=1e-3)
    parser.add_argument("--fkem-chi-min", type=float, default=1e-6)
    parser.add_argument("--id-bessel", choices=("raw", "lensing"), default="lensing",
                        help="Use CCL's documented j_l(x)/x^2 representation for radial-order-zero terms")
    args = parser.parse_args()
    if args.out.suffix != ".npz":
        parser.error("--out must end in .npz")
    if args.out.exists() or args.out.with_suffix(".json").exists():
        parser.error("Output exists. Choose a new filename to preserve artifacts.")
    if not 0 < args.fkem_chi_min < args.chi_min - args.kernel_edge_width:
        parser.error("Require 0 < fkem_chi_min < chi_min - kernel_edge_width")
    ell = np.unique([int(x) for x in args.ells.split(",")])
    if ell.min() < 2 or args.n_chi < 32 or args.fkem_nchi < 64:
        parser.error("Require ell >= 2 and adequate radial grids")
    begin = time.monotonic()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        chi_s, chi, lc, opz, vectors, _, source_metadata = folded.setup(args)
        ccl_args = SimpleNamespace(
            canoes_root=folded.CANOES, power_table=Path(source_metadata["power_table"]),
            background_zmax=1200., background_nz=5001, power_na=256, z_s=5.,
            kernel_samples=4096, chi_min=args.chi_min,
            kernel_edge_width=args.kernel_edge_width, fkem_chi_min=args.fkem_chi_min,
        )
        cosmo, standard_kappa, standard_shear, _, _, standard_metadata = standard.prepare(ccl_args)
        if abs(standard_metadata["chi_source_mpc"] - chi_s) > 1e-7:
            raise RuntimeError("The two source distances differ")
        tracers, radial_chi, mapping, kernels = build_tracers(
            args, cosmo, chi_s, chi, lc, opz, vectors,
        )
        tracers.update(standard_kappa=[standard_kappa], standard_shear=[standard_shear])
        pairs = {"cl_pp": ("phi", "phi"), "cl_pw": ("phi", "weyl"),
                 "cl_ww": ("weyl", "weyl"), "cl_tt": ("transverse", "transverse"),
                 "cl_tw": ("transverse", "weyl"),
                 "cl_standard_kk": ("standard_kappa", "standard_kappa")}
        payload = {"ell": ell, "radial_chi_mpc": radial_chi}
        timings = {}
        pair_cache = {}
        for name, (left, right) in pairs.items():
            start = time.monotonic()
            print(f"Computing {name}: {ell.tolist()}", flush=True)
            payload[name] = np.zeros(ell.size)
            # Installed CCL FKEM cannot broadcast get_transfer for a multi-term
            # Tracer. Sum the exact bilinear expansion of single-term pairs.
            # Symmetric pair caching avoids duplicate projections.
            for tr_left in tracers[left]:
                for tr_right in tracers[right]:
                    key = tuple(sorted((id(tr_left), id(tr_right))))
                    if key not in pair_cache:
                        pair_cache[key] = np.asarray(ccl.angular_cl(
                            cosmo, tr_left, tr_right, ell.astype(float),
                            l_limber=int(ell.max()) + 1, non_limber_integration_method="FKEM",
                            fkem_Nchi=args.fkem_nchi, fkem_chi_min=args.fkem_chi_min,
                        ))
                    payload[name] += pair_cache[key]
            if not np.all(np.isfinite(payload[name])):
                raise RuntimeError(f"Nonfinite projection in {name}")
            timings[name] = time.monotonic() - start
            print(f"{name}: {payload[name].tolist()}, {timings[name]:.3f}s", flush=True)
        spin_ratio = standard_shear.get_f_ell(ell)[0] / standard_kappa.get_f_ell(ell)[0]
        payload["cl_standard_ee"] = payload["cl_standard_kk"] * spin_ratio**2
        payload["cl_standard_ke"] = payload["cl_standard_kk"] * spin_ratio
        payload["covariance_ratio_pw2_over_ppww"] = payload["cl_pw"]**2 / (payload["cl_pp"] * payload["cl_ww"])
        payload.update({"kernel_" + key: value for key, value in kernels.items()})
        hashes = dict(standard_metadata["source_hashes"])
        ccl_root = Path(ccl.__file__).resolve().parent
        for path in [Path(__file__), Path(folded.__file__),
                     folded.CANOES / "src/canoes/nuell/transfer/basis.py",
                     folded.CANOES / "src/canoes/nuell/transfer/_sachs_decomp.py",
                     folded.CANOES / "src/canoes/nuell/cls/_compute_cl_windowed_sachs.py",
                     ccl_root / "tracers.py", ccl_root / "cells.py",
                     ccl_root / "nonlimber/_nonlimber_FKEM.py"]:
            hashes[str(path.resolve())] = fingerprint(path)
        metadata = {
            "python_executable": sys.executable, "python_version": sys.version,
            "pyccl_version": ccl.__version__, "pyccl_path": str(Path(ccl.__file__).resolve()),
            "source_fold": source_metadata,
            "standard_ccl": standard_metadata, "operator_mapping": mapping,
            "mapping_policy": "Each registered SachsCoefTerm projected independently by CCL, then all single-term cross-spectra summed. M[a] is retained and FKEM matter growth cancelled in the radial kernel. The shared Poisson amplitude enters once per transfer.",
            "multi_term_workaround": "Installed FKEM _chi_integrands has a get_transfer broadcasting error for multi-term Tracers. The exact bilinear pair expansion uses only supported single-term API calls, without package edits.",
            "boundary_policy": f"Physical chi_min cut with {args.kernel_edge_width} Mpc transition and radial zero padding. CCL Bessel derivatives act on the field. No integration-by-parts boundary omission.",
            "comparison_scope": "Five-multipole same-operator numerical diagnostic, not a converged figure or acceptance of the full Ricci t-form.",
            "settings": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            "channel_seconds": timings, "elapsed_seconds": time.monotonic() - begin,
            "source_hashes": hashes,
            "warnings": [{"category": item.category.__name__, "message": str(item.message)} for item in caught],
        }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, **payload, metadata_json=json.dumps(metadata, sort_keys=True))
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"out": str(args.out), "seconds": metadata["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
