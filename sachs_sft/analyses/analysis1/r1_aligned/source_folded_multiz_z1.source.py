"""Compute the Gaussian Sachs projection after folding each source leg first.

The source integral is the existing manuscript's effective affine kernel
(sections/cosmology.tex, labels `eq: K affine kernel` and `eq: xi kappa folded`).
It is evaluated before the angular covariance, avoiding interpolation of the
21-source C table. This is a reordering of linear integrals, with no amplitude
calibration. Fresh coefficient registrations supply the external weight to
canoes without altering its installed source or the physical lightcone.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import sys
import time
import warnings
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.interpolate import CubicSpline
from scipy.special import eval_legendre

HERE = Path(__file__).resolve().parent
SACHS = HERE.parents[2]
CANOES = Path("/Users/zzhang/projects/angular_statistics/canoes")
OMEGA_M = 0.3160919980475834
H = 0.6711
DEFAULT_SOURCE_MAP = SACHS / "analyses/r1_sft061/source_map.json"


def fingerprint(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def setup(args):
    sys.path.insert(0, str(CANOES / "src"))
    from canoes.cosmo.lightcone import build_lightcone, chi_of_lambda, z_of_chi
    from canoes.sachs.sft_input.corr_op.build import _load_pk_phi
    from canoes.nuell.transfer._sachs_decomp import (
        SachsCoefTerm, get_sachs_decomposition, register_sachs_decomposition,
    )
    from canoes.nuell.transfer.sachs_driver import COEFS_PHI00, COEFS_PSI0

    source_path = Path(getattr(args, "source_map", DEFAULT_SOURCE_MAP))
    mapping = json.loads(source_path.read_text())
    redshift = float(getattr(args, "z_source", 5.0))
    matches = np.flatnonzero(np.isclose(mapping["z"], redshift, rtol=0, atol=1e-10))
    if matches.size != 1:
        raise ValueError(f"Source redshift must identify one source_map entry: {mapping['z']}")
    mapping_cosmo = mapping["mapping_provenance"]["cosmology"]
    if not np.isclose(mapping_cosmo["Omega_m"], OMEGA_M, rtol=0, atol=1e-14) or not np.isclose(mapping_cosmo["h"], H, rtol=0, atol=1e-14):
        raise ValueError("The source map uses a different cosmology")
    index = int(matches[0])
    chi_s = float(mapping["chi_mpc"][index])
    lambda_s = float(mapping["lambda_mpc"][index])
    cosmo = SimpleNamespace(Omega_m=OMEGA_M, h=H)
    recovered_z = float(z_of_chi(cosmo, np.array([chi_s]), distance_unit="Mpc", z_max=1200)[0])
    recovered_chi = float(chi_of_lambda(cosmo, np.array([lambda_s]), convention="project", distance_unit="Mpc", z_max=1200)[0])
    if abs(recovered_z - redshift) > 1e-8 or abs(recovered_chi - chi_s) > 1e-7:
        raise ValueError("The source map does not match the live public canoes coordinate maps")
    if not 0 < args.chi_min < chi_s:
        raise ValueError("chi_min must lie strictly between the observer and selected source")
    source = {"z_source_requested": redshift, "z_source_recovered": recovered_z,
              "chi_source_mpc": chi_s, "lambda_source_mpc": lambda_s,
              "mapping_provenance": mapping["mapping_provenance"],
              "chi_roundtrip_abs_error_mpc": abs(recovered_chi - chi_s)}
    chi = np.geomspace(args.chi_min, chi_s, args.n_chi)
    lightcone, opz = build_lightcone(cosmo, chi, z_max=1200, distance_unit="Mpc")
    d_s = chi_s / float(opz[-1])

    # This is the continuous outer source integral in the cited affine kernel.
    # The integrand uses the same a and Dbar definitions as the existing C table.
    dense_chi = np.geomspace(args.chi_min, chi_s, args.source_nchi)
    dense_z = z_of_chi(cosmo, dense_chi, distance_unit="Mpc", z_max=1200)
    a = 1.0 / (1.0 + dense_z)
    d_bar = a * dense_chi
    source_integrand = a**2 / d_bar**2
    # Integrate from the source backwards to avoid subtracting nearby totals.
    tail = -cumulative_trapezoid(source_integrand[::-1], dense_chi[::-1], initial=0)[::-1]
    tail_spline = CubicSpline(dense_chi, tail)

    def external_weight(query):
        q = np.asarray(query, dtype=float)
        value = d_s**2 * tail_spline(np.clip(q, args.chi_min, chi_s))
        return np.where((q >= args.chi_min) & (q <= chi_s), value, 0.0)

    # Check the numerical source integral against independent adaptive quadrature.
    checks = []
    for lower in np.geomspace(args.chi_min, 0.99 * chi_s, 7):
        def integrand(distance):
            z = float(z_of_chi(cosmo, np.asarray(distance), distance_unit="Mpc", z_max=1200))
            scale = 1.0 / (1.0 + z)
            return scale**2 / (scale * distance)**2

        direct, error = quad(integrand, float(lower), chi_s, epsabs=1e-13, epsrel=2e-11)
        actual = float(tail_spline(lower))
        checks.append({"lower_chi_mpc": float(lower), "relative_error": actual / direct - 1.0,
                       "quad_error": error})
    if max(abs(row["relative_error"]) for row in checks) > 2e-7:
        raise RuntimeError("Continuous source integration failed its quadrature check")

    def register_folded(original, name, transverse=False):
        original_terms = get_sachs_decomposition(original)
        if transverse:
            # The angular-Laplacian term already present in PHI00's decomposition.
            original_terms = ((original_terms[0][0],), (), (), (), (), ())
            original = (lambda l2, x, hh, hp: l2 / x**2,) + tuple(
                lambda l2, x, hh, hp: 0.0 for _ in range(5)
            )
        folded = tuple(
            (lambda l2, x, hh, hp, fn=fn: fn(l2, x, hh, hp) * float(external_weight(x)))
            for fn in original
        )
        terms = tuple(tuple(
            SachsCoefTerm(
                label=f"{name}_{term.label}",
                ell_prefactor_fn=term.ell_prefactor_fn,
                chi_kernel_fn=lambda x, lc, fn=term.chi_kernel_fn: fn(x, lc) * external_weight(x),
            ) for term in basis) for basis in original_terms)
        register_sachs_decomposition(folded, terms)
        return folded

    vectors = {
        "phi": register_folded(COEFS_PHI00, "source_folded_phi"),
        "weyl": register_folded(COEFS_PSI0, "source_folded_weyl"),
        "transverse": register_folded(COEFS_PHI00, "source_folded_transverse", transverse=True),
    }
    power_path = CANOES / "examples/data/PCAMB_pyccl_stf_fid_z0.txt"
    pk = _load_pk_phi(power_path, omega_m=OMEGA_M, h_param=H, n_s=0.97, distance_unit="Mpc")
    metadata = {
        "source": source, "source_metadata_sha256": fingerprint(source_path),
        "source_map": str(source_path.resolve()), "source_redshift": redshift,
        "source_chi_mpc": chi_s, "source_lambda_mpc": lambda_s,
        "physical_chi_min_mpc": args.chi_min,
        "source_fold": "Numerical continuous outer affine integral before angular covariance, manuscript K affine kernel",
        "source_floor": "Observer affine zero, no artificial 406 Mpc source floor",
        "registration": "Fresh CoefVec plus SachsCoefTerm external weights, original lightcone M/calH/calHp unchanged",
        "derivative_order": "External window multiplies each field-derivative term. The operator acts on the field, not on the window. Existing low-ell IBP differentiates the complete weighted integrand.",
        "source_quadrature_check": checks,
        "power_table_sha256": fingerprint(power_path),
        "power_table": str(power_path),
        "script_sha256": fingerprint(Path(__file__)),
        "operator_policy": "Full registered Ricci and Weyl below split. ID-basis Limber above split. Transverse Ricci-only output is also saved separately.",
        "low_ell_boundary_caveat": "Existing low-ell canoes IBP assumes negligible observer-boundary terms. The finite 50 Mpc inner cutoff needs a separate sensitivity check for Ricci derivatives.",
    }
    return chi_s, chi, lightcone, opz, vectors, pk, metadata


def high_ell(ell, ca, cb, chi_s, chi, lc, opz, pk, shift):
    """Reuse canoes ID windows and its high-ell single-distance quadrature."""
    from canoes.nuell.transfer._sachs_decomp import get_sachs_decomposition
    from canoes.sachs.sft_input.corr_op._limber import _id_basis_windows, _safe_pk

    def windows(coefs):
        return _id_basis_windows(
            decomp=get_sachs_decomposition(coefs), chi=chi, lam=np.array([chi_s]),
            opz_chi=opz, opz_lam=opz[-1:], m_id=lc.M[0], lightcone=lc,
            l2_arr=ell.astype(float) * (ell + 1.0),
        )
    ga, pa = windows(ca)
    gb, pb = windows(cb)
    out = np.empty(ell.size)
    for index, multipole in enumerate(ell):
        wa = np.einsum("t,tc->c", pa[:, index], ga[:, 0, :])
        wb = np.einsum("t,tc->c", pb[:, index], gb[:, 0, :])
        out[index] = np.trapezoid(_safe_pk(pk, (multipole + shift) / chi) * wa * wb / chi**2, chi)
    return out


def multipole_grid(choice, split, sparse_nodes):
    if choice == "full":
        return np.arange(2, 5001)
    if choice in {"sparse", "withheld"}:
        nodes = np.unique(np.rint(np.geomspace(31, split, sparse_nodes)).astype(int))
        if choice == "withheld":
            return np.unique(np.rint(np.sqrt(nodes[:-1] * nodes[1:])).astype(int))
        return np.unique(np.concatenate((np.arange(2, 31), nodes,
                                         np.arange(split + 1, 5001))))
    return np.unique([int(value) for value in choice.split(",")])


def cached_direct_module(include_observer_boundary=False):
    """Load an isolated helper namespace and memoize exact Bessel arguments.

    The installed module and scipy functions are untouched. This changes only
    repeated evaluation cost, retaining the helper's original quadrature and
    derivative operations. Arrays are keyed by their full binary content.
    """
    path = CANOES / "src/canoes/nuell/cls/_compute_cl_windowed_sachs_low_ell.py"
    spec = importlib.util.spec_from_file_location("_r1_direct_sachs_snapshot", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original = module.spherical_jn
    cache = OrderedDict()
    stats = {"hits": 0, "misses": 0, "bytes": 0, "max_bytes": 512 * 1024**2,
             "helper_source": str(path), "helper_sha256": fingerprint(path)}

    def memoized(order, arguments, derivative=False):
        array = np.ascontiguousarray(arguments)
        digest = hashlib.sha256(memoryview(array).cast("B")).digest()
        key = (int(order), bool(derivative), array.shape, array.dtype.str, digest)
        if key in cache:
            stats["hits"] += 1
            cache.move_to_end(key)
            return cache[key]
        stats["misses"] += 1
        value = np.asarray(original(order, array, derivative=derivative))
        value.setflags(write=False)
        while cache and stats["bytes"] + value.nbytes > stats["max_bytes"]:
            _, old = cache.popitem(last=False)
            stats["bytes"] -= old.nbytes
        if value.nbytes <= stats["max_bytes"]:
            cache[key] = value
            stats["bytes"] += value.nbytes
        return value

    module.spherical_jn = memoized
    if include_observer_boundary:
        original_transfer = module._build_full_basis_per_leg_factor

        def transfer_with_lower_boundary(**kwargs):
            """Evaluate the cited helper's IBP boundary at both radial endpoints."""
            result = original_transfer(**kwargs)
            chi = kwargs["chi"]
            lam = kwargs["lam"]
            opz = kwargs["opz_chi"]
            amp = ((chi[None, :] / opz[None, :]) /
                   (lam[:, None] / kwargs["opz_lam"][:, None]))**2 * opz[None, :]**2
            k = kwargs["k_nodes"]
            for basis, (time_order, radial_order) in enumerate(module.BASIS_TO_AB):
                if radial_order == 0:
                    continue
                for term in kwargs["decomp"][basis]:
                    g = (np.asarray(term.chi_kernel_fn(chi, kwargs["lightcone"]))[None, :]
                         * kwargs["M_table"][time_order][None, :] * amp)
                    derivative = np.gradient(g, chi, axis=1)[:, 0] if radial_order == 2 else None
                    for index, ell in enumerate(kwargs["ell_arr"]):
                        angular = float(term.ell_prefactor_fn(kwargs["L2_arr"][index]))
                        j = memoized(int(ell), k * chi[0])
                        if radial_order == 1:
                            boundary = -j[:, None] * g[:, 0][None, :]
                        else:
                            j_prime = memoized(int(ell), k * chi[0], derivative=True)
                            boundary = (-k[:, None] * j_prime[:, None] * g[:, 0][None, :]
                                        + j[:, None] * derivative[None, :])
                        result[index] += angular * boundary
            return result

        module._build_full_basis_per_leg_factor = transfer_with_lower_boundary
    stats["observer_boundary"] = "included" if include_observer_boundary else "existing helper omission"
    return module, stats


def complete_angular_products(payload, vectors, args):
    """Transform a complete or explicitly interpolated angular spectrum."""
    from canoes.nuell.transfer._sachs_decomp import get_sachs_decomposition
    from ccl_reference import transform

    ell = payload["ell"]
    derived_transverse = {"cl_tt", "cl_tw"}.difference(payload)
    if "cl_ww" in payload:
        l2 = ell.astype(float) * (ell + 1)
        transverse = get_sachs_decomposition(vectors["transverse"])[0][0]
        weyl = get_sachs_decomposition(vectors["weyl"])[0][0]
        ratio = np.asarray([transverse.ell_prefactor_fn(x) / weyl.ell_prefactor_fn(x)
                            for x in l2])
        # Identical registered radial windows permit exact angular-factor reuse.
        payload.setdefault("cl_tt", payload["cl_ww"] * ratio**2)
        payload.setdefault("cl_tw", payload["cl_ww"] * ratio)
    if args.ells == "sparse":
        full_ell = np.arange(2, 5001)
        payload["sampled_ell"] = ell.copy()
        for key in list(payload):
            if not key.startswith("cl_"):
                continue
            values = payload[key]
            payload["sampled_" + key] = values.copy()
            low = ell <= args.split
            sign = float(np.sign(values[0]))
            if np.any(values * sign <= 0):
                raise ValueError(f"Log interpolation requires a fixed nonzero sign: {key}")
            curve = CubicSpline(np.log(ell[low]), np.log(np.abs(values[low])))
            dense = np.empty(full_ell.size)
            dense_low = full_ell <= args.split
            dense[dense_low] = sign * np.exp(curve(np.log(full_ell[dense_low])))
            dense[~dense_low] = values[~low]
            payload[key] = dense
        payload["ell"] = full_ell
        ell = full_ell
        if "cl_ww" in payload:
            full_l2 = ell.astype(float) * (ell + 1)
            full_ratio = np.asarray([transverse.ell_prefactor_fn(x) / weyl.ell_prefactor_fn(x)
                                     for x in full_l2])
            # Preserve the exact angular relation between, not only at, nodes.
            if "cl_tt" in derived_transverse:
                payload["cl_tt"] = payload["cl_ww"] * full_ratio**2
            if "cl_tw" in derived_transverse:
                payload["cl_tw"] = payload["cl_ww"] * full_ratio
    if args.ells not in {"full", "sparse"}:
        return
    angle_path = HERE / "outputs/order0_table_integral/order0_observables.npz"
    with np.load(angle_path, allow_pickle=False) as source:
        gamma = np.asarray(source["theta_arcmin"])
    if set(args.channels.split(",")) == {"cl_pp"}:
        # Same full-sky scalar sum and taper as ccl_reference.transform,
        # without evaluating unused spin-two angular kernels.
        taper = np.ones_like(ell, dtype=float)
        taper_start = float(ell[0] + 0.8 * (ell[-1] - ell[0]))
        tapered = ell > taper_start
        taper[tapered] = 0.5 * (1.0 + np.cos(np.pi * (ell[tapered] - taper_start) / (ell[-1] - taper_start)))
        weight = (2.0 * ell + 1.0) * taper / (4.0 * np.pi)
        kk = np.asarray([np.sum(weight * payload["cl_pp"] * eval_legendre(ell, np.cos(angle)))
                         for angle in np.deg2rad(gamma / 60.0)])
        payload.update(theta_arcmin=gamma, gamma_arcmin=gamma, kk=kk,
                       xi_kappa=kk, cl_kk=payload["cl_pp"], ell_taper=taper)
        return
    if args.observable_operator == "full":
        convergence, cross = payload["cl_pp"], -payload["cl_pw"]
    else:
        convergence, cross = payload["cl_tt"], -payload["cl_tw"]
    correlations, taper = transform(ell, convergence, payload["cl_ww"], cross, gamma, 0.2)
    payload.update(correlations, theta_arcmin=gamma, gamma_arcmin=gamma,
                   kk=correlations["xi_kappa"], xip=correlations["xi_plus"],
                   xim=correlations["xi_minus"], kgt=correlations["xi_kappa_gamma"],
                   cl_kk=convergence, cl_ee=payload["cl_ww"], cl_ke=cross,
                   ell_taper=taper)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--z-source", type=float, default=5.0,
                        help="Source redshift selected from the validated source map")
    parser.add_argument("--source-map", type=Path, default=DEFAULT_SOURCE_MAP)
    parser.add_argument("--ells", default="2,10,100,500,800",
                        help="Comma-separated multipoles, full, sparse, or withheld")
    parser.add_argument("--n-chi", type=int, default=512)
    parser.add_argument("--source-nchi", type=int, default=100001)
    parser.add_argument("--chi-min", type=float, default=50.0)
    parser.add_argument("--n-fftlog", type=int, default=256)
    parser.add_argument("--n-t", type=int, default=32)
    parser.add_argument("--low-ell-n-k", type=int, default=256)
    parser.add_argument("--low-ell-threshold", type=int, default=5)
    parser.add_argument("--kmin", type=float, default=1e-4)
    parser.add_argument("--kmax", type=float, default=50.0)
    parser.add_argument("--low-ell-kmin", type=float, default=1e-4)
    parser.add_argument("--low-ell-kmax", type=float, default=10.0)
    parser.add_argument("--split", type=int, default=800)
    parser.add_argument("--ell-block", type=int, default=32)
    parser.add_argument("--limber-shift", type=float, default=0.5)
    parser.add_argument("--channels", default="cl_pp,cl_pw,cl_ww,cl_tt,cl_tw",
                        help="Comma-separated output channels for bounded convergence probes")
    parser.add_argument("--sparse-nodes", type=int, default=49)
    parser.add_argument("--observable-operator", choices=("full", "transverse"), default="full")
    parser.add_argument("--cached-direct", action="store_true",
                        help="Use isolated exact-content Bessel caching in the existing direct-k helper")
    parser.add_argument("--observer-boundary", action="store_true",
                        help="Retain the existing helper's IBP lower-end terms at finite chi_min")
    args = parser.parse_args()
    if args.observer_boundary and not args.cached_direct:
        parser.error("--observer-boundary requires --cached-direct")
    if args.ell_block < 1 or args.n_chi < 3 or args.low_ell_n_k < 2 or args.sparse_nodes < 4:
        parser.error("Grid sizes and block size are too small")
    if not np.isfinite(args.z_source) or args.z_source <= 0:
        parser.error("--z-source must be finite and positive")
    channel_keys = set(args.channels.split(","))
    if channel_keys == {"cl_pp"} and args.observable_operator != "full":
        parser.error("The cl_pp-only transform requires --observable-operator full")
    if args.ells in {"full", "sparse"} and channel_keys != {"cl_pp"}:
        required = {"cl_pp", "cl_pw", "cl_ww"} if args.observable_operator == "full" else {"cl_ww"}
        if not required.issubset(channel_keys):
            parser.error(f"The selected angular transform requires channels {sorted(required)}")
    if args.out.suffix != ".npz":
        parser.error("--out must end in .npz")
    if any(args.out.with_suffix(suffix).exists() for suffix in (".npz", ".json", ".source.py")):
        parser.error("Output exists. Choose a new filename to preserve earlier artifacts.")
    script_bytes = Path(__file__).read_bytes()
    start = time.monotonic()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        chi_s, chi, lc, opz, vectors, pk, metadata = setup(args)
        from canoes.nuell.cls._compute_cl_windowed_sachs import compute_cl_windowed_sachs
        ell = multipole_grid(args.ells, args.split, args.sparse_nodes)
        channels = {"cl_pp": ("phi", "phi"), "cl_pw": ("phi", "weyl"),
                    "cl_ww": ("weyl", "weyl"), "cl_tt": ("transverse", "transverse"),
                    "cl_tw": ("transverse", "weyl")}
        requested = args.channels.split(",")
        if not requested or any(key not in channels for key in requested):
            parser.error("Unknown or empty channel list")
        channels = {key: channels[key] for key in requested}
        payload = {"ell": ell, "chi_mpc": chi}
        timings = {key: 0.0 for key in channels}
        low = ell <= args.split
        direct_module, cache_stats = cached_direct_module(args.observer_boundary) if args.cached_direct else (None, None)
        for key in channels:
            payload[key] = np.empty(ell.size)
        indices = np.flatnonzero(low)
        # Keep all channels of a small multipole block together, so the exact
        # Bessel cache is reused across both operators while memory stays bounded.
        for offset in range(0, indices.size, args.ell_block):
            subset = indices[offset:offset + args.ell_block]
            for key, (left, right) in channels.items():
                begin = time.monotonic()
                values = payload[key]
                print(f"{key}: full operators ell={ell[subset].tolist()}", flush=True)
                if direct_module is not None:
                    values[subset] = direct_module.compute_cl_windowed_sachs_low_ell(
                        ell[subset], np.array([chi_s]), None,
                        coefs_a=vectors[left], coefs_b=vectors[right], lightcone=lc,
                        one_plus_z_chi=opz, pk=pk, atom_mode="full",
                        n_k=args.low_ell_n_k,
                        k_range=(args.low_ell_kmin, args.low_ell_kmax),
                    )[:, 0, 0]
                    timings[key] += time.monotonic() - begin
                    continue
                values[subset] = compute_cl_windowed_sachs(
                    ell_arr=ell[subset], lambda_grid_a=np.array([chi_s]),
                    coefs_a=vectors[left], coefs_b=vectors[right], lightcone=lc,
                    one_plus_z_chi=opz, pk=pk, atom_mode="full", n_t=args.n_t,
                    n_fftlog=args.n_fftlog, low_ell_n_k=args.low_ell_n_k,
                    ell_lo_threshold=args.low_ell_threshold,
                    kmin=args.kmin, kmax=args.kmax,
                    low_ell_k_range=(args.low_ell_kmin, args.low_ell_kmax),
                    apply_density_correction=True,
                )[:, 0, 0]
                timings[key] += time.monotonic() - begin
        for key, (left, right) in channels.items():
            begin = time.monotonic()
            values = payload[key]
            if np.any(~low):
                values[~low] = high_ell(ell[~low], vectors[left], vectors[right],
                                       chi_s, chi, lc, opz, pk, args.limber_shift)
            if not np.all(np.isfinite(values)):
                raise RuntimeError(f"Non-finite output in {key}")
            timings[key] += time.monotonic() - begin
            print(f"{key}: seconds={timings[key]:.3f}, values={values.tolist() if ell.size<20 else 'saved'}", flush=True)
        complete_angular_products(payload, vectors, args)
        metadata.update(settings={key: str(value) if isinstance(value, Path) else value
                                  for key, value in vars(args).items()},
                        interpolation="log-CubicSpline of fixed-sign spectra below split" if args.ells == "sparse" else "none",
                        observable_operator=args.observable_operator,
                        direct_bessel_cache=cache_stats,
                        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024),
                        elapsed_seconds=time.monotonic() - start, channel_seconds=timings,
                        warnings=[{"category": item.category.__name__, "message": str(item.message)} for item in caught])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, **payload, metadata_json=json.dumps(metadata, sort_keys=True))
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    args.out.with_suffix(".source.py").write_bytes(script_bytes)
    print(json.dumps({"out": str(args.out), "elapsed_seconds": metadata["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
