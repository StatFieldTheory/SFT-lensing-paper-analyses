"""Assemble a deployed-format kappa3 vertex table from cached pieces.

Takes one LOW piece and a set of HIGH band pieces, sums the bands whose
upper edge does not exceed the requested cutoff, optionally restricts the
result to a subset of rows, and writes a table in exactly the layout the
production callable reads.

Because rows are independent and the bands are disjoint, every table this
produces is bit-identical to what a direct single-shot build of the same
(rows, cutoff, model) would have produced, given the same quadrature
windows. The one thing it cannot fake is a different window partition: a
single (60, 1000] rule with n_ell = 96 is a different quadrature from four
octave rules covering the same range, and both are legitimate. Pass the
band set you mean.

Usage::

    python assemble.py --low ../../products/pieces/low_r4.npz \
        --bands '../../products/pieces/band_tree_r4_*.npz' \
        --cutoff 1000 --subset-npz ../../products/triples_full_r1.npz \
        --out ../../products/table_tree_r1_cut1000.npz
"""

from __future__ import annotations

import argparse
import glob
from dataclasses import replace
from pathlib import Path

import numpy as np

from _common import CHANNELS, bootstrap_paths, dump_meta, load_meta

MOD_KEYS = ("bmod", "dmod")


def _as_output(piece, cls):
    """Rebuild a Kappa3Output from a saved piece."""
    return cls(
        cosine_triples=piece["cosine_triples"],
        lambda_shells=piece["lambda_shells"],
        chi_shells=piece["chi_shells"],
        z_shells=piece["z_shells"],
        zeta_TTT=piece["zeta_TTT"], zeta_TTP=piece["zeta_TTP"],
        zeta_TPP=piece["zeta_TPP"], zeta_PPP=piece["zeta_PPP"],
        ell_max=int(piece["ell_max"]),
        cosmo_meta=load_meta(piece["cosmo_meta"]),
    )


def _row_subset(built: np.ndarray, wanted: np.ndarray) -> np.ndarray:
    """Indices of `wanted` rows inside `built`, failing loudly if any is absent."""
    key = {tuple(np.round(r, 12)): i for i, r in enumerate(built)}
    try:
        return np.asarray([key[tuple(np.round(r, 12))] for r in wanted], dtype=int)
    except KeyError as exc:
        raise SystemExit(
            f"subset row {exc.args[0]} is not in the built table; the subset "
            "must be a subset of the rows the pieces were built on") from exc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--low", type=Path, required=True)
    ap.add_argument("--bands", nargs="+", required=True,
                    help="band NPZ paths or globs")
    ap.add_argument("--cutoff", type=int, required=True,
                    help="keep bands with upper edge <= this")
    ap.add_argument("--subset-npz", type=Path, default=None,
                    help="restrict output rows to this triple set")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    bootstrap_paths()
    from canoes.sachs import kappa3_combine_low_high
    from canoes.sachs.kappa3 import Kappa3Output

    paths: list[Path] = []
    for pattern in args.bands:
        hits = sorted(glob.glob(pattern))
        paths.extend(Path(h) for h in (hits if hits else [pattern]))

    low = np.load(args.low, allow_pickle=False)
    fingerprint = load_meta(low["fingerprint"])

    kept, models, quadrature = [], set(), set()
    for path in paths:
        piece = np.load(path, allow_pickle=False)
        build = load_meta(piece["build"])
        if load_meta(piece["fingerprint"]) != fingerprint:
            raise SystemExit(f"{path.name}: built on a different grid than the LOW piece")
        if int(build["hi"]) > args.cutoff:
            continue
        kept.append((int(build["lo"]), int(build["hi"]), piece, build))
        models.add(build["b_model"])
        quadrature.add((int(build["n_ell"]), int(build["n_phi"])))
    if not kept:
        raise SystemExit(f"no band with upper edge <= {args.cutoff}")
    if len(models) != 1:
        raise SystemExit(f"bands mix bispectrum models: {sorted(models)}")
    kept.sort()

    # The kept windows must tile (ell_cut, cutoff] exactly: a gap silently
    # drops multipoles, an overlap silently double counts them.
    edge = int(load_meta(low["build"])["ell_cut"])
    for lo, hi, _, _ in kept:
        if lo != edge:
            raise SystemExit(f"bands do not tile: expected a window starting at {edge}, got ({lo}, {hi}]")
        edge = hi
    if edge != args.cutoff:
        raise SystemExit(f"bands stop at {edge}, not at the requested cutoff {args.cutoff}")

    summed = {f"zeta_{c}": sum(p[f"zeta_{c}"] for _, _, p, _ in kept) for c in CHANNELS}
    summed.update({k: sum(p[k] for _, _, p, _ in kept) for k in MOD_KEYS})
    first = kept[0][2]

    high = Kappa3Output(
        cosine_triples=first["cosine_triples"],
        lambda_shells=first["lambda_shells"], chi_shells=first["chi_shells"],
        z_shells=first["z_shells"], ell_max=int(args.cutoff),
        cosmo_meta=load_meta(first["cosmo_meta"]), **summed_channels(summed))
    combined = kappa3_combine_low_high(_as_output(low, Kappa3Output), high)

    zeta_Bmod = low["bmod"] + summed["bmod"]
    zeta_Dmod = low["dmod"] + summed["dmod"]

    rows = slice(None)
    triples_out = np.asarray(combined.cosine_triples)
    if args.subset_npz is not None:
        wanted = np.asarray(np.load(args.subset_npz)["cosine_triples"], float)
        rows = _row_subset(triples_out, wanted)
        triples_out = triples_out[rows]
        combined = replace(combined, cosine_triples=triples_out, **{
            f"zeta_{c}": np.asarray(getattr(combined, f"zeta_{c}"))[rows]
            for c in CHANNELS})
        zeta_Bmod, zeta_Dmod = zeta_Bmod[rows], zeta_Dmod[rows]

    meta = dict(combined.cosmo_meta)
    meta.update(
        producer="driver_field_emulators/code/rebuild/assemble.py",
        artifact="equal_time_limber_lambda_NON_R_contracted",
        equal_time=True, limber_high=True,
        already_R_convolved=False, already_R_contracted=False,
        radial_measure=load_meta(low["build"])["radial_measure"],
        radial_density_measure=load_meta(low["build"])["radial_measure"],
        radial_kernel_convention="lambda_three_leg_native_density",
        sachs_prefactor_per_leg="(1+z)^4 bare driving field (SACHS_GEOMETRIC_POWER=4)",
        poisson_policy="P_matter input; tree_phi applies A(a) per leg internally",
        query_coordinate="lambda_project_mpc",
        ell_cut=int(load_meta(low["build"])["ell_cut"]),
        ell_high_max=int(args.cutoff),
        b_delta_model=sorted(models)[0],
        band_windows=[(lo, hi) for lo, hi, _, _ in kept],
        band_quadrature=sorted(quadrature),
        n_rows=int(triples_out.shape[0]),
        assembled_from=dict(low=args.low.name,
                            bands=[f"({lo},{hi}]" for lo, hi, _, _ in kept]),
        has_modulus_channels=True,
    )
    combined = replace(combined, cosmo_meta=meta)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    combined.save_npz(args.out)
    existing = dict(np.load(args.out, allow_pickle=True))
    np.savez(args.out, **existing, zeta_Bmod=zeta_Bmod, zeta_Dmod=zeta_Dmod,
             has_modulus_channels_flag=np.bool_(True))
    print(f"[asm] {triples_out.shape[0]} rows, model={sorted(models)[0]}, "
          f"cutoff={args.cutoff}, windows={[(lo, hi) for lo, hi, _, _ in kept]}")
    print(f"[asm] -> {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    return 0


def summed_channels(summed: dict) -> dict:
    return {c: summed[c] for c in summed if c.startswith("zeta_")}


if __name__ == "__main__":
    raise SystemExit(main())
