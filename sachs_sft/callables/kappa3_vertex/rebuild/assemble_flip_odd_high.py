"""Assemble a kappa3 vertex table with a checked spin-2 sign convention.

canoes before 810b133 built the HIGH (Limber) branch with
``response_psi = -response_phi``. In the exp(2 i phi_l) basis the HIGH branch
works in, Psi0(l) = +exp(2 i phi_l) Phi00(l), so those bands carry zeta_TTP,
zeta_PPP and dmod (one or three spin-2 legs) with the wrong sign; zeta_TTT,
zeta_TPP and bmod are right. Bands built with canoes 810b133 or later carry the
correct sign. The two kinds must never be summed together, and ``assemble.py``
cannot tell them apart.

This script repeats the summation of ``assemble.py`` in the same order and reads
the convention of every band from ``cosmo_meta["potential_convention"]``:

* ``--flip-odd-high``: every band must be of the old convention; the summed
  HIGH ``zeta_TTP``, ``zeta_PPP`` and ``dmod`` are multiplied by -1 before LOW
  and HIGH are combined. LOW arrays and the even channels are untouched.
  This is the control table of paper task T-001 (sign fix only, stored
  quadrature).
* without the flag: every band must be of the new convention; the sum is the
  plain ``assemble.py`` sum.

A band whose ``potential_convention`` is absent or unrecognised stops the run.
The LOW piece (exact branch) is not tested; it has no such key. The output
keeps the layout and metadata of ``assemble.py`` and adds
``potential_convention`` and ``spin2_sign_convention``. An existing ``--out``
is refused.

Usage::

    python assemble_flip_odd_high.py --flip-odd-high \
        --low products/pieces_nphi512/low_permclosed_r4_np512.npz \
        --bands 'products/pieces_nphi512/band_tree_permclosed_*_np512.npz' \
        --cutoff 15360 --out products/pieces_nphi512/table_permclosed_np512_spin2fix.npz
"""

from __future__ import annotations

import argparse
import glob
from dataclasses import replace
from pathlib import Path

import numpy as np

from _common import (CHANNELS, bootstrap_paths, build_leg_order, build_pk_identity,
                     load_meta)
from assemble import MOD_KEYS, _as_output, _row_subset, summed_channels

OLD_MARK = "Psi0=-A(a)"
NEW_MARK = "Psi0=+A(a)"
#: HIGH arrays with an odd number of spin-2 legs.
ODD_HIGH = ("zeta_TTP", "zeta_PPP", "dmod")
POTENTIAL_CONVENTION = (
    "Phi=Psi, W=2*Phi; high-ell Sachs responses are "
    "Phi00=A(a)*(1+z)^4*delta and Psi0=+A(a)*(1+z)^4*delta "
    "times explicit flat-sky spin phases"
)


def band_convention(name: str, cosmo_meta: dict) -> str:
    """'old' or 'new' from a band's potential_convention; stop on anything else."""
    text = cosmo_meta.get("potential_convention")
    if not isinstance(text, str):
        raise SystemExit(f"{name}: no potential_convention in cosmo_meta")
    has_old, has_new = OLD_MARK in text, NEW_MARK in text
    if has_old == has_new:
        raise SystemExit(f"{name}: unrecognised potential_convention {text!r}")
    return "old" if has_old else "new"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--low", type=Path, required=True)
    ap.add_argument("--bands", nargs="+", required=True, help="band NPZ paths or globs")
    ap.add_argument("--cutoff", type=int, required=True,
                    help="keep bands with upper edge <= this")
    ap.add_argument("--subset-npz", type=Path, default=None,
                    help="restrict output rows to this triple set")
    ap.add_argument("--flip-odd-high", action="store_true",
                    help="bands are of the old convention; reverse their odd channels")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit(f"{args.out} exists; refusing to overwrite")

    bootstrap_paths()
    from canoes.sachs import kappa3_combine_low_high
    from canoes.sachs.kappa3 import Kappa3Output

    paths: list[Path] = []
    for pattern in args.bands:
        hits = sorted(glob.glob(pattern))
        paths.extend(Path(h) for h in (hits if hits else [pattern]))

    low = np.load(args.low, allow_pickle=False)
    fingerprint = load_meta(low["fingerprint"])
    low_table, low_sha = build_pk_identity(load_meta(low["build"]))
    wanted_convention = "old" if args.flip_odd_high else "new"

    kept, models, quadrature, pk_identities, leg_orders = [], set(), set(), set(), set()
    for path in paths:
        piece = np.load(path, allow_pickle=False)
        build = load_meta(piece["build"])
        if load_meta(piece["fingerprint"]) != fingerprint:
            raise SystemExit(f"{path.name}: built on a different grid than the LOW piece")
        if int(build["hi"]) > args.cutoff:
            continue
        convention = band_convention(path.name, load_meta(piece["cosmo_meta"]))
        if convention != wanted_convention:
            raise SystemExit(
                f"{path.name}: {convention} spin-2 convention, expected {wanted_convention} "
                f"({'with' if args.flip_odd_high else 'without'} --flip-odd-high)")
        kept.append((int(build["lo"]), int(build["hi"]), piece, build))
        models.add(build["b_model"])
        quadrature.add((int(build["n_ell"]), int(build["n_phi"])))
        pk_identities.add(build_pk_identity(build))
        leg_orders.add(build_leg_order(build))
    if not kept:
        raise SystemExit(f"no band with upper edge <= {args.cutoff}")
    if len(models) != 1:
        raise SystemExit(f"bands mix bispectrum models: {sorted(models)}")
    if len(pk_identities) != 1:
        raise SystemExit(f"bands mix P(k) table identities: {sorted(pk_identities, key=repr)}")
    high_table, high_sha = next(iter(pk_identities))
    if len(leg_orders) != 1:
        raise SystemExit(f"bands mix leg orders: {sorted(leg_orders)}")
    kept.sort()

    edge = int(load_meta(low["build"])["ell_cut"])
    for lo, hi, _, _ in kept:
        if lo != edge:
            raise SystemExit(f"bands do not tile: expected a window starting at {edge}, got ({lo}, {hi}]")
        edge = hi
    if edge != args.cutoff:
        raise SystemExit(f"bands stop at {edge}, not at the requested cutoff {args.cutoff}")

    summed = {f"zeta_{c}": sum(p[f"zeta_{c}"] for _, _, p, _ in kept) for c in CHANNELS}
    summed.update({k: sum(p[k] for _, _, p, _ in kept) for k in MOD_KEYS})
    if args.flip_odd_high:
        summed.update({k: -summed[k] for k in ODD_HIGH})
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
        producer="SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/rebuild/assemble_flip_odd_high.py",
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
        pk_table_high=high_table,
        pk_table_low=low_table,
        leg_order_high=sorted(leg_orders)[0],
        band_windows=[(lo, hi) for lo, hi, _, _ in kept],
        band_quadrature=sorted(quadrature),
        n_rows=int(triples_out.shape[0]),
        assembled_from=dict(low=args.low.name,
                            bands=[f"({lo},{hi}]" for lo, hi, _, _ in kept]),
        has_modulus_channels=True,
        potential_convention=POTENTIAL_CONVENTION,
        spin2_sign_convention=(
            "plus_exp2iphi (flipped from stored HIGH pieces, 2026-10)" if args.flip_odd_high
            else "plus_exp2iphi (HIGH pieces built with canoes 810b133 or later)"),
    )
    for branch, digest in (("high", high_sha), ("low", low_sha)):
        key = f"pk_table_{branch}_sha256"
        meta.pop(key, None)
        if digest is not None:
            meta[key] = digest
    combined = replace(combined, cosmo_meta=meta)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    combined.save_npz(args.out)
    # re-read of the file written on the line above (object metadata), as in assemble.py
    existing = dict(np.load(args.out, allow_pickle=True))
    np.savez(args.out, **existing, zeta_Bmod=zeta_Bmod, zeta_Dmod=zeta_Dmod,
             has_modulus_channels_flag=np.bool_(True))
    print(f"[asm] {triples_out.shape[0]} rows, model={sorted(models)[0]}, cutoff={args.cutoff}, "
          f"convention={wanted_convention}, flip_odd_high={args.flip_odd_high}")
    print(f"[asm] -> {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
