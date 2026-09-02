"""Vertex amplitude as a function of the multipole cutoff, and in k-space.

Reads the octave band pieces and forms partial sums. Because the windows
are disjoint in max-leg, the partial sum up to an octave edge IS the vertex
built with that cutoff, so one band set gives the whole curve without
rebuilding anything.

Two outputs, and the second is the one that turns a convergence statement
into a model-error statement:

* zeta(ell_max) at chosen separations and shells, for each bispectrum model
* the same contributions expressed in physical wavenumber. Under the Limber
  branch a multipole ell at a shell of comoving distance chi_h maps to
  k = ell / chi_h, and the hard pair's third leg reaches k = 2 ell / chi_h.
  Comparing that against the range where the nonlinear bispectrum fit is
  calibrated, and against the scale where baryons stop being a small
  correction, says how much of the converged answer rests on modes the
  model does not control.

Usage::

    python analyse_cutoff.py --pieces products/pieces \
        --out products/cutoff_study.npz
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np

from _common import load_meta

ARCMIN = np.pi / (180.0 * 60.0)

#: Separations to report, arcmin. Chosen from the production grid.
GAMMAS = (0.5, 5.30, 21.88, 114.27)
#: BiHalofit is fitted to simulations resolving up to roughly this wavenumber.
BIHALOFIT_CALIBRATION_KMAX = 10.0
#: Above roughly this wavenumber baryonic feedback is not a small correction.
BARYON_SCALE_K = 1.0


def load_bands(pieces: Path, model: str) -> list[dict]:
    """Octave bands for one model, as a chain that tiles without overlap."""
    out = []
    for path in sorted(glob.glob(str(pieces / f"band_{model}_r4_*.npz"))):
        piece = np.load(path, allow_pickle=False)
        build = load_meta(piece["build"])
        out.append(dict(lo=int(build["lo"]), hi=int(build["hi"]),
                        path=Path(path), data=piece))
    out.sort(key=lambda b: (b["lo"], b["hi"]))

    # The pieces directory also holds windows built for other purposes, such
    # as the (960, 1000] quadrature control, whose range is contained in an
    # octave. Summing those alongside the octaves double counts multipoles
    # silently. Keep the chain that tiles: at each edge take the band that
    # continues it, preferring the octave hi = 2 lo.
    chain, edge = [], None
    for band in out:
        if edge is None:
            edge = band["lo"]
        if band["lo"] != edge:
            continue
        same = [b for b in out if b["lo"] == edge]
        pick = next((b for b in same if b["hi"] == 2 * b["lo"]), same[-1])
        if pick not in chain:
            chain.append(pick)
            edge = pick["hi"]
    dropped = [f"({b['lo']},{b['hi']}]" for b in out if b not in chain]
    if dropped:
        print(f"[cut] {model}: not part of the octave chain, excluded: "
              f"{', '.join(dropped)}")
    return chain


def collapsed_rows(triples: np.ndarray, gammas=GAMMAS) -> dict[float, int]:
    """Row index of (1, cos g, cos g) for each requested separation."""
    rows = {}
    for g in gammas:
        target = np.array([1.0, np.cos(g * ARCMIN), np.cos(g * ARCMIN)])
        distance = np.linalg.norm(triples - target, axis=1)
        index = int(np.argmin(distance))
        if distance[index] > 1e-7:
            raise SystemExit(f"gamma = {g}' is not a table row (nearest {distance[index]:.2e})")
        rows[g] = index
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pieces", type=Path, required=True)
    ap.add_argument("--low", type=Path, default=None,
                    help="LOW piece; added to every partial sum if given")
    ap.add_argument("--models", nargs="+", default=["tree", "bihalofit"])
    ap.add_argument("--channel", default="zeta_TTT")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    low = np.load(args.low, allow_pickle=False) if args.low else None
    record = {}
    for model in args.models:
        bands = load_bands(args.pieces, model)
        if not bands:
            print(f"[cut] no bands for {model}, skipping")
            continue
        triples = bands[0]["data"]["cosine_triples"]
        # chi_shells is COMOVING DISTANCE IN Mpc. The bispectrum is evaluated
        # in h units, so the multipole-to-wavenumber map needs chi in Mpc/h.
        # Getting this wrong misplaces every mode by a factor 1/h = 1.49.
        meta = load_meta(bands[0]["data"]["cosmo_meta"])
        h = float(meta["h"])
        chi_mpc = np.asarray(bands[0]["data"]["chi_shells"], float)
        chi_h = chi_mpc * h
        rows = collapsed_rows(triples)

        running = (np.asarray(low[args.channel], float) if low is not None
                   else np.zeros_like(np.asarray(bands[0]["data"][args.channel], float)))
        edges, curves = [], []
        for band in bands:
            running = running + np.asarray(band["data"][args.channel], float)
            edges.append(band["hi"])
            curves.append(running.copy())
        edges = np.asarray(edges)
        curves = np.stack(curves)          # (n_band, n_row, n_shell)

        print(f"\n=== {model}: {args.channel}, partial sums up to each cutoff ===")
        print(f"    shells chi_h [Mpc/h]: {chi_h[0]:.0f} (near) ... {chi_h[-1]:.0f} (source)")
        for g, row in rows.items():
            reference = curves[np.argmin(np.abs(edges - 960)), row, -1]
            line = "  ".join(
                f"{e:6d}:{curves[i, row, -1] / reference:7.2f}"
                for i, e in enumerate(edges))
            print(f"  gamma={g:7.2f}'  source shell, ratio to cutoff 960: {line}")

        record[f"{model}_edges"] = edges
        record[f"{model}_curves"] = curves
        record[f"{model}_chi_h"] = chi_h
        record[f"{model}_chi_mpc"] = chi_mpc
        record[f"{model}_h"] = np.float64(h)
        record[f"{model}_rows"] = np.asarray(list(rows.values()))

        print("  wavenumbers [h/Mpc]. soft = ell / chi_h is the mode a leg at "
              "the cutoff carries;")
        print("  hard = 2 ell / chi_h is the largest wavenumber the P(k) table "
              "must supply.")
        for label, index in (("near shell", 0), ("source shell", -1)):
            soft = edges / chi_h[index]
            hard = 2 * edges / chi_h[index]
            print(f"    {label:12s} chi={chi_mpc[index]:7.1f} Mpc "
                  f"({chi_h[index]:7.1f} Mpc/h)")
            print("      soft " + "  ".join(f"{e}:{k:6.2f}" for e, k in zip(edges, soft)))
            print("      hard " + "  ".join(f"{e}:{k:6.2f}" for e, k in zip(edges, hard)))
        print(f"    (BiHalofit calibrated to k < {BIHALOFIT_CALIBRATION_KMAX} h/Mpc; "
              f"baryons matter above k ~ {BARYON_SCALE_K} h/Mpc)")

    if args.out and record:
        record["gammas"] = np.asarray(GAMMAS)
        np.savez(args.out, **record)
        print(f"\n[cut] -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
