"""Read a folded sweep NPZ and form the four plotted observables.

The sweep stores one row per (component pair, order, separation, source
distance) in the real (Phi_00, Re Psi_0, Im Psi_0) basis. The paper's
figures plot four combinations of those component pairs. The definitions
here are copied from the figure generator
`analyses/analysis3/plot_analysis3_nlo_decomposition.py:76-85`, so a curve
built with this module is the same curve the paper plots.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

ARCMIN = 180.0 * 60.0 / np.pi

#: name -> list of (component pair, coefficient).
OBSERVABLES: dict[str, list[tuple[tuple[int, int], float]]] = {
    "xi_kappa": [((0, 0), +1.0)],
    "xi_plus": [((1, 1), +1.0), ((2, 2), +1.0)],
    "xi_minus": [((1, 1), +1.0), ((2, 2), -1.0)],
    "xi_kappa_gamma_t": [((0, 1), -1.0)],
}


def _gamma_arcmin(x, y) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    x = x / np.linalg.norm(x)
    y = y / np.linalg.norm(y)
    return float(np.arccos(np.clip(x @ y, -1.0, 1.0)) * ARCMIN)


def load_pairs(path: Path, order: int, t_final: float | None = None):
    """(gamma_arcmin, {(a, b): values}) for one order and source distance.

    `allow_pickle` is required because the sweep stores the direction
    vectors as object arrays. These files are written by the local sft-wick
    run; they are never downloaded or user supplied.
    """
    with np.load(path, allow_pickle=True) as d:
        a, b, o = np.asarray(d["a"]), np.asarray(d["b"]), np.asarray(d["order"])
        v = np.asarray(d["value"], float)
        x, y = d["x"], d["y"]
        t = np.asarray(d["t_final"], float)

    keep = o == order
    if t_final is not None:
        keep &= np.isclose(t, t_final, rtol=1e-9)
    if not keep.any():
        raise ValueError(f"{path.name}: no rows at order={order}, t_final={t_final}")

    grouped, gamma = {}, None
    for pair in {(int(ai), int(bi)) for ai, bi in zip(a[keep], b[keep])}:
        m = keep & (a == pair[0]) & (b == pair[1])
        idx = np.flatnonzero(m)
        g = np.array([_gamma_arcmin(x[i], y[i]) for i in idx])
        s = np.argsort(g)
        if gamma is None:
            gamma = g[s]
        grouped[pair] = v[m][s]
    return gamma, grouped


def combine(grouped: dict, name: str) -> np.ndarray | None:
    """One named observable from the component pairs, or None if incomplete."""
    total = None
    for pair, sign in OBSERVABLES[name]:
        if pair not in grouped:
            return None
        term = sign * grouped[pair]
        total = term if total is None else total + term
    return total


def load_observables(path: Path, order: int, t_final: float | None = None):
    """(gamma_arcmin, {observable name: values})."""
    gamma, grouped = load_pairs(path, order, t_final)
    return gamma, {name: combine(grouped, name) for name in OBSERVABLES}


def source_distances(path: Path) -> np.ndarray:
    with np.load(path, allow_pickle=True) as d:
        return np.unique(np.asarray(d["t_final"], float))
