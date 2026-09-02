"""Build the validation figure for the revision: FF markers only, corrected FK line.

Option A of the fig-6 decision: the FK channel is validated by the deterministic
kernel cross-check (code/rebuild/fk_kernel_crosscheck.py), not by Monte-Carlo
markers, so the figure shows the three analytic channels with the FF Monte-Carlo
overlay at the published 15-point gamma grid.  The FF values come from the
existing cache (seed-deterministic, identical settings to the published figure);
the FK line is rebound to the corrected converged fold, exactly as for the other
regenerated figures.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[2]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fk-xi", type=Path, required=True,
                    help="the corrected folded FK sweep")
    ap.add_argument("--fk-markers", type=Path, default=None,
                    help="pooled FK Monte-Carlo markers from "
                         "code/mc_fk_complete/_markers_pooled.npz; omit for "
                         "the FF-markers-only version")
    args = ap.parse_args()

    sys.path.insert(0, str(_MC))
    spec = importlib.util.spec_from_file_location(
        "fig_xi_channels", _MC / "fig_xi_channels.py")
    fig = importlib.util.module_from_spec(spec)
    sys.modules["fig_xi_channels"] = fig
    spec.loader.exec_module(fig)

    original = fig._load_kk

    def patched(name: str):
        if name != "C_corr_op_K_limber_FK":
            return original(name)
        d = np.load(args.fk_xi, allow_pickle=True)
        m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
        gammas = []
        for xi, yi in zip(d["x"][m], d["y"][m]):
            xi = np.asarray(xi, float); yi = np.asarray(yi, float)
            gammas.append(math.degrees(math.acos(float(np.clip(
                np.dot(xi / np.linalg.norm(xi), yi / np.linalg.norm(yi)),
                -1, 1)))) * 60)
        order = np.argsort(gammas)
        return np.asarray(gammas)[order], np.asarray(d["value"], float)[m][order]

    fig._load_kk = patched
    print(f"[val-fig] FK line -> {args.fk_xi.name}")

    g, o0 = fig._load_kk("C_corr_op_O0")
    _, ff_mom = fig._load_kk("C_corr_op_K_limber_FF")
    _, fk3 = fig._load_kk("C_corr_op_K_limber_FK")

    cache = dict(np.load(fig._cache_path(), allow_pickle=True))
    g_c = np.asarray(cache["g_mc"], float)
    ff_c = np.asarray(cache["ffc_mc"], float)
    se_c = np.asarray(cache["ffc_se"], float)
    g_mc, ffc_mc, ffc_se = [], [], []
    for gv in fig.GAMMA_MC:
        j = int(np.argmin(np.abs(g_c - gv)))
        if abs(g_c[j] - gv) > 1e-6:
            raise SystemExit(f"FF value for gamma={gv}' not in the cache")
        g_mc.append(gv); ffc_mc.append(ff_c[j]); ffc_se.append(se_c[j])
    print(f"[val-fig] FF markers at {len(g_mc)} gammas "
          f"({g_mc[0]}'..{g_mc[-1]}'), from the cached seed-11 runs")

    D = dict(g=g, o0=o0, ff_mom=ff_mom, fk3=fk3, g_mc=np.array(g_mc),
             ffc_mc=np.array(ffc_mc), ffc_se=np.array(ffc_se))
    if args.fk_markers is not None:
        m = np.load(args.fk_markers)
        D.update(g_fk=np.asarray(m["gamma"], float),
                 fk_mc=np.asarray(m["fk"], float),
                 fk_se=np.asarray(m["err"], float))
        print(f"[val-fig] FK markers at {len(D['g_fk'])} gammas "
              f"({D['g_fk'][0]}'..{D['g_fk'][-1]}'), pooled over "
              f"{int(m['n_blocks'])} independent seed blocks x "
              f"{int(m['seeds_per_block'])} seeds; errors are the block scatter")
    np.savez(fig._cache_path(), **D)
    print(f"[val-fig] cache rewritten -> {fig._cache_path()}")
    fig._plot(D)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
