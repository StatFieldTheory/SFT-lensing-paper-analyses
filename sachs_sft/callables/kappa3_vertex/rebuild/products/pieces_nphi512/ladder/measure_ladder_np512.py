"""Cutoff ladder FK/Order-0 at 0.5', old (n_phi=64) vs new (n_phi=512).

Measurement and extrapolation are verbatim the deployed recipe in
SFT-lensing-paper-analyses/reproduce/check_numbers.py (lines 108-129):
  value = 100 * xi_kappa_FK[0] / xi_kappa_O0[0]   (grid point 0 == 0.5')
  inc = diff(v); r_last = inc[-1]/inc[-2]
  extrap = v[-1] + inc[-1]*r_last/(1-r_last)
"""
import math, sys
from pathlib import Path
import numpy as np

REPO = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing")
SACHS = REPO / "SFT-lensing-paper-analyses" / "sachs_sft"
RUNS = SACHS / "sftwick_outputs" / "2PCF"
O0P = RUNS / "C_corr_op_O0" / "xi_C_corr_op_O0.npz"
OLD = RUNS / "cutoff_ladder"
NEWD = SACHS / "callables/kappa3_vertex/rebuild/products/pieces_nphi512/ladder"
NEWTOP = SACHS / "callables/kappa3_vertex/rebuild/products/pieces_nphi512/table_permclosed_np512_xi.npz"
ARCMIN = 180 * 60 / np.pi
OBS = {"xi_kappa": [((0, 0), 1.0)], "xi_plus": [((1, 1), 1.0), ((2, 2), 1.0)],
       "xi_minus": [((1, 1), 1.0), ((2, 2), -1.0)], "xi_kappa_gamma_t": [((0, 1), -1.0)]}


def load_pairs(path, order):
    d = np.load(path, allow_pickle=True)
    a, b, o = np.asarray(d["a"]), np.asarray(d["b"]), np.asarray(d["order"])
    v = np.asarray(d["value"], float); x, y = d["x"], d["y"]
    keep = o == order
    grouped, gamma = {}, None
    for pair in {(int(ai), int(bi)) for ai, bi in zip(a[keep], b[keep])}:
        m = keep & (a == pair[0]) & (b == pair[1]); idx = np.flatnonzero(m)
        g = np.array([math.acos(float(np.clip(np.dot(
            np.asarray(x[i], float) / np.linalg.norm(np.asarray(x[i], float)),
            np.asarray(y[i], float) / np.linalg.norm(np.asarray(y[i], float))), -1, 1))) * ARCMIN
            for i in idx])
        s = np.argsort(g)
        if gamma is None:
            gamma = g[s]
        grouped[pair] = v[m][s]
    return gamma, grouped


def obs(path, order):
    g, gr = load_pairs(path, order)
    return g, {n: sum(sg * gr[p] for p, sg in combos) for n, combos in OBS.items()
               if all(p in gr for p, _ in combos)}


CUTS = (960, 1920, 3840, 7680, 15360)
kk = "xi_kappa"
g, o0 = obs(O0P, 0)
assert abs(g[0] - 0.5) < 1e-6, g[0]
BAND = np.array([2.0, 5.0, 8.0, 12.0])
lin = lambda gg, arr, x: float(np.interp(x, gg, arr))

def ladder(paths):
    v, bands = [], []
    for p in paths:
        gg, ob = obs(p, 2)
        v.append(100 * ob[kk][0] / o0[kk][0])
        bands.append([100 * lin(gg, ob[kk], x) / lin(g, o0[kk], x) for x in BAND])
    return np.array(v), np.array(bands)

def extrap(v):
    inc = np.diff(v)
    r_last = inc[-1] / inc[-2]
    return v[-1] + inc[-1] * r_last / (1 - r_last), r_last, inc


old_paths = [OLD / f"table_tree_cut{c}_r4_xi.npz" for c in CUTS]
new_paths = [NEWD / f"table_tree_np512_cut{c}_xi.npz" for c in CUTS[:-1]] + [NEWTOP]

for label, paths in (("DEPLOYED n_phi=64 (r4 rows)", old_paths),
                     ("REBUILT  n_phi=512 (permclosed)", new_paths)):
    if not all(Path(p).exists() for p in paths):
        missing = [p.name for p in paths if not Path(p).exists()]
        print(f"\n### {label}: MISSING {missing}")
        continue
    v, bv = ladder(paths)
    e, r_last, inc = extrap(v)
    print(f"\n### {label}")
    print("  ladder : " + "  ".join(f"{c}:{x:.3f}%" for c, x in zip(CUTS, v)))
    print("  increments        : " + "  ".join(f"{x:+.4f}" for x in inc))
    print("  successive v-ratios: " + "  ".join(f"{x:.3f}" for x in (v[1:] / v[:-1])))
    print("  increment ratios   : " + "  ".join(f"{x:.3f}" for x in (inc[1:] / inc[:-1])))
    print(f"  r_last = inc[-1]/inc[-2] = {r_last:.4f}")
    print(f"  geometric extrap ell_max->inf = {e:.3f}%  "
          f"(= {100*(e/v[-1]-1):.1f}% above the {v[-1]:.3f}% cutoff-15360 value)")
    print(f"  final doubling adds {100*(v[-1]/v[-2]-1):.2f}% at 0.5'")
    print("  final doubling across 2'-12': " +
          ", ".join(f"{x:.0f}':{100*(bv[-1][i]/bv[-2][i]-1):.2f}%" for i, x in enumerate(BAND)))
    print("  band FK/O0 at ell_max=15360: " +
          ", ".join(f"{x:.0f}':{bv[-1][i]:.3f}%" for i, x in enumerate(BAND)))
