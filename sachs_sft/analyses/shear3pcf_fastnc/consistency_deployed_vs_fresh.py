"""Consistency cross-check: the targeted FRESH query reproduces the DEPLOYED
equal-time kappa3 callable table at the cosine nodes they SHARE.

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

WHY THIS EXISTS
---------------
STAGE 0/1 do NOT re-generate the deployed callable; they call the SAME canoes
functions that BUILT it (compute_kappa3_zeta_table / compute_kappa3_sigma3_high /
compute_kappa3_mod_*_equal_time) at the specific COLLAPSED triples (cos t3, cos g,
cos g) the validation needs.  Those triples live in the cos->1 squeezed corner
that the deployed table's 15-node cosine grid (linspace(-1,1,15); only node with
cos>0.99 is cos=1.0) cannot resolve -- hence a targeted fresh query rather than a
read+interpolate.

This script proves the targeted query is FAITHFUL to the deployed callable: it
re-runs the fresh query at a handful of deployed grid nodes (exact same triples,
same lambda shells, same build parameters) and checks zeta_TTT and zeta_Dmod
agree with the stored table.  LOW (exact Wigner-3j) should match to ~machine
precision; HIGH (Limber quadrature) to <~1e-2.  This closes the "why not just read
the npz" question with a number rather than an assertion.

Deployed-build parameters (from build_equal_time_limber_table.py argparse defaults
+ cosmo_meta): Nmax=128, ell_cut=60, ell_high_max=1000, ell_min=1e-3, n_ell=96,
n_phi=64, radial_measure='lambda', units='physical', lambda_convention='project'.
The lambda grid is the deployed lambda_shells_Mpc (physical Mpc) * h.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

_HERE = Path(__file__).resolve().parent
_ETL = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
            "callables/kappa3_vertex/equal_time_limber")
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

_DEPLOYED_NPZ = _ETL / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz"

# How many deployed nodes to re-query (kept small: each is exact + a Limber pass).
_N_NODES = 3
# Targets we WANT represented among the chosen deployed nodes (nearest match used).
_TARGET_TRIPLES = np.array([
    [1.000000, 0.857143, 0.857143],   # collapsed node nearest our cos->1 regime
    [1.000000, 1.000000, 1.000000],   # fully degenerate (gamma->0 collapsed limit)
    [0.857143, 0.857143, 0.714286],   # a generic non-degenerate small triangle
])


def _pick_deployed_nodes(cosine_triples: np.ndarray) -> np.ndarray:
    """Return row indices of the deployed triples nearest each target (unique)."""
    idx = []
    for tgt in _TARGET_TRIPLES:
        j = int(np.argmin(np.sum((cosine_triples - tgt[None, :]) ** 2, axis=1)))
        if j not in idx:
            idx.append(j)
        if len(idx) >= _N_NODES:
            break
    return np.array(idx, dtype=int)


def main() -> int:
    from canoes.sachs import (  # noqa: E402
        compute_kappa3_sigma3_high,
        compute_kappa3_zeta_table,
        kappa3_combine_low_high,
    )
    from canoes.sachs.kappa3 import (  # noqa: E402
        compute_kappa3_mod_zeta_equal_time,
        compute_kappa3_mod_sigma3_high,
    )

    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)

    # allow_pickle=False: this script reads only plain numeric arrays
    # (cosine_triples / lambda_shells_Mpc / zeta_TTT / zeta_Dmod). The deployed
    # npz is our own trusted build artifact, and we never touch its object-dtype
    # cosmo_meta here, so no unpickling is needed.
    dep = np.load(_DEPLOYED_NPZ, allow_pickle=False)
    ct = np.asarray(dep["cosine_triples"], dtype=np.float64)         # (1671, 3)
    lam_phys = np.asarray(dep["lambda_shells_Mpc"], dtype=np.float64)  # (16,) physical Mpc
    lam_h = lam_phys * cosmo.h
    zTTT_dep = np.asarray(dep["zeta_TTT"], dtype=np.float64)          # (1671, 16)
    zDmod_dep = np.asarray(dep["zeta_Dmod"], dtype=np.float64)        # (1671, 16)

    rows = _pick_deployed_nodes(ct)
    triples = ct[rows]                                               # (k, 3) EXACT deployed nodes
    print(f"[consistency] testing {len(rows)} deployed nodes (rows={rows.tolist()}):")
    for r in rows:
        print(f"    cos_triple = {np.array2string(ct[r], precision=6)}")
    print(f"[consistency] lambda shells (physical Mpc): "
          f"[{lam_phys[0]:.1f} .. {lam_phys[-1]:.1f}], n={lam_phys.size}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # --- TTT (LOW exact Wigner-3j + HIGH Limber), EXACT deployed params -----
        low = compute_kappa3_zeta_table(
            triples, lam_h, pk_matter_today=pk, cosmo=cosmo,
            ell_max=60, Nmax=128, spt_kind="tree_phi", units="physical",
            lambda_convention="project", radial_discretization="sample",
            radial_measure="lambda", channels=("TTT", "TTP", "TPP", "PPP"),
            parallelism="sequential", shell_block_size=4, triple_chunk_size=4096,
            n_workers=1, fuse_tree_phi_channels=True)
        high = compute_kappa3_sigma3_high(
            triples, lam_h, pk_matter_today=pk, cosmo=cosmo,
            ell_cut=60, ell_high_max=1000, ell_min=1e-3, n_ell=96, n_phi=64,
            channels=("TTT", "TTP", "TPP", "PPP"), units="physical",
            lambda_convention="project", radial_discretization="sample",
            radial_measure="lambda")
        combined = kappa3_combine_low_high(low, high)
        zTTT_fresh = np.asarray(combined.zeta_TTT, dtype=np.float64)

        # --- Dmod (dominant shear modulus channel) ------------------------------
        bmod_low, dmod_low = compute_kappa3_mod_zeta_equal_time(
            triples, lam_h, pk_matter_today=pk, cosmo=cosmo,
            ell_max=60, Nmax=128, spt_kind="tree_phi", units="physical",
            lambda_convention="project", radial_measure="lambda")
        bmod_high, dmod_high = compute_kappa3_mod_sigma3_high(
            triples, lam_h, pk_matter_today=pk, cosmo=cosmo,
            ell_cut=60, ell_high_max=1000, ell_min=1e-3, n_ell=96, n_phi=64,
            units="physical", lambda_convention="project", radial_measure="lambda")
        zDmod_fresh = np.asarray(dmod_low + dmod_high, dtype=np.float64)

    def _reldiff(fresh, stored):
        denom = np.maximum(np.abs(stored), 1e-300)
        return np.abs(fresh - stored) / denom

    rTTT = _reldiff(zTTT_fresh, zTTT_dep[rows])
    rDmod = _reldiff(zDmod_fresh, zDmod_dep[rows])

    print("\n[consistency] === zeta_TTT  fresh vs deployed (per node, max over shells) ===")
    for i, r in enumerate(rows):
        print(f"    node {r}: max rel diff = {np.max(rTTT[i]):.3e}  "
              f"(deployed[0]={zTTT_dep[r, 0]:+.4e}, fresh[0]={zTTT_fresh[i, 0]:+.4e})")
    print("[consistency] === zeta_Dmod fresh vs deployed (per node, max over shells) ===")
    for i, r in enumerate(rows):
        print(f"    node {r}: max rel diff = {np.max(rDmod[i]):.3e}  "
              f"(deployed[0]={zDmod_dep[r, 0]:+.4e}, fresh[0]={zDmod_fresh[i, 0]:+.4e})")

    out = _HERE / "outputs" / "consistency_deployed_vs_fresh.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out, rows=rows, triples=triples, lam_phys=lam_phys,
             zTTT_fresh=zTTT_fresh, zTTT_dep=zTTT_dep[rows],
             zDmod_fresh=zDmod_fresh, zDmod_dep=zDmod_dep[rows],
             reldiff_TTT=rTTT, reldiff_Dmod=rDmod)

    maxTTT, maxDmod = float(np.max(rTTT)), float(np.max(rDmod))
    print(f"\n[consistency] OVERALL max rel diff: TTT={maxTTT:.3e}, Dmod={maxDmod:.3e}")
    ok = (maxTTT < 1e-2) and (maxDmod < 1e-2)
    print(f"[consistency] VERDICT: {'PASS' if ok else 'FAIL'} "
          f"(fresh query is {'faithful to' if ok else 'INCONSISTENT with'} the deployed callable)")
    print(f"[consistency] saved -> {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
