"""LOW-only exact-Wigner-3j kappa3 build (no flat-sky HIGH branch).

Produces a pure exact-3j non-contracted reference table for the clean
windowed-bug isolation (the flat-sky HIGH branch grossly undershoots the
exact-3j at the contact apex / low ell, so it is omitted entirely).

Reuses the build_equal_time_limber_table.py setup helpers; only LOW
(compute_kappa3_zeta_table, ell_max=48) is computed + saved. Metadata mirrors
the canonical build so the equal_time callable consumes it unchanged.
"""
from __future__ import annotations

import sys
import time
import warnings
from dataclasses import replace
from pathlib import Path

import numpy as np

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

import build_equal_time_limber_table as B  # noqa: E402  setup helpers
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.sachs import compute_kappa3_zeta_table  # noqa: E402

ELL_MAX = 48
L2NPZ = Path(
    "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/callables/"
    "kappa3_vertex/R_contracted_lcap48/"
    "windowed_squeezed_kappa3_z5covgrid16_omega0316_h06711_physMpc.npz"
)
OUT = Path(
    "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/callables/"
    "kappa3_vertex/equal_time_exact3j/"
    "equal_time_exact3j_kappa3_z5covgrid16_omega0316_h06711.npz"
)


def main() -> int:
    cosmo = FiducialCosmology()
    lam_phys, lam_h = B._load_l2_lambda_grid_h(L2NPZ, cosmo.h)
    pk_delta = load_pk_delta(n_s=cosmo.n_s)
    triples = B._triangle_filtered_triples(15)
    channels = ("TTT", "TTP", "TPP", "PPP")
    print(f"[exact3j] cosmo Om={cosmo.Omega_m}, h={cosmo.h}; lam n={lam_phys.size} "
          f"[{lam_phys[0]:.1f}..{lam_phys[-1]:.1f}]; triples={triples.shape[0]}; "
          f"ell_max(exact-3j)={ELL_MAX}; measure=lambda")
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        low = compute_kappa3_zeta_table(
            triples, lam_h, pk_matter_today=pk_delta, cosmo=cosmo,
            ell_max=ELL_MAX, Nmax=128, spt_kind="tree_phi",
            units="physical", lambda_convention="project",
            radial_discretization="sample", radial_measure="lambda",
            channels=channels, parallelism="sequential", shell_block_size=4,
            triple_chunk_size=4096, n_workers=1, fuse_tree_phi_channels=True,
            verbose=True)
    meta = dict(low.cosmo_meta)
    meta.update(
        producer="build_exact3j_lowonly.py",
        artifact="equal_time_EXACT3J_lambda_NON_R_contracted",
        equal_time=True, limber_high=False,
        already_R_convolved=False, already_R_contracted=False,
        radial_measure="lambda", radial_density_measure="lambda",
        radial_kernel_convention="lambda_three_leg_native_density",
        sachs_prefactor_per_leg="(1+z)^4 bare driving field (SACHS_GEOMETRIC_POWER=4)",
        poisson_policy="P_matter input; tree_phi applies A(a) per leg internally",
        query_coordinate="lambda_project_mpc",
        l2_lambda_grid_npz=L2NPZ.name,
        ell_cut=ELL_MAX, ell_high_max=ELL_MAX, n_cosine=15,
        high_branch="OMITTED (pure exact-Wigner-3j LOW only)",
    )
    low = replace(low, cosmo_meta=meta)
    for ch in channels:
        arr = np.asarray(getattr(low, f"zeta_{ch}"))
        print(f"[exact3j] zeta_{ch}: shape={arr.shape}, "
              f"nonzero={int(np.count_nonzero(arr))}/{arr.size}, "
              f"|max|={np.max(np.abs(arr)):.3e}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    low.save_npz(OUT)
    print(f"[exact3j] total {time.perf_counter()-t0:.1f}s -> {OUT} "
          f"({OUT.stat().st_size/1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
