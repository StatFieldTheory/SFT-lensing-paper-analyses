"""Refold the multi-z FK planes with the rebuilt n_phi=512 kappa3 table.

Mirrors analyses/revision_2026-08/regenerate_multiz.py exactly (same
materialise_variant, same multiz_sweep.run_vertex_sweep, same perm-aware
callable) but

  * adds the z_s = 5 baseline plane lambda = 2313.0288751857356 so the six
    planes come out of ONE sweep, and
  * writes to a private npz so the deployed
    analysis3/outputs/multiz_kappa_2pcf_5z.npz is not touched.

Only FK is computed; Order-0 and FF carry no three-point vertex.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

SACHS = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft")
A3 = SACHS / "analyses" / "analysis3"
REBUILD = SACHS / "callables" / "kappa3_vertex" / "rebuild"
PERMAWARE = SACHS / "callables" / "kappa3_vertex" / "equal_time_limber_cut15360_permaware"

TABLE = REBUILD / "products" / "pieces_nphi512" / "table_permclosed_np512.npz"
CALLABLE = PERMAWARE / "perm_aware_kappa3_callable.py"
WORK = REBUILD / "products" / "pieces_nphi512" / "_variant_multiz_np512"
OUT = REBUILD / "products" / "pieces_nphi512" / "multiz_fk_np512_6planes.npz"

LAM_BASELINE_Z5 = 2313.0288751857356
Z5 = np.array([1.0, 1.7, 2.5, 3.2, 4.0])

sys.path.insert(0, str(A3))
sys.path.insert(0, str(REBUILD))
import multiz_sweep as R          # noqa: E402
import run_fk_variant as rfv      # noqa: E402


def main() -> int:
    lam_of_z = R.build_lambda_of_z()
    lam5 = np.asarray(lam_of_z(Z5), dtype=float)
    lam = np.concatenate([lam5, [LAM_BASELINE_Z5]])
    print("[refold] z      :", Z5.tolist(), "+ z_s=5 baseline plane", flush=True)
    print("[refold] lambda :", np.round(lam, 4).tolist(), flush=True)

    _, config = rfv.materialise_variant(
        TABLE.resolve(), WORK.resolve(), n_jobs=6,
        callable_src=CALLABLE.resolve())
    print(f"[refold] config: {config}", flush=True)

    t0 = time.perf_counter()
    gamma, xi = R.run_vertex_sweep(config, lam, "fk")
    print(f"[refold] sweep done in {time.perf_counter() - t0:.1f} s", flush=True)

    np.savez(OUT, z=Z5, lam=lam, gamma=gamma, fk=xi,
             lam_baseline_z5=LAM_BASELINE_Z5,
             table=str(TABLE), callable_src=str(CALLABLE))
    print(f"[refold] saved {OUT}  fk{xi.shape}", flush=True)
    for i, l in enumerate(lam):
        print(f"   lam={l:.4f}  FK(0.5')={xi[i,0]:.6e}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
