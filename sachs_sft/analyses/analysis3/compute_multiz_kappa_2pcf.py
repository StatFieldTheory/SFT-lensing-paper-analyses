"""Compute the convergence 2PCF xi_kappa(gamma) decomposed into O0 / FF / FK at
the five source redshifts of the multi-z kappa figure (z_s = 1, 1.7, 2.5, 3.2,
4.0), reusing the validated analysis-3 sft-wick configs.

Only the per-z source distance lambda(z_s) (= sweep t_final) changes; the
propagator / kappa3 callable tables are z-independent and reused.  The heavy
lifting (run_vertex_sweep, z->lambda mapping) is imported from ``multiz_sweep.py``
in this folder (the talk repository's ``run_multiz_components.py`` until 2026-09-02).

Run with the sft-wick env:
  <sft-wick python> compute_multiz_kappa_2pcf.py
Output: outputs/multiz_kappa_2pcf_5z.npz  (z, gamma, o0, ff, fk)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import multiz_sweep as R  # noqa: E402  (moved here from talk/scripts/run_multiz_components.py, 2026-09-02)

Z = np.array([1.0, 1.7, 2.5, 3.2, 4.0])
OUT = HERE / "outputs" / "multiz_kappa_2pcf_5z.npz"


def main() -> int:
    lam_of_z = R.build_lambda_of_z()
    lam = np.asarray(lam_of_z(Z), dtype=float)
    print("z grid :", Z.tolist(), flush=True)
    print("lambda :", np.round(lam, 3).tolist(), flush=True)

    comp: dict[str, np.ndarray] = {}
    gamma_ref = None
    for label in ("o0", "ff", "fk"):
        gamma, xi = R.run_vertex_sweep(R.CONFIGS[label], lam, label)
        if gamma_ref is None:
            gamma_ref = gamma
        elif not np.allclose(gamma, gamma_ref, rtol=1e-6, atol=1e-6):
            raise RuntimeError(f"gamma grid mismatch for {label}")
        comp[label] = xi

    assert gamma_ref is not None
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, z=Z, gamma=gamma_ref, lam=lam,
             o0=comp["o0"], ff=comp["ff"], fk=comp["fk"])
    o0, fk = comp["o0"], comp["fk"]
    for i, z in enumerate(Z):
        print(f"  z={z:.2f}: O0(0.5')={o0[i,0]:.3e}  "
              f"|FK|/O0(0.5')={abs(fk[i,0])/o0[i,0]:.4f}", flush=True)
    print(f"saved {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
