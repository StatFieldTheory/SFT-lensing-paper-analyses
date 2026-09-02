"""Real per-redshift FF sweep for the multi-z figure, serial (no loky).

The FF slices of the multi-z table have never been computed: the June build
used a talk-era placeholder because the parallel FF order-2 sweep deadlocked
under loky (reproduced 2026-08-26: 1h47m hang with a dead-worker warning,
while the FK sweep with identical settings completed). This runner performs
the genuine computation with sweep.n_jobs=1, which removes loky entirely, one
source distance at a time so progress is visible and a hang loses at most one
slice.

Output: sftwick_outputs/2PCF/multiz/multiz_ff_real.npz  (lam, gamma, ff)
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(_REPO / "talk" / "scripts"))
import run_multiz_components as R  # noqa: E402

LAM = [1822.721, 2094.894, 2216.793, 2266.582, 2297.289]
OUT = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "sftwick_outputs" / "2PCF" / "multiz" / "multiz_ff_real.npz"


def main() -> int:
    rows, gamma_ref = [], None
    for ls in LAM:
        t0 = time.time()
        print(f"[ff-real] lam_source = {ls} ...", flush=True)
        gamma, xi = R.run_vertex_sweep(R.CONFIGS["ff"], np.array([ls]), "ff",
                                       n_jobs=1)
        if gamma_ref is None:
            gamma_ref = gamma
        rows.append(np.asarray(xi[0], float))
        print(f"[ff-real] lam_source = {ls} done in {time.time()-t0:.0f} s; "
              f"xi(0.5')={rows[-1][0]:+.4e}  xi(5000')={rows[-1][-1]:+.4e}",
              flush=True)
        np.savez(OUT, lam=np.array(LAM[:len(rows)]), gamma=gamma_ref,
                 ff=np.array(rows))
    print(f"[ff-real] all saved -> {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
