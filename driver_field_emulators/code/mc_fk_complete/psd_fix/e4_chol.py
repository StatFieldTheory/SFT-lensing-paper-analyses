"""How often the eigen-floor branch of ``sachs_mc_core._cholesky_psd`` fires.

The original module is untouched; a counting wrapper is bound onto the module
object at runtime, exactly where ``fk_complete_core`` looks it up
(``_CORE._cholesky_psd``), and the production call pattern
``[_cholesky_psd(V[k]) for k in range(N)]`` is replayed.
"""
from __future__ import annotations
import sys
import numpy as np

HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, HERE); sys.path.insert(0, HERE + "/psd_fix")
import fk_expect_exact as ex
import sigma2_repaired as R
_CORE = ex._CORE
_ORIG = _CORE._cholesky_psd


class Counter:
    def __init__(self): self.calls = self.fb = 0; self.eig = []; self.lam = []
    def bind(self, lam): self._lam = lam
    def __call__(self, M):
        self.calls += 1
        try:
            return np.linalg.cholesky(M)
        except np.linalg.LinAlgError:
            self.fb += 1
            w, V = np.linalg.eigh(M)
            self.eig.append((float(w.min()), float(np.median(np.abs(w)))))
            self.lam.append(float(self._lam[self.calls - 1]))
            return V * np.sqrt(np.maximum(w, 0.0))[None, :]


SIG = 8.0
print(f"{'N':>6} {'dlam':>7} {'builder':>9} {'chol fallbacks':>15} "
      f"{'worst min-eig':>14} {'rel to median|eig|':>19} {'lam of fallbacks'}")
for N in (200, 500, 1000, 1500, 2000):
    for gam in (1.0,):
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        for tag, cls in (("base", None), ("R1", R.R1Knots)):
            grid = ex.build_grid(n_lambda=N)
            if cls is not None:
                grid.builder = cls(background=grid.bg, apply_c0=False)
            V = np.array([_CORE._assemble_6x6(grid.builder, cosg, float(l))
                          for l in grid.lam]) / (2.0 * SIG)
            c = Counter(); c.bind(grid.lam)
            _CORE._cholesky_psd = c
            try:
                _ = [_CORE._cholesky_psd(V[k]) for k in range(N)]
            finally:
                _CORE._cholesky_psd = _ORIG
            if c.eig:
                w = min(c.eig, key=lambda t: t[0])
                s = f"{w[0]:>14.3e} {w[0]/w[1]:>19.2e}"
            else:
                s = f"{'-':>14} {'-':>19}"
            print(f"{N:>6} {grid.dlam:>7.3f} {tag:>9} {c.fb:>7d}/{c.calls:<7d} "
                  f"{s} {np.round(c.lam, 2).tolist()}")
