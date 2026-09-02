"""V8: the EXACT proposed patch, transcribed as a subclass, tested.
Nothing outside psd_fix/ is modified."""
import sys, math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline
ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/psd_fix")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import sigma2_repaired as R
import corr_op_C_callable as _corr_op
N_COMP = ds.N_COMP
_N1 = ds._N1
_n2_of_cos = ds._n2_of_cos

_TABLE_LAM = np.load(_corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)

class Patched(ds.Sigma2Builder):
    """Body of the proposed patched Sigma2Builder, verbatim."""
    def __init__(self, background=None, n_nodes: int = 160,
                 lam_lo: float = 406.0, lam_hi: float = 2328.0,
                 apply_c0: bool = True) -> None:
        self.bg = background if background is not None else bgm.Background()
        self.lam_lo = float(lam_lo); self.lam_hi = float(lam_hi)
        self.apply_c0 = bool(apply_c0)
        inner = _TABLE_LAM[(_TABLE_LAM > self.lam_lo + 1e-9)
                           & (_TABLE_LAM < self.lam_hi - 1e-9)]
        self._breaks = np.concatenate([[self.lam_lo], inner, [self.lam_hi]])
        dlam = (self.lam_hi - self.lam_lo) / max(int(n_nodes) - 1, 1)
        self._cell_nodes = [
            np.linspace(lo, hi, max(6, int(np.ceil((hi - lo) / dlam)) + 1))
            for lo, hi in zip(self._breaks[:-1], self._breaks[1:])
        ]
        self._cache = {}

    def _density_splines(self, cos_gamma: float):
        key = round(float(cos_gamma), 12)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        n2 = _n2_of_cos(cos_gamma)
        per_cell = []
        for xs in self._cell_nodes:
            c_diag = np.empty((xs.size, N_COMP, N_COMP), dtype=float)
            for i, t in enumerate(xs):
                c_diag[i] = _corr_op.C_fn(_N1, float(t), n2, float(t))
            f = (np.asarray(self.bg.D(xs), dtype=float) ** 4)[:, None, None] * c_diag
            per_cell.append([CubicSpline(xs, f[:, a, b])
                             for a in range(N_COMP) for b in range(N_COMP)])
        self._cache[key] = per_cell
        return per_cell

    def matrix(self, cos_gamma: float, lam: float):
        per_cell = self._density_splines(cos_gamma)
        lam_c = float(np.clip(lam, self.lam_lo, self.lam_hi))
        j = int(np.clip(np.searchsorted(self._breaks, lam_c) - 1, 0, len(per_cell) - 1))
        flat = per_cell[j]
        d4 = float(self.bg.D(lam_c)) ** 4
        out = np.empty((N_COMP, N_COMP), dtype=float)
        for a in range(N_COMP):
            for b in range(N_COMP):
                out[a, b] = flat[a * N_COMP + b](lam_c, 1) / d4
        out = 0.5 * (out + out.T)
        if self.apply_c0:
            out *= ds.ANCHOR_C0
        return out

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
P  = Patched(background=bg, apply_c0=False)
R1 = R.R1Knots(background=bg, apply_c0=False)
S  = ds.Sigma2Builder(background=bg, apply_c0=False)

# 1. agreement with R1 and node count
print("cells:", len(P._cell_nodes), " total sample points:",
      sum(x.size for x in P._cell_nodes), " (stock uses 160)")
lam = np.linspace(406.0, 2328.0, 601)
for g in (1.0, 17.3):
    cg = math.cos(math.radians(g/60.0))
    dp = np.array([P.matrix(cg, float(l))[0,0] for l in lam])
    dr = np.array([R1.matrix(cg, float(l))[0,0] for l in lam])
    rel = np.abs(dp/dr - 1)
    print(f"  g={g:5.2f}'  patched vs R1Knots: median {np.median(rel):.2e}  max {rel.max():.2e}")

# 2. PSD over a fine scan
def six(b, cg, l):
    w = b.matrix(1.0, l); c = b.matrix(cg, l)
    M = np.zeros((6,6)); M[:3,:3]=w; M[3:,3:]=w; M[:3,3:]=c; M[3:,:3]=c.T
    return M
scan = np.linspace(406.5, 2327.5, 1500)
for g in (0.5, 1.0, 17.3, 114.3):
    cg = math.cos(math.radians(g/60.0))
    bad_s = bad_p = 0; worst_s = worst_p = np.inf
    for l in scan:
        e = np.linalg.eigvalsh(six(S,cg,float(l))); r = e.min()/e.max()
        if e.min() <= 0: bad_s += 1
        worst_s = min(worst_s, r)
        e = np.linalg.eigvalsh(six(P,cg,float(l))); r = e.min()/e.max()
        if e.min() <= 0: bad_p += 1
        worst_p = min(worst_p, r)
    print(f"  g={g:6.2f}'  stock nonPSD {bad_s:4d}/{len(scan)} worst {worst_s:+.3e} | "
          f"patched {bad_p:4d}/{len(scan)} worst {worst_p:+.3e}")
