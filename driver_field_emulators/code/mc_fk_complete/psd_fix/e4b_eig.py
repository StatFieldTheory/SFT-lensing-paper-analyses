"""Conditioning of the assembled 6x6 along the ray: how negative, and where."""
from __future__ import annotations
import sys
import numpy as np
HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, HERE); sys.path.insert(0, HERE + "/psd_fix")
import fk_expect_exact as ex
import sigma2_repaired as R
_CORE = ex._CORE

grid = ex.build_grid(n_lambda=2)
bg = grid.bg
lam = np.linspace(407.0, 2327.0, 4000)          # 0.48 Mpc scan, builder-independent
for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    print(f"=== gamma = {gam}'   fine scan, {lam.size} lambda, 0.48 Mpc step ===")
    for tag, cls in (("baseline", None), ("R1Knots", R.R1Knots)):
        b = (ex._DS.Sigma2Builder(background=bg, apply_c0=False) if cls is None
             else cls(background=bg, apply_c0=False))
        M = np.array([_CORE._assemble_6x6(b, cosg, float(l)) for l in lam])
        w = np.linalg.eigvalsh(M)
        mn, mx = w[:, 0], w[:, -1]
        neg = mn < 0
        rel = mn / mx
        print(f"  {tag:<9} nodes with min-eig<0 : {int(neg.sum()):>4d}/{lam.size} "
              f"({100*neg.mean():.2f}%)   worst min/max eig = {rel.min():+.3e}")
        if neg.any():
            idx = np.where(neg)[0]
            brk = np.where(np.diff(idx) > 1)[0]
            groups = np.split(idx, brk + 1)
            print(f"  {'':<9} windows: " + ", ".join(
                f"[{lam[g[0]]:.1f},{lam[g[-1]]:.1f}] ({lam[g[-1]]-lam[g[0]]+0.48:.1f} Mpc)"
                for g in groups))
