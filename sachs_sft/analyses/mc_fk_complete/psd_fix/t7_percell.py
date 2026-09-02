"""Which lambda cells carry the Order-0 shift, and how validated is each one.

Applies the corrected (PCHIP-in-log through the 21 exact node values) diagonal
in ONE cell at a time and reports the resulting change in order0_mc.  The two
cells with a dense rebuild behind them (1300-1447, 1700-1900) are marked; the
sum over all cells is compared with the all-cells-at-once number so the
linearity of the decomposition is visible.
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402

ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import sigma2_repaired as R  # noqa: E402
from t5_order0_impact import CorrectedDiag, TABLE_LAM  # noqa: E402

VALIDATED = {(1300.0, 1447.0), (1700.0, 1900.0), (550.0, 700.0), (2160.0, 2250.0)}


def main() -> None:
    bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
    base = R.R0Reference(background=bg, apply_c0=False)
    allc = CorrectedDiag(background=bg, apply_c0=False)
    for gam in (1.0, 17.3):
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        o_base = float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE,
                                    builder=base, n_gauss=256))
        o_all = float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE,
                                   builder=allc, n_gauss=256))
        print(f"\n=== gamma = {gam}'  O0(R0, production C) = {o_base:.6e} ===")
        print(f"{'cell':>16} {'h':>6} {'dO0/O0 [%]':>11}  validated")
        tot = 0.0
        for lo, hi in zip(TABLE_LAM[:-1], TABLE_LAM[1:]):
            b = CorrectedDiag(background=bg, apply_c0=False,
                              cells=[(float(lo), float(hi))])
            o = float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE,
                                   builder=b, n_gauss=256))
            d = 100.0 * (o / o_base - 1.0)
            tot += d
            tag = "  DENSE-REBUILT" if (float(lo), float(hi)) in VALIDATED else ""
            print(f"{lo:7.1f}-{hi:7.1f} {hi-lo:6.1f} {d:11.3f}{tag}")
        print(f"{'sum of cells':>23} {tot:11.3f}")
        print(f"{'all cells at once':>23} {100*(o_all/o_base-1):11.3f}")


if __name__ == "__main__":
    main()
