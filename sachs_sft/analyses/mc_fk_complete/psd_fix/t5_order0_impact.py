"""Order-0 impact of the 21-node lambda grid, compared with the +1.61% ringing fix.

Three Sigma2 sources, all differentiated the SAME way (dense central differences,
cell-aware -- the R0Reference recipe), so the only difference is the C diagonal:

  PROD  : C(lam,lam) as the production callable returns it (bilinear -> Bezier)
  CORR  : C(lam,lam) reconstructed from the 21 EXACT node values with the 1-D
          model validated against the dense rebuilds by t4_diag_model.py

and, for context, the stock ringing spline (`current`) and `R1Knots`.

order0_mc(gamma) = int Sigma2_00(cos gamma, la) G(la)^2 dla  is the anchor
integral, so the reported median ratio to the analysis-3 O0 curve is exactly the
quantity ANCHOR_RATIO = 0.881 calibrates.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402

ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import sigma2_repaired as R  # noqa: E402
import corr_op_C_callable as _corr_op  # noqa: E402

TABLE_LAM = np.load(_corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)
N_COMP = ds.N_COMP


class CorrectedDiag(R.R0Reference):
    """R0 differentiation, but the C diagonal comes from the 1-D node model.

    Only the (0,0) entry is modelled (order0_mc uses Sigma2[0,0]); the other
    entries keep the production values so the object stays a full (3,3).
    """

    def __init__(self, *a, model: str = "pchiplog", cells=None, **kw):
        super().__init__(*a, **kw)
        self.model = model
        self.cells = cells          # None = correct everywhere
        self._diagcache: dict[float, object] = {}

    def _model_fn(self, cos_gamma: float):
        key = round(float(cos_gamma), 12)
        f = self._diagcache.get(key)
        if f is not None:
            return f
        n2 = ds._n2_of_cos(cos_gamma)
        node = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))[0, 0]
                         for t in TABLE_LAM])
        if self.model == "pchiplog" and np.all(node > 0):
            g = PchipInterpolator(TABLE_LAM, np.log(node))
            f = lambda x: np.exp(g(x))                       # noqa: E731
        else:
            g = PchipInterpolator(TABLE_LAM, node)
            f = lambda x: g(x)                               # noqa: E731
        self._diagcache[key] = f
        return f

    def _sample_f(self, cos_gamma: float, lam):
        out = super()._sample_f(cos_gamma, lam)              # D^4 * C, production
        f = self._model_fn(cos_gamma)
        lam = np.asarray(lam, float)
        d4 = np.asarray(self.bg.D(lam), float) ** 4
        new = d4 * np.asarray(f(lam), float)
        if self.cells is None:
            out[:, 0, 0] = new
        else:
            m = np.zeros(lam.shape, bool)
            for lo, hi in self.cells:
                m |= (lam >= lo) & (lam <= hi)
            out[m, 0, 0] = new[m]
        return out


def main() -> None:
    bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
    cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
    d = np.load(cache, allow_pickle=True)
    g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
    GAM = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3, 293.9, 596.9])

    builders = {
        "current": ds.Sigma2Builder(background=bg, apply_c0=False),
        "R1 knots": R.R1Knots(background=bg, apply_c0=False),
        "R0 (prod C)": R.R0Reference(background=bg, apply_c0=False),
        "CORR diag": CorrectedDiag(background=bg, apply_c0=False),
        "CORR cellA+B": CorrectedDiag(background=bg, apply_c0=False,
                                      cells=[(1300.0, 1447.0), (1700.0, 1900.0)]),
    }
    rows = {k: [] for k in builders}
    print(f"{'gamma':>8} {'O0 a3':>12} " + " ".join(f"{k:>13}" for k in builders))
    for gam in GAM:
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        a3 = float(np.interp(gam, g_a3, o0_a3))
        vals = []
        for k, b in builders.items():
            v = float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE,
                                   builder=b, n_gauss=256))
            rows[k].append(v / a3)
            vals.append(v / a3)
        print(f"{gam:>8.2f} {a3:>12.4e} " + " ".join(f"{v:>13.5f}" for v in vals))
    med = {k: float(np.median(rows[k])) for k in rows}
    print("\n" + f"{'median':>21} " + " ".join(f"{med[k]:>13.5f}" for k in builders))
    base = med["current"]
    print(f"{'shift vs current':>21} " +
          " ".join(f"{med[k]/base-1:>+12.3%} " for k in builders))
    r0 = med["R0 (prod C)"]
    print(f"{'shift vs R0':>21} " +
          " ".join(f"{med[k]/r0-1:>+12.3%} " for k in builders))


if __name__ == "__main__":
    main()
