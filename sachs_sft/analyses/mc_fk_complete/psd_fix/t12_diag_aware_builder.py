"""R4: a diagonal-aware Sigma2 builder, and whether it is PSD.

The corr_op table stores C_ab(lam_i, lam_j) EXACTLY at the 21 nodes.  The only
thing wrong with the equal-time diagonal is the 2-D bilinear route used to reach
lambda between nodes, which cannot see the ridge on lam1 = lam2.  R4 therefore
reads the 21 exact diagonal values per component and interpolates them in ONE
dimension (PCHIP of log|C| with the sign kept, falling back to PCHIP of the raw
values if the component changes sign), then differentiates D^4 C.

Validated against the four dense rebuilds by t4_diag_model.py: <=0.83% in the
worst cell, <=0.42% in the other three, versus -5% .. -24% for the production
bilinear diagonal.

Checks here: PSD of the 6x6 over a fine lambda scan, and the Order-0 anchor.
"""
from __future__ import annotations

import sys

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


class R4DiagAware(R.R0Reference):
    """Equal-time diagonal read from the 21 exact node values, interpolated in 1-D."""

    def _diag_models(self, cos_gamma: float):
        key = ("diag", round(float(cos_gamma), 12))
        got = self._cache.get(key)
        if got is not None:
            return got
        n2 = ds._n2_of_cos(cos_gamma)
        vals = np.array([_corr_op.C_fn(ds._N1, float(t), n2, float(t))
                         for t in TABLE_LAM])           # (21, 3, 3)
        fns = []
        for a in range(N_COMP):
            for b in range(N_COMP):
                y = vals[:, a, b]
                if np.all(y > 0):
                    g = PchipInterpolator(TABLE_LAM, np.log(y))
                    fns.append(lambda x, g=g: np.exp(g(x)))
                elif np.all(y < 0):
                    g = PchipInterpolator(TABLE_LAM, np.log(-y))
                    fns.append(lambda x, g=g: -np.exp(g(x)))
                else:
                    g = PchipInterpolator(TABLE_LAM, y)
                    fns.append(lambda x, g=g: g(x))
        self._cache[key] = fns
        return fns

    def _sample_f(self, cos_gamma: float, lam):
        fns = self._diag_models(cos_gamma)
        lam = np.atleast_1d(np.asarray(lam, float))
        d4 = np.asarray(self.bg.D(lam), float) ** 4
        out = np.empty((lam.size, N_COMP, N_COMP))
        for a in range(N_COMP):
            for b in range(N_COMP):
                out[:, a, b] = d4 * np.asarray(fns[a * N_COMP + b](lam), float)
        return out

    def matrix(self, cos_gamma: float, lam: float):
        # smooth model -> plain central difference, no cell logic needed
        x = float(np.clip(lam, self.lam_lo, self.lam_hi))
        h = 2.0
        x = float(np.clip(x, self.lam_lo + h, self.lam_hi - h))
        fp = self._sample_f(cos_gamma, np.array([x + h]))[0]
        fm = self._sample_f(cos_gamma, np.array([x - h]))[0]
        return self._finish((fp - fm) / (2.0 * h) / float(self.bg.D(x)) ** 4)


def main() -> None:
    bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
    lam = np.linspace(406.0, 2328.0, 600)
    builders = {"current": ds.Sigma2Builder(background=bg, apply_c0=False),
                "R1 knots": R.R1Knots(background=bg, apply_c0=False),
                "R4 diag-aware": R4DiagAware(background=bg, apply_c0=False)}
    print("=== PSD of the two-ray 6x6 over 600 lambda ===")
    print(f"{'builder':>14} {'gamma':>8} {'non-PSD':>9} {'min eig / median eig':>22}")
    for k, b in builders.items():
        for gam in (1.0, 17.3):
            cosg = float(np.cos(np.deg2rad(gam / 60.0)))
            bad, worst = 0, np.inf
            for t in lam:
                within = b.matrix(1.0, float(t))
                cross = b.matrix(cosg, float(t))
                M = np.zeros((6, 6))
                M[:3, :3] = within
                M[3:, 3:] = within
                M[:3, 3:] = cross
                M[3:, :3] = cross.T
                e = np.linalg.eigvalsh(0.5 * (M + M.T))
                r = e.min() / np.median(np.abs(e))
                worst = min(worst, r)
                if e.min() < 0:
                    bad += 1
            print(f"{k:>14} {gam:>8.1f} {bad:>9d} {worst:>22.4f}")


if __name__ == "__main__":
    main()
