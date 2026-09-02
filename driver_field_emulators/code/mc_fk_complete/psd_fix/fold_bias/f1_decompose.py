"""Decompose the Order-0 fold over its 24x24 Gauss-Legendre grid.

The fold measured empirically (bit-exact against the production npz) is

    xi_ab(gamma) = int_0^{lam_f} dlam1 int_0^{lam_f} dlam2  C_ab(n1,lam1; n2,lam2)

with a tensor-product 24-node Gauss-Legendre rule on [0, lam_f]^2.
The corr_op table spans lam in [406, 2328]; RegularGridInterpolator is called
with bounds_error=False, fill_value=None, i.e. LINEAR EXTRAPOLATION outside.
"""
from __future__ import annotations
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, "/Users/zzhang/projects/angular_statistics/canoes/src")
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
                   "sachs_sft/callables/C_propagator/corr_op")
import corr_op_C_callable as co  # noqa: E402

LF = 2313.0288751857356
G21 = np.load(co.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)
N1 = np.array([0.0, 0.0, 1.0])

GAMMAS = {  # (arcmin, y-direction)
    "0.50'":   np.array([0.0001454441, 0.0, 0.9999999894]),
    "17.28'":  np.array([0.0050252287, 0.0, 0.9999873735]),
    "293.90'": np.array([0.0853881735, 0.0, 0.9963477605]),
}


def gl_nodes(n: int = 24):
    x, w = leggauss(n)
    return (x + 1) / 2 * LF, w / 2 * LF


def main() -> None:
    lam, jw = gl_nodes(24)
    lo, hi = G21[0], G21[-1]
    below = lam < lo
    print(f"table span            : [{lo:.0f}, {hi:.0f}]  ({G21.size} nodes)")
    print(f"GL span               : [0, {LF:.2f}]  (24 nodes)")
    print(f"GL nodes below {lo:.0f}   : {below.sum()} of 24 -> "
          f"{np.round(lam[below],2)}")
    print(f"1-D quadrature weight below {lo:.0f}: "
          f"{jw[below].sum()/jw.sum():.4f} of the length")
    print(f"2-D nodes with >=1 leg extrapolated: "
          f"{1 - ((~below).sum()/24)**2:.4f}")
    print()

    for tag, n2 in GAMMAS.items():
        C = np.array([[co.C_fn(N1, float(a), n2, float(b))[0, 0] for b in lam]
                      for a in lam])
        W = np.outer(jw, jw)
        contrib = W * C
        tot = contrib.sum()
        ext = below[:, None] | below[None, :]
        onlow = below[:, None] & below[None, :]
        print(f"--- gamma = {tag}   xi_00^O0 = {tot:.6e}")
        print(f"    from nodes with >=1 leg extrapolated (<406): "
              f"{contrib[ext].sum()/tot*100:+8.3f} %")
        print(f"    from nodes with BOTH legs extrapolated     : "
              f"{contrib[onlow].sum()/tot*100:+8.3f} %")
        print(f"    from the 24 strictly-diagonal nodes        : "
              f"{np.trace(contrib)/tot*100:+8.3f} %")
        # where does the mass sit relative to |lam1-lam2| ?
        d = np.abs(lam[:, None] - lam[None, :])
        for cut in (0.0, 100.0, 300.0, 1000.0):
            m = d <= cut if cut > 0 else np.eye(24, dtype=bool)
            print(f"      |lam1-lam2| <= {cut:6.0f}: "
                  f"{contrib[m].sum()/tot*100:7.3f} %  "
                  f"({m.sum():4d} nodes)")
        print()

    # what the extrapolation actually returns near lam -> 0
    print("=== C_00 on the equal-time diagonal, extrapolated region ===")
    print(f"{'lam':>9} {'C_00(lam,lam)':>15}   note")
    for t in (5.566, 29.2268, 71.3864, 131.3621, 208.1703, 300.5501,
              406.0, 406.9848, 550.0):
        v = co.C_fn(N1, float(t), N1, float(t))[0, 0]
        note = "EXTRAPOLATED" if t < G21[0] else ("table node" if t in G21 else "")
        print(f"{t:9.4f} {v:15.6e}   {note}")


if __name__ == "__main__":
    main()
