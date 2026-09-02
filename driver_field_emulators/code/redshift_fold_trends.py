"""Source-redshift trends of the FK fold, measured from the dense tables.

The equal-time FK fold (fk_analytic.py eq. FK, validated against the
production sft-wick fold to 0.1%) integrates the vertex over shells up to
the source affine distance lambda_f. The dense tables cover shells to
2326.6 Mpc >= lambda_f(z_s = 5) = 2318, so the SAME table can be folded at
any lambda_f(z_s <= 5): the z_s trend is a measurement, not an estimate.

For the kk pair only zeta_TTT + zeta_Bmod enter (fk_analytic docstring:
G_00 = -2 (Z_000 + Z_011 + Z_022) = -2 (zeta_TTT + zeta_Bmod)), so this
script implements the fold directly on the table rows, with PCHIP-in-
lambda interpolation of sign * ln|G_00| (the callable's linear-in-lambda
interpolation is a systematic at low z_s where the shell grid built for
z_s = 5 is sparse; both interpolations are run and compared).

Measurements:
  T1  control: fold the tree table at the production t_final and compare
      with the archived production fold (expect the 0.1%-level agreement
      of the mc cross-check);
  T2  xi_FK(gamma; z_s) for tree (floor 50 Mpc), tree (floor 397) and
      BiHalofit (floor 397); the floor-matched NL/tree ratio; the
      lambda_min effect vs z_s (floor 397 / floor 50);
  T3  FK/O0 (gamma; z_s) with the exact production O0 from the talk's
      multiz_components.npz (validated against production at z = 5);
      the deployed-table FK trend from the same npz for comparison;
  T4  fold-weight profile over shells at gamma = 0.5', the fold-weighted
      cutoff undercount <U>(z_s) (per-shell U from redshift_pershell.py),
      and the closure check <rnl>_w vs the measured NL fold ratio.

Run (canoes venv, from driver_field_emulators/code):
    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python redshift_fold_trends.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline, PchipInterpolator

HERE = Path(__file__).resolve().parent
PRODUCTS = HERE.parent / "products"
OUT_DIR = PRODUCTS / "redshift"
OUT_DIR.mkdir(exist_ok=True)

_SCRIPTS = (HERE.parent.parent / "SFT-lensing-paper-analyses" / "sachs_sft"
            / "scripts")
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
import D_callable as Dc  # noqa: E402

from mc_crosscheck_dense import load_xi  # noqa: E402  (audited reader)

MULTIZ = (HERE.parent.parent / "talk" / "assets" / "figures" / "_data"
          / "multiz_components.npz")
T_FINAL_PROD = 2313.029          # production z_s = 5 baseline
ARCMIN = np.pi / (180.0 * 60.0)

# ---------------------------------------------------------------- tables


class CollapsedTable:
    """kk-channel collapsed-family vertex G_00(lambda; gamma) from a dense
    table, with per-row interpolation in lambda."""

    def __init__(self, path: Path):
        d = np.load(path)
        self.name = path.stem
        cts = np.asarray(d["cosine_triples"], float)
        self.g_rows = np.degrees(np.arccos(np.clip(cts[:, 1], -1, 1))) * 60.0
        self.lam = np.asarray(d["lambda_shells_Mpc"], float)
        self.g00 = -2.0 * (np.asarray(d["zeta_TTT"], float)
                           + np.asarray(d["zeta_Bmod"], float))

    def row(self, gamma_arcmin: float) -> int:
        return int(np.argmin(np.abs(self.g_rows - gamma_arcmin)))

    def interp(self, gamma_arcmin: float, kind: str = "logpchip"):
        """Return G_00(lambda) interpolant over the table's lambda span.

        Outside the shell span the vertex is taken as zero, matching the
        production callable. Rows with mixed signs fall back to direct
        PCHIP; single-sign rows use PCHIP of sign * ln|G_00|.
        """
        lam = self.lam
        y = self.g00[self.row(gamma_arcmin)]
        in_span = lambda x: ((np.asarray(x) >= lam[0])  # noqa: E731
                             & (np.asarray(x) <= lam[-1]))
        if kind == "linear":
            return lambda x: np.where(
                in_span(x), np.interp(np.asarray(x), lam, y), 0.0)
        if np.any(y == 0.0) or len(set(np.sign(y))) > 1:
            f = PchipInterpolator(lam, y)
            return lambda x: np.where(
                in_span(x), f(np.clip(x, lam[0], lam[-1])), 0.0)
        s = float(np.sign(y[0]))
        f = PchipInterpolator(lam, np.log(np.abs(y)))
        return lambda x: np.where(
            in_span(x), s * np.exp(f(np.clip(x, lam[0], lam[-1]))), 0.0)


# ------------------------------------------------------------------ fold


class FoldEngine:
    """Windows and quadrature of the equal-time FK fold at one lambda_f.

    Transcribes fk_analytic.fk_analytic (validated to 0.1% against the
    production sft-wick fold): xi_FK = INT dlam W u^4 Hresp G_00 with
    W(la) = INT_la^{lam_f} (D(la)/D(t))^2 dt, u = D/D(lam_f),
    Hresp(la) = INT_la^{lam_f} W u^-4 dt, D the Jacobi map a*chi.
    """

    def __init__(self, lam_f: float, lam_lo: float, n_lam: int = 96,
                 n_win: int = 64):
        self.lam_f = float(lam_f)
        self.lam_lo = float(lam_lo)
        dense = np.linspace(max(lam_lo, 1.0), lam_f, 600)
        D_dense = np.asarray(Dc.D_at(dense), float)
        Df = float(np.asarray(Dc.D_at(np.array([lam_f])), float)[0])
        u_dense = D_dense / Df
        xg, wg = leggauss(n_win)
        W = np.empty_like(dense)
        for i, la in enumerate(dense):
            t = 0.5 * (lam_f - la) * xg + 0.5 * (lam_f + la)
            jw = 0.5 * (lam_f - la) * wg
            W[i] = float(np.sum(jw * (D_dense[i]
                                      / np.asarray(Dc.D_at(t), float)) ** 2))
        self._W = CubicSpline(dense, W)
        h_spline = CubicSpline(dense, W * u_dense ** (-4))
        anti = h_spline.antiderivative()
        self._Hresp = lambda la: float(anti(lam_f)) - anti(la)
        self._u = CubicSpline(dense, u_dense)
        x, w = leggauss(n_lam)
        self.nodes = 0.5 * (lam_f - self.lam_lo) * x \
            + 0.5 * (lam_f + self.lam_lo)
        self.jac = 0.5 * (lam_f - self.lam_lo) * w
        self.weight_geom = np.asarray(
            [self._W(la) * self._u(la) ** 4 * self._Hresp(la)
             for la in self.nodes])

    def fold(self, g00_fn) -> float:
        return float(np.sum(self.jac * self.weight_geom * g00_fn(self.nodes)))


def main() -> None:
    tree = CollapsedTable(PRODUCTS / "collapsed_dense_tree_cut1000.npz")
    tree_lo = CollapsedTable(
        PRODUCTS / "collapsed_dense_tree_cut1000_lowshells.npz")
    bih = CollapsedTable(PRODUCTS / "collapsed_dense_bihalofit_cut1000.npz")

    mz = np.load(MULTIZ)
    z_all = np.asarray(mz["z"], float)
    lam_all = np.asarray(mz["lam"], float)
    g_mz = np.asarray(mz["gamma"], float)
    o0_mz = np.asarray(mz["o0"], float)
    fk_mz = np.asarray(mz["fk"], float)

    # z_s picks (exact multiz grid points; lambda_f from the same npz)
    picks = [0, 2, 5, 8, 11, 15, 19]
    z_list = z_all[picks]
    lam_list = lam_all[picks]

    # Gammas must exist in BOTH the dense-table rows AND the production /
    # multiz 40-gamma grid: the dense grid interleaves half-step rows, and
    # snapping a table row to the nearest production gamma silently pairs
    # values 11% apart in angle, which biased the first version of the T1
    # and T3 tables (review round of 2026-08-26).
    common = [g for g in g_mz
              if np.min(np.abs(tree.g_rows - g)) < 1e-3 * max(g, 1.0)]
    gam_act = [min(common, key=lambda c: abs(c - t))
               for t in (0.5, 4.0, 20.0, 60.0, 145.0)]
    print("[gammas] matched table<->production grid:",
          np.round(gam_act, 3))

    # ---- T1: control against the archived production fold ----------------
    print("[T1] control: this fold vs archived production fold "
          f"(tree, t_final = {T_FINAL_PROD}):")
    gam_prod, _cos, xi_prod, _tf = load_xi(
        PRODUCTS / "collapsed_dense_tree_cut1000_xi.npz")
    eng = FoldEngine(T_FINAL_PROD, tree.lam[0])
    for g in gam_act:
        i = int(np.argmin(np.abs(gam_prod - g)))
        mine = eng.fold(tree.interp(g))
        print(f"  gamma {g:7.2f}'   mine {mine:+.5e}   "
              f"prod {xi_prod[i]:+.5e}   ratio {mine / xi_prod[i]:7.4f}")

    # interpolation sensitivity at the lowest z_s
    eng05 = FoldEngine(lam_list[0], tree_lo.lam[0])
    for g in (0.5, 20.0):
        a = eng05.fold(tree_lo.interp(g, "logpchip"))
        b = eng05.fold(tree_lo.interp(g, "linear"))
        print(f"  [interp @ z_s=0.5] gamma {g}': logpchip {a:+.4e} "
              f"linear {b:+.4e}  ratio {a / b:.3f}")

    # ---- T2 + T3: trends -------------------------------------------------
    pershell = np.load(OUT_DIR / "pershell.npz")
    U_lam = PchipInterpolator(pershell["lam"], pershell["U"])

    res = {k: np.zeros((len(z_list), len(gam_act)))
           for k in ("fk_tree50", "fk_tree397", "fk_nl397")}
    Uw = np.zeros(len(z_list))
    rnl_w = np.zeros(len(z_list))
    for iz, (zs, lf) in enumerate(zip(z_list, lam_list)):
        e50 = FoldEngine(lf, tree_lo.lam[0])
        e397 = FoldEngine(lf, tree.lam[0])
        for ig, g in enumerate(gam_act):
            res["fk_tree50"][iz, ig] = e50.fold(tree_lo.interp(g))
            res["fk_tree397"][iz, ig] = e397.fold(tree.interp(g))
            res["fk_nl397"][iz, ig] = e397.fold(bih.interp(g))
        # weight profile at 0.5' on the floor-397 tree fold
        w = np.abs(e397.jac * e397.weight_geom
                   * tree.interp(gam_act[0])(e397.nodes))
        w /= w.sum()
        Uw[iz] = float(np.sum(w * U_lam(np.clip(
            e397.nodes, pershell["lam"][0], pershell["lam"][-1]))))
        rnl_nodes = (bih.interp(gam_act[0])(e397.nodes)
                     / tree.interp(gam_act[0])(e397.nodes))
        rnl_w[iz] = float(np.sum(w * rnl_nodes))
        if iz in (0, len(z_list) - 1):
            zl = np.asarray([float(Dc.z_of_lambda(x)) for x in e397.nodes])
            order = np.argsort(e397.nodes)
            cw = np.cumsum(w[order])
            for frac in (0.25, 0.5, 0.75):
                jq = int(np.searchsorted(cw, frac))
                print(f"  [weights z_s={zs:.2f}] {int(frac*100)}% of the "
                      f"0.5' FK fold lies below shell z = {zl[order][jq]:.2f}")

    print("\n[T2] xi_FK trends (kk):")
    print("   z_s   " + "".join(f"   g={g:>6.1f}'" for g in gam_act))
    for tag in ("fk_tree50", "fk_tree397", "fk_nl397"):
        print(f"  -- {tag}")
        for iz, zs in enumerate(z_list):
            print(f"  {zs:5.2f} " + "".join(
                f"  {res[tag][iz, ig]:+.2e}" for ig in range(len(gam_act))))
    print("\n  lambda_min effect (floor 397 / floor 50):")
    for iz, zs in enumerate(z_list):
        print(f"  z_s = {zs:5.2f} " + "".join(
            f"  {res['fk_tree397'][iz, ig] / res['fk_tree50'][iz, ig]:7.3f}"
            for ig in range(len(gam_act))))
    print("\n  NL / tree (floor-matched):")
    for iz, zs in enumerate(z_list):
        print(f"  z_s = {zs:5.2f} " + "".join(
            f"  {res['fk_nl397'][iz, ig] / res['fk_tree397'][iz, ig]:7.3f}"
            for ig in range(len(gam_act))))

    # ---- T3: against exact O0 and the deployed trend ---------------------
    print("\n[T3] FK/O0 (exact multiz O0; FK = tree floor-50 fold):")
    print("   z_s   " + "".join(f"   g={g:>6.1f}'" for g in gam_act))
    ratio_fo = np.zeros((len(z_list), len(gam_act)))
    for iz, (iz_all, zs) in enumerate(zip(picks, z_list)):
        for ig, g in enumerate(gam_act):
            jg = int(np.argmin(np.abs(g_mz - g)))
            ratio_fo[iz, ig] = res["fk_tree50"][iz, ig] / o0_mz[iz_all, jg]
        print(f"  {zs:5.2f} " + "".join(
            f"  {100*ratio_fo[iz, ig]:7.3f}%" for ig in range(len(gam_act))))
    print("  deployed-table FK/O0 from the multiz npz (artifact-affected):")
    for iz, iz_all in enumerate(picks):
        cells = []
        for g in gam_act:
            jg = int(np.argmin(np.abs(g_mz - g)))
            cells.append(f"  {100*fk_mz[iz_all, jg]/o0_mz[iz_all, jg]:7.3f}%")
        print(f"  {z_list[iz]:5.2f} " + "".join(cells))

    # normalized amplitude trends at 0.5'
    print("\n  amplitude trends normalized to z_s = 5 (gamma = 0.5'):")
    jg = int(np.argmin(np.abs(g_mz - gam_act[0])))
    o0_trend = o0_mz[picks, jg] / o0_mz[picks[-1], jg]
    fk_trend = res["fk_tree50"][:, 0] / res["fk_tree50"][-1, 0]
    fkd_trend = fk_mz[picks, jg] / fk_mz[picks[-1], jg]
    print("   z_s     O0      FK(tree)  FK(deployed)  FK/O0(tree)")
    for iz, zs in enumerate(z_list):
        print(f"  {zs:5.2f}  {o0_trend[iz]:7.4f}  {fk_trend[iz]:8.4f}  "
              f"{fkd_trend[iz]:10.4f}  {fk_trend[iz]/o0_trend[iz]:9.3f}")

    # low-z_s local slope of FK/O0 vs chi_s (power-law anchor comparison)
    chi_s = np.asarray([float(Dc.chi_of_lambda(x)) for x in lam_list])
    r0 = ratio_fo[:, 0]
    slope = np.gradient(np.log(r0), np.log(chi_s))
    print("\n  d ln(FK/O0) / d ln chi_s at 0.5':",
          " ".join(f"{s:.2f}" for s in slope))
    print("  (power-law low-z anchor: -(n+1); n is scale- and z-mixed here)")

    # ---- T4: fold-weighted undercount and NL closure ---------------------
    print("\n[T4] fold-weighted cutoff undercount and NL closure (0.5'):")
    for iz, zs in enumerate(z_list):
        meas = res["fk_nl397"][iz, 0] / res["fk_tree397"][iz, 0]
        print(f"  z_s = {zs:5.2f}   <U>_w = {Uw[iz]:5.2f}   "
              f"<rnl>_w = {rnl_w[iz]:6.3f}   measured NL fold = {meas:6.3f}")

    np.savez(OUT_DIR / "fold_trends.npz",
             z_list=z_list, lam_list=lam_list, gammas_arcmin=np.array(gam_act),
             **res, ratio_fk_o0=ratio_fo, U_w=Uw, rnl_w=rnl_w,
             o0_trend=o0_trend, fk_trend=fk_trend, fk_deployed_trend=fkd_trend)
    print(f"\n[saved] {OUT_DIR / 'fold_trends.npz'}")


if __name__ == "__main__":
    main()
