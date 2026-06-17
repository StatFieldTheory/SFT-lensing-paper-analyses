"""STAGE 1 v2 comparison: ours_v2 (derived map) vs fastnc tree-level.

Interpreter: any python with numpy (e.g.
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python).

READ-ONLY on outputs/stage1_fastnc.npz (the parallel fastnc tree-level run).
Writes outputs/stage1_compare_v2.npz with per-component ratios + phases and the
frame-invariant ratio across the (gamma, phi) grid.  Reports the PRIMARY point
first.  NO fitting: this only measures closure.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_OUT = _HERE / "outputs"


def _frame_invariant(d) -> np.ndarray:
    return np.sqrt(sum(np.abs(d[k]) ** 2 for k in ["Gamma0", "Gamma1", "Gamma2", "Gamma3"]))


def main() -> int:
    ours = np.load(_OUT / "stage1_ours_v2.npz", allow_pickle=True)
    fnc = np.load(_OUT / "stage1_fastnc.npz", allow_pickle=True)
    g = ours["gamma_arcmin"]
    p = ours["phi_deg"]
    assert np.allclose(g, fnc["gamma_arcmin"]) and np.allclose(p, fnc["phi_deg"]), \
        "grid mismatch between ours_v2 and fastnc"

    comps = ["Gamma0", "Gamma1", "Gamma2", "Gamma3"]
    Iou = _frame_invariant(ours)
    Ifn = _frame_invariant(fnc)
    finv_ratio = Iou / Ifn

    # per-component magnitude ratio and phase difference (deg)
    mag_ratio = {k: np.abs(ours[k]) / np.abs(fnc[k]) for k in comps}
    phase_diff = {k: np.degrees(np.angle(ours[k] / fnc[k])) for k in comps}

    ig0 = int(np.argmin(np.abs(g - 10.0)))
    ip0 = int(np.argmin(np.abs(p - 60.0)))

    print("=" * 72)
    print("STAGE 1 v2 vs fastnc tree-level  (PRIMARY point gamma=10', phi=60deg)")
    print("=" * 72)
    for k in comps:
        ov = ours[k][ig0, ip0]
        fv = fnc[k][ig0, ip0]
        print(f"  {k}: ours={ov:+.4e}  fnc={fv:+.4e}  "
              f"|ratio|={abs(ov)/abs(fv):.3f}  dphase={np.degrees(np.angle(ov/fv)):+.1f}deg")
    print(f"  frame-invariant: ours={Iou[ig0, ip0]:.4e}  fnc={Ifn[ig0, ip0]:.4e}  "
          f"ratio={finv_ratio[ig0, ip0]:.3f}")

    print("\n--- frame-invariant ratio across grid (rows=gamma, cols=phi) ---")
    for ig in range(g.size):
        print(f"  g={g[ig]:6.1f}  " + " ".join(f"{finv_ratio[ig, ip]:6.2f}" for ip in range(p.size)))
    print(f"  range {finv_ratio.min():.2f} - {finv_ratio.max():.2f}, "
          f"median {np.median(finv_ratio):.2f}, scatter(std/mean) {finv_ratio.std()/finv_ratio.mean():.2f}")

    print("\n--- per-component |ratio| (median over grid) ---")
    for k in comps:
        print(f"  {k}: median |ratio| = {np.median(mag_ratio[k]):.2f}  "
              f"(range {mag_ratio[k].min():.2f}-{mag_ratio[k].max():.2f})")

    print("\n--- per-component phase diff (median over grid, deg) ---")
    for k in comps:
        print(f"  {k}: median dphase = {np.median(phase_diff[k]):+.1f}  "
              f"(range {phase_diff[k].min():+.1f} .. {phase_diff[k].max():+.1f})")

    closed = (0.5 <= np.median(finv_ratio) <= 2.0) and (finv_ratio.std() / finv_ratio.mean() < 0.5)
    verdict = "CLOSED (O(1), structured-residual-free)" if closed else \
        "NOT CLOSED (structured residual remains)"
    print(f"\nVERDICT: {verdict}")

    np.savez(
        _OUT / "stage1_compare_v2.npz",
        gamma_arcmin=g, phi_deg=p,
        frame_invariant_ours=Iou, frame_invariant_fnc=Ifn,
        frame_invariant_ratio=finv_ratio,
        **{f"mag_ratio_{k}": mag_ratio[k] for k in comps},
        **{f"phase_diff_deg_{k}": phase_diff[k] for k in comps},
        primary_ig=ig0, primary_ip=ip0,
        verdict=verdict,
        note=(
            "ours_v2 (derived great-circle->centroid map) vs fastnc tree. "
            "frame_invariant_ratio = sqrt(sum|G^mu|^2)_ours / _fnc. mag_ratio_* "
            "are phase-invariant per-component magnitude ratios; phase_diff_* are "
            "arg(ours/fnc) in deg. No fitting; closure measured only."
        ),
    )
    print(f"\nsaved -> {_OUT / 'stage1_compare_v2.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
