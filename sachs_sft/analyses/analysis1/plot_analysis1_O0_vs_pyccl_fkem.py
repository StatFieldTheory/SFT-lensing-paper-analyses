"""Analysis 1 — lensing Order-0 2-point statistics: SFT corr_op (full
non-Limber) vs PyCCL **FKEM non-Limber**.

Four observables (kappa-kappa, xi_+, xi_-, kappa-gamma_t), one 2x2 grid, NO
ratio panels -- styled exactly like analysis3 via the shared _plot_style.py.

Per panel, two curves of the SAME observable:
  * "PyCCL FKEM"  — thick gray sign-aware BACKGROUND line. C_L^{kk} from a
                    z_s=5 CMBLensingTracer computed with the FKEM non-Limber
                    integrator (ell<250) + Limber (ell>=250, exact there),
                    then summed full-sky with P_L / Wigner-d^L_{m,m'}.
  * "SFT O0"      — sft-wick FF Order-0 with the corr_op (full non-Limber)
                    C propagator. Coloured sign-aware markers.

Both sides are now NON-LIMBER, so this is an apples-to-apples validation: the
SFT corr_op Order-0 should track PyCCL FKEM across ALL gamma (the earlier
Limber reference only matched at small gamma / high ell). Residual departures
appear only at large gamma where each observable crosses zero (ratio
ill-defined) and on the tiny negative tail.

Inputs (consumed, not built here):
  ../../sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz   (SFT, non-Limber)
  ./pyccl_xi_shear_reference_fkem.npz                           (PyCCL FKEM)

Output:
  outputs/analysis1_O0_vs_pyccl_fkem.{png,pdf}
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # _plot_style sibling

from _plot_style import (  # noqa: E402
    GAMMA_X_MAX_ARCMIN,
    PALETTE,
    annotate_sign_legend,
    apply_rcparams,
    clip_gamma_axis,
    plot_signed_line,
    plot_signed_markers,
)

RUNS = HERE.parent.parent / "sftwick_outputs" / "2PCF"
O0_NPZ = RUNS / "C_corr_op_O0" / "xi_C_corr_op_O0.npz"
PYCCL_REF = HERE / "pyccl_xi_shear_reference_fkem.npz"
OUT_STEM = HERE / "outputs" / "analysis1_O0_vs_pyccl_fkem"

COLOR_SFT, MARKER_SFT = PALETTE[0], "o"

# Shared y-range across all four panels (data span within the gamma x-clip).
Y_LIM = (1e-9, 2e-3)

# allow_pickle: x/y are object arrays of unit vectors in TRUSTED local sft-wick
# outputs we generated; no untrusted-source deserialization.
LOAD_KW: dict[str, Any] = {"allow_pickle": True}

# (name, title, Sachs-scalar component combos, PyCCL-reference key).
# Component->observable signs follow analysis3 (cosmology.tex 'kappa from
# Dsachs1' / 'gamma from Dsachs23'): kappa=-int X1, gamma_+=-int X2. Pairwise
# (-1) signs cancel except for an odd count of gamma_t legs, so kk/xi_+/xi_-
# carry +1 and kappa_gamma_t carries -1 on the raw (0,1) leg.
OBSERVABLE_SPECS = (
    ("kappa-kappa", r"$\xi_\kappa(\gamma)\;:\;\langle\kappa\kappa\rangle$",
     [((0, 0), +1.0)], "xi_kappa"),
    ("xi_plus", r"$\xi_+(\gamma)\;:\;\langle\gamma_+\gamma_+\rangle+\langle\gamma_\times\gamma_\times\rangle$",
     [((1, 1), +1.0), ((2, 2), +1.0)], "xi_plus"),
    ("xi_minus", r"$\xi_-(\gamma)\;:\;\langle\gamma_+\gamma_+\rangle-\langle\gamma_\times\gamma_\times\rangle$",
     [((1, 1), +1.0), ((2, 2), -1.0)], "xi_minus"),
    ("kappa_gamma_t", r"$\xi_{\kappa\gamma_t}(\gamma)\;:\;\langle\kappa\gamma_t\rangle$",
     [((0, 1), -1.0)], "xi_kappa_gamma"),
)


def _gamma_arcmin(xi, yi):
    xi = np.asarray(xi, float) / np.linalg.norm(xi)
    yi = np.asarray(yi, float) / np.linalg.norm(yi)
    return np.degrees(np.arccos(np.clip(np.dot(xi, yi), -1, 1))) * 60.0


def load_sweep_order(path: Path, order: int):
    """Return (gamma_sorted, {(a,b): value_sorted}) for rows at `order`."""
    with np.load(path, **LOAD_KW) as d:
        a, b, o = np.asarray(d["a"]), np.asarray(d["b"]), np.asarray(d["order"])
        v, x, y = np.asarray(d["value"], float), d["x"], d["y"]
    grouped: dict[tuple[int, int], np.ndarray] = {}
    g_ref = None
    for (pa, pb) in sorted({(int(ai), int(bi)) for ai, bi in zip(a, b)}):
        m = (a == pa) & (b == pb) & (o == order)
        if not m.any():
            continue
        idx = np.flatnonzero(m)
        g = np.array([_gamma_arcmin(x[i], y[i]) for i in idx])
        s = np.argsort(g)
        if g_ref is None:
            g_ref = g[s]
        grouped[(pa, pb)] = v[m][s]
    if g_ref is None:
        raise ValueError(f"no rows at order={order} in {path}")
    return g_ref, grouped


def _combine(grouped, combos):
    out = None
    for (pa, pb), sign in combos:
        if (pa, pb) not in grouped:
            return None
        term = sign * grouped[(pa, pb)]
        out = term if out is None else out + term
    return out


def main() -> int:
    for p in (O0_NPZ, PYCCL_REF):
        if not p.exists():
            print(f"[analysis1] MISSING required input: {p}")
            return 2

    gamma, o0g = load_sweep_order(O0_NPZ, 0)
    ref = np.load(PYCCL_REF)
    g_pc = np.asarray(ref["gamma_arcmin"], float)
    if not np.allclose(gamma, g_pc, rtol=1e-6):
        print("[analysis1] NOTE: gamma grids differ; interpolating PyCCL onto SFT grid")

    apply_rcparams()
    import matplotlib.pyplot as plt

    panel_rc = {"axes.labelsize": 18, "axes.titlesize": 16,
                "xtick.labelsize": 14.5, "ytick.labelsize": 14.5,
                "legend.fontsize": 14.5}
    saved = {k: plt.rcParams[k] for k in panel_rc}
    plt.rcParams.update(panel_rc)

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.4),
                             sharex=True, sharey=True)
    legend_handles = None
    print("\n=== Analysis 1 Order-0: SFT corr_op (non-Limber) vs PyCCL FKEM ===")
    for i, (name, title, combos, refkey) in enumerate(OBSERVABLE_SPECS):
        ax = axes[i // 2, i % 2]
        sft = _combine(o0g, combos)
        pc = np.interp(gamma, g_pc, np.asarray(ref[refkey], float))
        if sft is None:
            ax.text(0.5, 0.5, "(component_pair missing)", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(pc != 0, sft / pc, np.nan)
        print(f"{name}:")
        print(f"  SFT  range = [{sft.min():+.3e}, {sft.max():+.3e}]  SFT(0.5')={sft[0]:+.4e}")
        print(f"  PyCCL range = [{pc.min():+.3e}, {pc.max():+.3e}]  PyCCL(0.5')={pc[0]:+.4e}")
        print(f"  SFT/PyCCL @ small gamma = {r[0]:.4f}, {r[1]:.4f}, {r[2]:.4f}")

        plot_signed_line(ax, gamma, pc, label="PyCCL")
        plot_signed_markers(ax, gamma, sft, color=COLOR_SFT, marker=MARKER_SFT,
                            label="SFT O0")
        ax.set_title(title)
        if i // 2 == 1:  # bottom row only (x shared)
            ax.set_xlabel(r"$\gamma\;[\mathrm{arcmin}]$")
        clip_gamma_axis(ax, gamma_min=gamma.min(), gamma_max=GAMMA_X_MAX_ARCMIN)
        if legend_handles is None:
            handles, labels = ax.get_legend_handles_labels()
            seen, uh, ul = set(), [], []
            for h, l in zip(handles, labels):
                if l not in seen:
                    seen.add(l); uh.append(h); ul.append(l)
            legend_handles = (uh, ul)
            annotate_sign_legend(ax)

    axes[0, 0].set_ylim(*Y_LIM)  # sharey -> applies to all panels
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    if legend_handles is not None:
        h, l = legend_handles
        fig.legend(h, l, loc="upper center", ncol=len(h),
                   bbox_to_anchor=(0.5, 0.995), frameon=False, fontsize=15.5)
    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf"):
        fig.savefig(OUT_STEM.with_suffix(ext))
    plt.close(fig)
    plt.rcParams.update(saved)
    print(f"\n-> {OUT_STEM.with_suffix('.png')}")
    print(f"-> {OUT_STEM.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
