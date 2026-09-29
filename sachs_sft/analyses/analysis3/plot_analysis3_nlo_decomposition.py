"""Analysis 3 — Order-2 NLO decomposition of the lensing 2-point statistics:
compare Order-0 (O0) vs the two Order-2 post-Born channels FF and FK.

Per the L2 expansion (and the archived analysis_3):

    Full = Order-0 + Order-2 (FF + FK)
    O0   : Born 2-point                       (C_corr_op_O0)
    FF   : two F-vertex insertions  -> nonlinear-propagation correction
    FK   : one F + one K(kappa3) vertex -> 3-point-cumulant contribution

(KK vanishes at Order-2: six tilde-phi legs cannot pair in MSR.)

Unlike the archived run (which bundled O0+FF in one "base" npz), the sachs_sft
runs are THREE separate sweeps — each `vertex_types`-isolated:
    O0 : C_corr_op_O0                  -> order 0 rows
    FF : C_corr_op_K_limber_FF (["F"]) -> order 2 rows (both-F)
    FK : C_corr_op_K_limber_FK (["FK"])-> order 2 rows (F x K)
All three use the **corr_op** C propagator + identical geometry (40-pt gamma,
t_final=2313.029, n_gauss=24), so they sum coherently.

Four observables (kappa-kappa, xi_+, xi_-, kappa-gamma_t); style ported from
the archived analysis_3 via the self-contained _plot_style.py sibling.

Output: outputs/analysis3_nlo_O0_FF_FK.{png,pdf}
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
FF_NPZ = RUNS / "C_corr_op_K_limber_FF" / "xi_C_corr_op_K_limber_FF.npz"
# The manuscript's FK: ell_max = 15360 with the permutation-aware vertex (2026-08-26).
# The June cut1000 sweep in C_corr_op_K_limber_FK/ is superseded and about 4.5x low at 0.5'.
FK_NPZ = RUNS / "C_corr_op_K_limber_FK_cut15360_permfix" / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz"
OUT_STEM = HERE / "outputs" / "analysis3_nlo_O0_FF_FK"

sys.path.insert(0, str(HERE.parents[2] / "reproduce"))
from product_paths import resolve_product  # noqa: E402

O0_NPZ = resolve_product("order0", O0_NPZ)
FF_NPZ = resolve_product("ff", FF_NPZ)
FK_NPZ = resolve_product("fk", FK_NPZ)

# allow_pickle: x/y are object arrays of unit vectors in TRUSTED local sft-wick
# outputs we just generated; no untrusted-source deserialization.
LOAD_KW: dict[str, Any] = {"allow_pickle": True}

COLOR_O0, COLOR_FF, COLOR_FK = PALETTE[0], PALETTE[1], PALETTE[2]
MARKER_O0, MARKER_FF, MARKER_FK = "o", "^", "s"

# Shared y-range matched to analysis1 (the 3j-suppressed FK in the spin-2
# panels falls below this floor; that suppression is discussed in the text).
Y_LIM = (1e-9, 2e-3)
# The separation range over which FK is quoted in the text: the range current
# cosmic-shear analyses retain after their small-scale cuts (KiDS-Legacy and
# DES Y3 cut at ~2', HSC Y3 at 7', UNIONS at 12'). Shaded in every panel so the
# quoted range is visible rather than asserted.
QUOTE_RANGE_ARCMIN = (2.0, 12.0)

# Sachs-scalar -> observable SIGN CONVENTION.
# The sft-wick output legs are the RAW Sachs-scalar LoS integrals ∫X_i (sign not
# yet applied). The observables are defined (cosmology.tex eqs 'kappa from Dsachs1',
# 'gamma from Dsachs23') as:
#     kappa   = -∫X_1            (component 0, Phi00/Ricci)
#     gamma_+ = -∫X_2 , gamma_x = -∫X_3   (components 1,2, Psi0/Weyl)
#     gamma_t = -gamma_+         (tangential rotation into the great-circle frame)
# For a 2-point ⟨O_a O_b⟩ the per-leg (-1) signs cancel pairwise EXCEPT where an
# ODD number of gamma_t legs appears. Net effect on the visible 2PCFs:
#   - kk, xi_+, xi_-  : sign^2 = +1  -> UNCHANGED (raw == physical)
#   - kappa-gamma_t   : kappa(-1) x gamma_t(=-gamma_+ -> +∫X_2)  ->  = -(raw (0,1))
# So kappa_gamma_t carries weight -1 on (0,1); everything else +1.
OBSERVABLE_SPECS = (
    ("kappa-kappa", r"$\xi_\kappa(\gamma)\;:\;\langle\kappa\kappa\rangle$",
     [((0, 0), +1.0)]),
    ("xi_plus", r"$\xi_+(\gamma)\;:\;\langle\gamma_+\gamma_+\rangle+\langle\gamma_\times\gamma_\times\rangle$",
     [((1, 1), +1.0), ((2, 2), +1.0)]),
    ("xi_minus", r"$\xi_-(\gamma)\;:\;\langle\gamma_+\gamma_+\rangle-\langle\gamma_\times\gamma_\times\rangle$",
     [((1, 1), +1.0), ((2, 2), -1.0)]),
    ("kappa_gamma_t", r"$\xi_{\kappa\gamma_t}(\gamma)\;:\;\langle\kappa\gamma_t\rangle$",
     [((0, 1), -1.0)]),   # Sachs-scalar sign: kappa_gamma_t = -(raw (0,1))
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
    pairs = {(int(ai), int(bi)) for ai, bi in zip(a, b)}
    for (pa, pb) in pairs:
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
    for p in (O0_NPZ, FF_NPZ, FK_NPZ):
        if not p.exists():
            print(f"[analysis3] MISSING input: {p}")
            return 2

    g0, o0g = load_sweep_order(O0_NPZ, 0)
    gff, ffg = load_sweep_order(FF_NPZ, 2)
    gfk, fkg = load_sweep_order(FK_NPZ, 2)
    for gx, nm in ((gff, "FF"), (gfk, "FK")):
        if not np.allclose(g0, gx, rtol=1e-6):
            print(f"[analysis3] WARNING: {nm} gamma grid differs from O0")
    gamma = g0

    apply_rcparams()
    import matplotlib.pyplot as plt

    panel_rc = {"axes.labelsize": 18, "axes.titlesize": 16,
                "xtick.labelsize": 14.5, "ytick.labelsize": 14.5,
                "legend.fontsize": 14.5}
    saved = {k: plt.rcParams[k] for k in panel_rc}
    plt.rcParams.update(panel_rc)

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.4),
                             sharex=True, sharey="row")
    legend_handles = None
    print("\n=== Analysis 3 NLO breakdown (O0 / FF / FK) ===")
    for i, (name, title, combos) in enumerate(OBSERVABLE_SPECS):
        ax = axes[i // 2, i % 2]
        o0 = _combine(o0g, combos)
        ff = _combine(ffg, combos)
        fk = _combine(fkg, combos)
        if o0 is None or ff is None or fk is None:
            ax.text(0.5, 0.5, "(component_pair missing)", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])
            continue
        full = o0 + ff + fk
        print(f"{name}:")
        print(f"  O0 range = [{o0.min():+.3e}, {o0.max():+.3e}]  O0(0.5')={o0[0]:+.4e}")
        print(f"  FF range = [{ff.min():+.3e}, {ff.max():+.3e}]  FF(0.5')={ff[0]:+.4e}  FF/O0(0.5')={ff[0]/o0[0]:+.3e}")
        print(f"  FK range = [{fk.min():+.3e}, {fk.max():+.3e}]  FK(0.5')={fk[0]:+.4e}  FK/O0(0.5')={fk[0]/o0[0]:+.3e}")

        ax.axvspan(*QUOTE_RANGE_ARCMIN, color="0.85", alpha=0.55, lw=0, zorder=0)
        plot_signed_line(ax, gamma, full, label="Full (O0+FF+FK)")
        plot_signed_markers(ax, gamma, o0, color=COLOR_O0, marker=MARKER_O0, label="Order-0")
        plot_signed_markers(ax, gamma, ff, color=COLOR_FF, marker=MARKER_FF, label="Order-2 FF")
        # Beyond ~1 degree the FK multipole sum stops converging (the low- and
        # high-ell branches cancel); those points are drawn faint and are not
        # quoted anywhere. See the cutoff subsection of the insights section.
        plot_signed_markers(ax, gamma, fk, color=COLOR_FK, marker=MARKER_FK,
                            label="Order-2 FK", faint_outside=(None, 60.0))
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

    # sharey="row": set each row independently. Top row (κκ, ξ_+) is dominated by
    # the O0/FF/FK plateaus, so raise its floor; bottom row keeps the full range.
    axes[0, 0].set_ylim(1e-8, 2e-3)   # top row: raised ymin
    axes[1, 0].set_ylim(*Y_LIM)       # bottom row: full (1e-9, 2e-3)
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
