"""MOCK-UP ONLY: what the validation figure would look like with FK markers.

Writes into this folder, never into figures/.  Uses the paper's own style module
and its cached analytic channels, with the FK line rebound to the corrected
fold exactly as code/figures_corrected/make_val_figure.py does.
"""
import math, sys
from pathlib import Path
import numpy as np
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, str(_bootstrap.MC_DIR))
from _plot_style import (PALETTE, annotate_sign_legend, apply_rcparams,
                         clip_gamma_axis, plot_signed_line)
apply_rcparams()
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = dict(np.load(_bootstrap.MC_DIR / "outputs" / "appendix_mc_curve.npz",
                 allow_pickle=True))
g = np.asarray(D["g"], float); o0 = np.asarray(D["o0"], float)
ff = np.asarray(D["ff_mom"], float); fk3 = np.asarray(D["fk3"], float)
g_mc = np.asarray(D["g_mc"], float)
ffc = np.asarray(D["ffc_mc"], float); ffe = np.asarray(D["ffc_se"], float)
M = np.load("_markers_pooled.npz")
gk, vk, ek = M["gamma"], M["fk"], M["err"]

C_O0, C_FF, C_FK = PALETTE[0], PALETTE[1], PALETTE[2]
fig, ax = plt.subplots(figsize=(6.4, 4.8), constrained_layout=True)
_smk = dict(sign_marker_size=3.2, sign_marker_alpha=0.45)
plot_signed_line(ax, g, o0, color=C_O0, lw=1.7, alpha=0.95, label="O0 (workflow)", **_smk)
plot_signed_line(ax, g, ff, color=C_FF, lw=1.7, alpha=0.95, label="FF (workflow)", **_smk)
plot_signed_line(ax, g, fk3, color=C_FK, lw=1.7, alpha=0.95, label="FK (workflow)",
                 faint_outside=(None, 60.0), **_smk)
ax.errorbar(g_mc, np.abs(ffc), yerr=ffe, fmt="s", color=C_FF, ms=6.5,
            mfc=C_FF, mec=C_FF, capsize=2.5, lw=1.1, zorder=5, label="FF (MC)")
ax.errorbar(gk, np.abs(vk), yerr=ek, fmt="o", color=C_FK, ms=6.5,
            mfc=C_FK, mec=C_FK, capsize=2.5, lw=1.1, zorder=6, label="FK (MC)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$\gamma\;[\mathrm{arcmin}]$")
ax.set_ylabel(r"$\xi_\kappa(\gamma)\;:\;\langle\kappa\kappa\rangle$")
clip_gamma_axis(ax, gamma_min=float(g.min()))
ax.legend(loc="lower left", ncol=1)
annotate_sign_legend(ax)
fig.savefig("mock_val_figure.pdf"); fig.savefig("mock_val_figure.png", dpi=200)

lo, hi = ax.get_ylim()
dec = math.log10(hi / lo)
bb = ax.get_window_extent()
h_in = bb.height / fig.dpi
print(f"y-axis spans {dec:.2f} decades over {h_in:.2f} inch  "
      f"({h_in/dec:.3f} inch per decade)")
for lbl, r in (("FK, 1% agreement", 1.011), ("FK error bar, 0.5%", 1.005),
               ("FF, 20% offset", 1.20)):
    pts = math.log10(r) * (h_in / dec) * 72.0
    print(f"  {lbl:<22} -> {pts:6.2f} pt on the page   "
          f"(marker diameter is 6.5 pt)")
