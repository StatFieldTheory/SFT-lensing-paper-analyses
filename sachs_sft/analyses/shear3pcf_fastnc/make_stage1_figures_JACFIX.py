"""
make_stage1_figures_JACFIX.py
Generate STAGE 1 synthesis figures for the shear 3PCF validation, JACFIX edition.
  - figures/stage1_threeway_gamma.pdf
  - figures/stage1_threeway_components.pdf
  - outputs/stage1_threeway_gamma.png
  - outputs/stage1_threeway_components.png
Reads from outputs/stage1_threeway_compare.npz (now the JACFIX threeway record)
and its source npz files.  The per-component figure reads the JACFIX canoes
output stage1_ours_v2_JACFIX.npz (the (1+z)^-4 Jacobian enters ONLY from the
canoes equal-shell vertex; no manual born3).
DO NOT recompute physics; read-only on all data files.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os

BASE = os.path.dirname(os.path.abspath(__file__))
FIG  = os.path.join(BASE, "figures")
OUT  = os.path.join(BASE, "outputs")

os.makedirs(FIG, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

# ── load data ────────────────────────────────────────────────────────────────
tw   = np.load(os.path.join(OUT, "stage1_threeway_compare.npz"), allow_pickle=True)
dref = np.load(os.path.join(OUT, "stage1_bkappa_phase_ref.npz"), allow_pickle=True)
dfnc = np.load(os.path.join(OUT, "stage1_fastnc.npz"),           allow_pickle=True)
dcan = np.load(os.path.join(OUT, "stage1_ours_v2_JACFIX.npz"),   allow_pickle=True)

gamma_tw = tw["gamma_arcmin"]   # shape (4,)
phi_tw   = tw["phi_deg"]        # shape (4,)
fi_ref   = tw["fi_ref"]         # (4,)  frame-invariant for ref
fi_fnc   = tw["fi_fnc"]         # (4,)
fi_can   = tw["fi_can"]         # (4,)
r_ref_fnc = tw["r_ref_fnc"]     # fi_fnc / fi_ref
r_ref_can = tw["r_ref_can"]     # fi_can / fi_ref

# config labels
labels = [rf"$\gamma={g:.0f}'$, $\phi={p:.0f}°$"
          for g, p in zip(gamma_tw, phi_tw)]

# ── FIGURE 1: frame-invariant Gamma (threeway) + ratio panels ────────────────
# Top panel: |Gamma_frame_inv| for the 3 pipelines
# Bottom panel: fi_can/fi_ref and fi_can/fi_fnc

x = np.arange(len(labels))
width = 0.25

fig, axes = plt.subplots(2, 1, figsize=(8, 6.5),
                         gridspec_kw={"height_ratios": [1.8, 1]})
fig.subplots_adjust(hspace=0.05)

# ---- top: absolute values ----
ax = axes[0]
bars_ref = ax.bar(x - width, fi_ref, width, label=r"$B_\kappa$-phase ref",
                  color="#4477AA", alpha=0.85, edgecolor="k", linewidth=0.5)
bars_fnc = ax.bar(x,          fi_fnc, width, label="fastnc tree",
                  color="#66CCEE", alpha=0.85, edgecolor="k", linewidth=0.5)
bars_can = ax.bar(x + width,  fi_can, width, label=r"canoes $\zeta_D$",
                  color="#EE6677", alpha=0.85, edgecolor="k", linewidth=0.5)

ax.set_yscale("log")
ax.set_ylabel(r"$\sqrt{\sum_\mu |\Gamma^\mu|^2}$", fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels([])
ax.legend(loc="upper right", fontsize=9)
ax.set_title(r"Three-way frame-invariant shear 3PCF: $B_\kappa$-phase ref vs fastnc vs canoes $\zeta_D$",
             fontsize=10, pad=6)
ax.tick_params(axis="y", which="both", labelsize=9)
ax.grid(axis="y", which="major", ls="--", lw=0.4, alpha=0.5)

# annotate canoes bars with the magnification factor
for i, (fc, fr) in enumerate(zip(fi_can, fi_ref)):
    factor = fc / fr
    ypos = fc * 1.05
    ax.text(x[i] + width, ypos, f"{factor:.1f}x",
            ha="center", va="bottom", fontsize=7.5, color="#CC2222",
            fontweight="bold")

# ---- bottom: ratio panels ----
ax2 = axes[1]
ratio_can_ref = fi_can / fi_ref   # JACFIX: ~0.3-1.4x (OOM agreement, median ~0.70x)
ratio_can_fnc = fi_can / fi_fnc   # similar

ax2.plot(x, ratio_can_ref, "o-", color="#EE6677", lw=1.6, ms=7,
         label=r"canoes / ref")
ax2.plot(x, ratio_can_fnc, "s--", color="#AA3355", lw=1.6, ms=7,
         label=r"canoes / fastnc")
ax2.plot(x, r_ref_fnc,     "D:",  color="#4477AA", lw=1.4, ms=6,
         label=r"fastnc / ref (self-check)")

ax2.axhline(1.0, color="k", lw=0.8, ls="--")
ax2.set_yscale("log")
ax2.set_ylabel("Ratio", fontsize=11)
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=8.5)
ax2.tick_params(axis="y", which="both", labelsize=9)
ax2.legend(loc="upper right", fontsize=8.5)
ax2.grid(axis="y", which="major", ls="--", lw=0.4, alpha=0.5)

# annotate median lines
med_cr = np.median(ratio_can_ref)
ax2.axhline(med_cr, color="#EE6677", lw=0.8, ls=":")
ax2.text(3.4, med_cr * 1.08,
         f"median {med_cr:.2f}x", color="#EE6677", fontsize=7.5, ha="right")

med_fc = np.median(r_ref_fnc)
ax2.axhline(med_fc, color="#4477AA", lw=0.8, ls=":")
ax2.text(3.4, med_fc * 0.88,
         f"median {med_fc:.2f}", color="#4477AA", fontsize=7.5, ha="right")

fig.tight_layout()
for ext, d in [(".pdf", FIG), (".png", OUT)]:
    fig.savefig(os.path.join(d, f"stage1_threeway_gamma{ext}"),
                dpi=200, bbox_inches="tight")
plt.close(fig)
print("Wrote stage1_threeway_gamma.pdf / .png")


# ── FIGURE 2: per-component bar chart at PRIMARY point (10', 60deg) ──────────
# PRIMARY point is config index 0  (gamma=10, phi=60)

# Extract per-component magnitudes for the three pipelines at (10', 60deg)
# ref: gamma_arcmin has shape (4,); index 0 = (10', 60deg)
comp_labels = [r"$|\Gamma^0|$", r"$|\Gamma^1|$",
               r"$|\Gamma^2|$", r"$|\Gamma^3|$"]
mu_idx = 0  # config index

ref_vals = np.array([np.abs(dref["Gamma0"][mu_idx]),
                     np.abs(dref["Gamma1"][mu_idx]),
                     np.abs(dref["Gamma2"][mu_idx]),
                     np.abs(dref["Gamma3"][mu_idx])])

# fastnc: find (10', 60deg) indices
ig_f = np.argmin(np.abs(dfnc["gamma_arcmin"] - 10))
ip_f = np.argmin(np.abs(dfnc["phi_deg"]   - 60))
fnc_vals = np.array([np.abs(dfnc["Gamma0"][ig_f, ip_f]),
                     np.abs(dfnc["Gamma1"][ig_f, ip_f]),
                     np.abs(dfnc["Gamma2"][ig_f, ip_f]),
                     np.abs(dfnc["Gamma3"][ig_f, ip_f])])

# canoes: same indices
ig_c = np.argmin(np.abs(dcan["gamma_arcmin"] - 10))
ip_c = np.argmin(np.abs(dcan["phi_deg"]   - 60))
can_vals = np.array([np.abs(dcan["Gamma0"][ig_c, ip_c]),
                     np.abs(dcan["Gamma1"][ig_c, ip_c]),
                     np.abs(dcan["Gamma2"][ig_c, ip_c]),
                     np.abs(dcan["Gamma3"][ig_c, ip_c])])

ncomp = 4
xc = np.arange(ncomp)
wc = 0.25

fig2, ax3 = plt.subplots(figsize=(7, 4.5))
ax3.bar(xc - wc, ref_vals, wc, label=r"$B_\kappa$-phase ref",
        color="#4477AA", alpha=0.85, edgecolor="k", linewidth=0.5)
ax3.bar(xc,       fnc_vals, wc, label="fastnc tree",
        color="#66CCEE", alpha=0.85, edgecolor="k", linewidth=0.5)
ax3.bar(xc + wc,  can_vals, wc, label=r"canoes $\zeta_D$",
        color="#EE6677", alpha=0.85, edgecolor="k", linewidth=0.5)

ax3.set_yscale("log")
ax3.set_ylabel(r"$|\Gamma^\mu|$", fontsize=12)
ax3.set_xticks(xc)
ax3.set_xticklabels(comp_labels, fontsize=11)
ax3.legend(fontsize=9)
ax3.set_title(r"Per-component $|\Gamma^\mu|$ at primary point "
              r"($\gamma=10'$, $\phi=60°$, $z_s=5$)", fontsize=10, pad=6)
ax3.grid(axis="y", which="major", ls="--", lw=0.4, alpha=0.5)

# annotate canoes/ref ratio per component
for i in range(ncomp):
    ratio_i = can_vals[i] / ref_vals[i]
    yp = max(can_vals[i], ref_vals[i]) * 1.12
    color = "#CC2222" if ratio_i > 1.5 else "#2266AA"
    ax3.text(xc[i] + wc, yp,
             f"can/ref={ratio_i:.2f}",
             ha="center", va="bottom", fontsize=7, color=color,
             fontweight="bold")

# annotate the order-of-magnitude agreement (post (1+z)^-4 Jacobian fix)
fi_can_pt = np.sqrt(np.sum(can_vals ** 2))
fi_ref_pt = np.sqrt(np.sum(ref_vals ** 2))
ax3.annotate(
    rf"canoes/ref (frame-inv) = {fi_can_pt / fi_ref_pt:.2f}$\times$"
    "\nOOM agreement (tree-level Limber)",
    xy=(xc[3] + wc, can_vals[3]),
    xytext=(xc[1] - 0.1, max(can_vals.max(), ref_vals.max()) * 0.18),
    fontsize=7.5, color="#2266AA",
    arrowprops=dict(arrowstyle="->", color="#2266AA", lw=0.9),
    ha="left",
)

fig2.tight_layout()
for ext, d in [(".pdf", FIG), (".png", OUT)]:
    fig2.savefig(os.path.join(d, f"stage1_threeway_components{ext}"),
                 dpi=200, bbox_inches="tight")
plt.close(fig2)
print("Wrote stage1_threeway_components.pdf / .png")
