#!/usr/bin/env python
"""Multi-faceted slices of the equal_time_limber kappa3 vertex.

Primary axis requested: angular difference gamma (the pairwise separation
cos gamma_ij = n_i . n_j that indexes the zeta table).

We sweep gamma along two physically interpretable one-parameter families that
both live on the native cosine grid:

  * SQUEEZED / isoceles : (cos12, cos23, cos31) = (1, cos g, cos g).
    Two directions coincide (apex), the third is separated by g. This is the
    family adjacent to the FK 2-point apex cos=(1,1,1).
  * EQUILATERAL          : (cos g, cos g, cos g).
    All three directions pairwise separated by g (exists only for g <= 120 deg,
    i.e. cos g >= -1/2).

For each family we plot all four driving-field channels (TTT, TTP, TPP, PPP),
overlaying:
  - raw grid nodes (markers): the callable queried EXACTLY at grid cosines, so
    dist~0 -> the table value with no cosine interpolation.
  - callable interpolation (lines): the callable at a dense gamma grid, i.e. the
    k=4 inverse-distance kNN surface that sft-wick actually consumes.

We also slice the OTHER axes:
  - vs redshift z (equiv. lambda shell) at fixed gamma,
  - 2D heatmap zeta(gamma, z),
  - an interpolation-fidelity zoom (callable vs nodes between two adjacent nodes).

Run in the PyCCL env:
  /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python slice_equal_time_limber.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

HERE = Path(__file__).resolve().parent
CALLABLE_DIR = HERE.parent
OUT = HERE
CHANNELS = ("zeta_TTT", "zeta_TTP", "zeta_TPP", "zeta_PPP")
CH_LABEL = {"zeta_TTT": r"$\zeta_{TTT}$ ($\phi^3$)",
            "zeta_TTP": r"$\zeta_{TTP}$",
            "zeta_TPP": r"$\zeta_{TPP}$",
            "zeta_PPP": r"$\zeta_{PPP}$ ($\psi^3$)"}
# index of each channel inside the (3,3,3) tensor (its representative E-mode slot)
CH_SLOT = {"zeta_TTT": (0, 0, 0), "zeta_TTP": (0, 0, 1),
           "zeta_TPP": (0, 1, 1), "zeta_PPP": (1, 1, 1)}
LINTHRESH = 1e-16  # below this is numerical noise (median |zeta| ~ 1e-18)


def load_callable():
    """Import the sibling callable module by path."""
    path = CALLABLE_DIR / "equal_time_limber_kappa3_callable.py"
    spec = importlib.util.spec_from_file_location("etl_callable", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["etl_callable"] = mod
    spec.loader.exec_module(mod)
    return mod


def squeezed_dirs(gamma_rad):
    """Apex pair coincident along z; third direction at separation gamma."""
    n0 = np.array([0.0, 0.0, 1.0])
    n2 = np.array([np.sin(gamma_rad), 0.0, np.cos(gamma_rad)])
    return [n0, n0, n2]


def equilateral_dirs(gamma_rad):
    """Three directions pairwise separated by gamma (cone, 120 deg azimuths).

    pairwise cos = cos^2(theta) - 1/2 sin^2(theta) = cos gamma
      => cos^2 theta = (cos gamma + 1/2) / (3/2).
    Defined only for cos gamma >= -1/2 (gamma <= 120 deg).
    """
    cg = np.cos(gamma_rad)
    cos2 = (cg + 0.5) / 1.5
    cos2 = np.clip(cos2, 0.0, 1.0)
    ct = np.sqrt(cos2)
    st = np.sqrt(max(0.0, 1.0 - cos2))
    out = []
    for az in (0.0, 2 * np.pi / 3, 4 * np.pi / 3):
        out.append(np.array([st * np.cos(az), st * np.sin(az), ct]))
    return out


def eval_family(mod, dirs_fn, gammas_rad, lam):
    """coupling_fn over a gamma sweep at fixed equal-time lambda -> dict[ch]->array."""
    vals = {ch: np.full(gammas_rad.size, np.nan) for ch in CHANNELS}
    for i, g in enumerate(gammas_rad):
        K = mod.coupling_fn(dirs_fn(g), [lam, lam, lam])
        for ch in CHANNELS:
            vals[ch][i] = K[CH_SLOT[ch]]
    return vals


def grid_nodes_cos(npz):
    return np.unique(np.round(np.asarray(npz["cosine_triples"], float), 6))


# ---------------------------------------------------------------------------
def main():
    mod = load_callable()
    npz = np.load(mod.TABLE_PATH, allow_pickle=False)
    lam_shells = np.asarray(npz["lambda_shells_Mpc"], float)
    z_shells = np.asarray(npz["z_shells"], float)
    cos_nodes = grid_nodes_cos(npz)            # 15 cosine nodes, -1..1

    # representative shells for the gamma-slice curves
    z_pick = [0.45, 1.0, 2.5, 5.0]
    shell_idx = [int(np.argmin(np.abs(z_shells - z))) for z in z_pick]
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, len(shell_idx)))

    dense_g_deg = np.linspace(0.5, 180.0, 400)
    dense_g = np.deg2rad(dense_g_deg)

    # ---- node gammas for each family (markers land exactly on grid) --------
    node_g_sq = np.arccos(np.clip(cos_nodes, -1, 1))            # 0..180 deg
    node_g_eq = np.arccos(np.clip(cos_nodes[cos_nodes >= -0.5], -1, 1))  # <=120

    families = {
        "squeezed": dict(dirs=squeezed_dirs, node_g=node_g_sq,
                         title="Squeezed / isoceles  (1, cos g, cos g)",
                         gmax=180.0),
        "equilateral": dict(dirs=equilateral_dirs, node_g=node_g_eq,
                            title="Equilateral  (cos g, cos g, cos g)",
                            gmax=120.0),
    }

    figdata = {}
    FLOOR = 1e-18  # log-|zeta| floor (median |zeta| ~ 1e-18 is numerical noise)

    # =====================================================================
    # FIGURE 1 (headline): log|zeta_TTT| vs gamma, squeezed family, ALL
    # z-shells colored by viridis. Exposes the contact spike at gamma=0
    # (apex cos=1,1,1) sitting ~2 decades above the wide-separation plateau,
    # and the monotone growth toward high z.
    # =====================================================================
    zcolors = plt.cm.viridis(np.linspace(0.0, 0.95, z_shells.size))
    fig, ax = plt.subplots(figsize=(11, 7))
    for c, iz in zip(zcolors, range(z_shells.size)):
        lam = lam_shells[iz]
        interp = eval_family(mod, squeezed_dirs, dense_g, lam)["zeta_TTT"]
        ax.plot(dense_g_deg, np.abs(interp) + FLOOR, "-", color=c, lw=1.3)
        node = eval_family(mod, squeezed_dirs, node_g_sq, lam)["zeta_TTT"]
        ax.plot(np.rad2deg(node_g_sq), np.abs(node) + FLOOR, "o", color=c, ms=3)
    ax.set_yscale("log")
    ax.set_xlabel(r"angular separation $\gamma$ [deg]")
    ax.set_ylabel(r"$|\zeta_{TTT}|$  (line=kNN interp, $\circ$=grid node)")
    ax.set_title("equal_time_limber kappa3 vertex: contact spike in gamma "
                 "(squeezed, color=z 0.1->5.7)")
    ax.axvspan(0, 5, color="0.85", alpha=0.5, zorder=0)
    ax.annotate("collapsed apex\ncos=(1,1,1)", xy=(2, 1e-11), xytext=(25, 1e-11),
                fontsize=9, va="center",
                arrowprops=dict(arrowstyle="->", color="0.4"))
    sm = plt.cm.ScalarMappable(cmap="viridis",
                               norm=plt.Normalize(z_shells[0], z_shells[-1]))
    fig.colorbar(sm, ax=ax, label="source redshift z")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    p = OUT / "slice_gamma_TTT_contact.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"[saved] {p}")

    # =====================================================================
    # FIGURE 2: all four channels vs gamma at the highest-signal shell
    # (z=5.7). Squeezed (solid) vs equilateral (dashed) overlaid, so the
    # config-dependence of each driving-field cumulant is visible.
    # =====================================================================
    iz_hi = z_shells.size - 1
    lam_hi = lam_shells[iz_hi]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    for ax, ch in zip(axes.flat, CHANNELS):
        for fam_key, fam, ls in (("squeezed", families["squeezed"], "-"),
                                 ("equilateral", families["equilateral"], "--")):
            mask = dense_g_deg <= fam["gmax"] + 1e-9
            gd, gd_deg = dense_g[mask], dense_g_deg[mask]
            interp = eval_family(mod, fam["dirs"], gd, lam_hi)[ch]
            ax.plot(gd_deg, np.abs(interp) + FLOOR, ls, lw=1.5,
                    color="C0" if fam_key == "squeezed" else "C3",
                    label=fam_key)
            node = eval_family(mod, fam["dirs"], fam["node_g"], lam_hi)[ch]
            ax.plot(np.rad2deg(fam["node_g"]), np.abs(node) + FLOOR, "o",
                    ms=4, color="C0" if fam_key == "squeezed" else "C3",
                    mfc="white", mew=1.0)
            figdata[f"ch_{ch}_{fam_key}_zhi"] = np.abs(interp)
        ax.set_yscale("log")
        ax.set_title(CH_LABEL[ch])
        ax.grid(alpha=0.25, which="both")
        if ch == "zeta_TTT":
            ax.legend(fontsize=9, title=f"config (z={z_shells[iz_hi]:.1f})")
    for ax in axes[-1]:
        ax.set_xlabel(r"angular separation $\gamma$ [deg]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$|\zeta|$")
    fig.suptitle("equal_time_limber kappa3 vertex: per-channel gamma slice "
                 f"at z={z_shells[iz_hi]:.1f} (squeezed vs equilateral)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p = OUT / "slice_gamma_channels_zhi.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"[saved] {p}")

    # =====================================================================
    # FIGURE 3: zeta vs redshift z at fixed gamma (squeezed family).
    # =====================================================================
    g_pick_deg = [20.0, 60.0, 90.0, 135.0]
    gcolors = plt.cm.plasma(np.linspace(0.0, 0.8, len(g_pick_deg)))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    for ax, ch in zip(axes.flat, CHANNELS):
        for c, gdeg in zip(gcolors, g_pick_deg):
            g = np.deg2rad(gdeg)
            curve = np.array([mod.coupling_fn(squeezed_dirs(g),
                              [lam, lam, lam])[CH_SLOT[ch]] for lam in lam_shells])
            ax.plot(z_shells, curve, "-o", color=c, ms=4, lw=1.5,
                    label=f"g={gdeg:.0f} deg")
            figdata[f"zslice_{ch}_g{gdeg:.0f}"] = curve
        ax.axhline(0, color="0.7", lw=0.6, zorder=0)
        ax.set_yscale("symlog", linthresh=LINTHRESH)
        ax.set_title(CH_LABEL[ch])
        ax.grid(alpha=0.25)
        if ch == "zeta_TTT":
            ax.legend(fontsize=8, title="squeezed sep.")
    for ax in axes[-1]:
        ax.set_xlabel("source-plane redshift  z")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\zeta$ (equal-time)")
    fig.suptitle("equal_time_limber kappa3 vertex  --  redshift slice "
                 "(squeezed config, fixed gamma)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p = OUT / "slice_redshift_squeezed.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"[saved] {p}")

    # =====================================================================
    # FIGURE 4: 2D heatmap zeta(gamma, z) for the squeezed family.
    # =====================================================================
    g_heat_deg = np.linspace(2.0, 178.0, 120)
    g_heat = np.deg2rad(g_heat_deg)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, ch in zip(axes.flat, CHANNELS):
        Z = np.zeros((z_shells.size, g_heat.size))
        for iz, lam in enumerate(lam_shells):
            row = np.array([mod.coupling_fn(squeezed_dirs(g),
                            [lam, lam, lam])[CH_SLOT[ch]] for g in g_heat])
            Z[iz] = row
        vmax = np.nanmax(np.abs(Z))
        norm = SymLogNorm(linthresh=max(LINTHRESH, vmax * 1e-4),
                          vmin=-vmax, vmax=vmax)
        im = ax.pcolormesh(g_heat_deg, z_shells, Z, cmap="RdBu_r",
                           norm=norm, shading="nearest")
        ax.set_title(CH_LABEL[ch])
        ax.set_xlabel(r"$\gamma$ [deg]")
        ax.set_ylabel("z")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        figdata[f"heat_{ch}"] = Z
    fig.suptitle("equal_time_limber kappa3 vertex  --  zeta(gamma, z) "
                 "(squeezed; RdBu symlog)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p = OUT / "heatmap_gamma_z_squeezed.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"[saved] {p}")

    # =====================================================================
    # FIGURE 5: interpolation fidelity -- callable vs grid nodes, TTT,
    # zoomed on the small-gamma region where the signal is largest.
    # =====================================================================
    fig, ax = plt.subplots(figsize=(10, 6))
    si = z_shells.size - 1   # highest-signal shell: contact spike is real here
    lam = lam_shells[si]
    fine_deg = np.linspace(0.2, 90.0, 800)
    fine = np.deg2rad(fine_deg)
    interp = eval_family(mod, squeezed_dirs, fine, lam)["zeta_TTT"]
    ax.plot(fine_deg, interp, "-", color="C0", lw=1.5,
            label="callable kNN interp (k=4 inv-dist)")
    nmask = node_g_sq <= np.deg2rad(90)
    node = eval_family(mod, squeezed_dirs, node_g_sq[nmask], lam)["zeta_TTT"]
    ax.plot(np.rad2deg(node_g_sq[nmask]), node, "s", color="C3", ms=7,
            mfc="white", mew=1.5, label="grid nodes (exact table)")
    for gnode in np.rad2deg(node_g_sq[nmask]):
        ax.axvline(gnode, color="0.85", lw=0.6, zorder=0)
    ax.set_xlabel(r"angular separation $\gamma$ [deg]")
    ax.set_ylabel(r"$\zeta_{TTT}$")
    ax.set_title(f"Interpolation fidelity (squeezed, z={z_shells[si]:.1f}): "
                 f"kNN surface vs native cosine nodes -- note the contact spike "
                 f"at the apex")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    p = OUT / "interp_fidelity_TTT.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"[saved] {p}")

    np.savez(OUT / "slice_figdata.npz",
             lam_shells=lam_shells, z_shells=z_shells, cos_nodes=cos_nodes,
             **{k: np.asarray(v) for k, v in figdata.items()})
    print(f"[saved] {OUT / 'slice_figdata.npz'}")

    # ---- concise numeric summary printed to stdout ----------------------
    print("\n=== NUMERIC SUMMARY (squeezed, raw grid nodes) ===")
    print("zeta_TTT[10^-12] over (gamma_node, z_shell):")
    hdr = "  g\\z  " + "".join(f"{z:7.2f}" for z in z_shells[shell_idx])
    print(hdr)
    for gnode in node_g_sq:
        vals = [mod.coupling_fn(squeezed_dirs(gnode),
                [lam_shells[si]] * 3)[CH_SLOT["zeta_TTT"]] * 1e12
                for si in shell_idx]
        print(f"  {np.rad2deg(gnode):5.0f} " + "".join(f"{v:7.3f}" for v in vals))


if __name__ == "__main__":
    main()
