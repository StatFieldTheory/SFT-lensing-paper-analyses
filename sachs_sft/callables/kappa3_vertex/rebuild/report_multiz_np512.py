"""Ladder + band table: deployed vs rebuilt multi-z FK, six source planes."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

A3 = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3")
RUNS = A3.parents[1] / "sftwick_outputs" / "2PCF"
PIECES = A3.parents[1] / "callables" / "kappa3_vertex" / "rebuild" / "products" / "pieces_nphi512"
sys.path.insert(0, str(A3))
import plot_analysis3_cl_decomposition as Cl  # noqa: E402

def xi_kappa(path, order):
    g, grouped = Cl.load_sweep_order(Path(path), order)
    return g, Cl._combine(grouped, [((0, 0), +1.0)])

dep = np.load(A3 / "outputs" / "multiz_kappa_2pcf_5z.npz", allow_pickle=True)
gam, o0_5, fk_dep5, z5 = dep["gamma"], dep["o0"], dep["fk"], dep["z"]

new = np.load(PIECES / "multiz_fk_np512_6planes.npz", allow_pickle=True)
gam_n, fk_new = new["gamma"], new["fk"]
assert np.allclose(gam, gam_n), "gamma grid mismatch"

# z_s = 5 plane: single-plane sweeps, O0 from the production Order-0 run.
g0, o0_z5 = xi_kappa(RUNS / "C_corr_op_O0" / "xi_C_corr_op_O0.npz", 0)
_, fk_dep_z5 = xi_kappa(
    RUNS / "C_corr_op_K_limber_FK_cut15360_permfix"
    / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz", 2)
assert np.allclose(g0, gam)

O0 = np.vstack([o0_5, o0_z5])
FKD = np.vstack([fk_dep5, fk_dep_z5])
FKN = fk_new                      # six planes, last one is the z_s=5 baseline
Z = list(z5) + [5.0]
LAM = list(new["lam"])

band = (gam >= 2.0) & (gam <= 12.0)

def pct(fk, o0):
    return 100.0 * np.abs(fk) / o0

print(f"{'z_s':>5} {'lambda':>11} | {'0.5 dep':>8} {'0.5 new':>8} {'ratio':>6} "
      f"| {'band 2-12 deployed':>20} | {'band 2-12 rebuilt':>20}")
for i, z in enumerate(Z):
    rd, rn = pct(FKD[i], O0[i]), pct(FKN[i], O0[i])
    print(f"{z:5.2f} {LAM[i]:11.4f} | {rd[0]:7.3f}% {rn[0]:7.3f}% "
      f"{rn[0]/rd[0]:6.3f} | {rd[band].min():8.3f}-{rd[band].max():.3f}%   "
      f"| {rn[band].min():8.3f}-{rn[band].max():.3f}%")

print()
print("FK(0.5') absolute, deployed -> rebuilt:")
for i, z in enumerate(Z):
    print(f"  z_s={z:4.2f}: {FKD[i,0]:+.6e} -> {FKN[i,0]:+.6e}  "
          f"({FKN[i,0]/FKD[i,0]:.4f}x)")
