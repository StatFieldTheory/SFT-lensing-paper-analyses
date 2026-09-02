"""FK as a fraction of Order-0 over 2'-12', the range current cosmic-shear
analyses actually use: xi_kappa at z_s = 5 and at the five multi-z planes.

Paths resolve relative to this file (the script moved into the analysis package
on 2026-09-02; it used to be run from the repository root with relative strings).
"""
import math
from pathlib import Path

import numpy as np

SACHS = Path(__file__).resolve().parents[2]
O0 = SACHS / "sftwick_outputs" / "2PCF" / "C_corr_op_O0" / "xi_C_corr_op_O0.npz"
FK = (SACHS / "sftwick_outputs" / "2PCF" / "C_corr_op_K_limber_FK_cut15360_permfix"
      / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz")
MULTIZ = SACHS / "analyses" / "analysis3" / "outputs" / "multiz_kappa_2pcf_5z.npz"
GS = np.array([2.0, 5.0, 8.0, 12.0])


def fold(path, order, a, b):
    # allow_pickle: the sweep npz files are this pipeline's own outputs, whose
    # direction vectors are stored as object arrays (every loader in the package
    # reads them this way); nothing here is downloaded or user supplied.
    d = np.load(path, allow_pickle=True)
    m = (d["a"] == a) & (d["b"] == b) & (d["order"] == order)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)), -1, 1)))) * 60)
    o = np.argsort(g)
    return np.asarray(g)[o], np.asarray(d["value"], float)[m][o]


print("=== z_s = 5, the main analysis: FK / Order-0 [%] ===")
go, o0 = fold(O0, 0, 0, 0)
gk, fk = fold(FK, 2, 0, 0)
print(f"{'gamma':>7} {'FK':>12} {'FK/O0 %':>9}")
r5 = []
for t in GS:
    a = float(np.interp(t, gk, fk)); b = float(np.interp(t, go, o0))
    r5.append(100 * a / b)
    print(f"{t:>7.1f} {a:>12.4e} {r5[-1]:>9.2f}")
print(f"  -> kappa-kappa over 2'-12': {min(r5):.2f}% to {max(r5):.2f}%")

print("\n=== multi-z, kappa-kappa: FK / Order-0 [%] at each source plane ===")
d = np.load(MULTIZ, allow_pickle=True)
z, g = np.asarray(d["z"], float), np.asarray(d["gamma"], float)
print(f"{'z_s':>5} {'2 arcmin':>10} {'5':>8} {'8':>8} {'12':>8}   {'range':>14}")
for i, zi in enumerate(z):
    r = [100 * float(np.interp(x, g, d["fk"][i])) / float(np.interp(x, g, d["o0"][i])) for x in GS]
    print(f"{zi:>5.1f} " + " ".join(f"{v:>8.2f}" for v in r) + f"   {min(r):.2f}-{max(r):.2f}%")
