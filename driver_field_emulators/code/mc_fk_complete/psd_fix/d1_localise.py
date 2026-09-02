"""Where does the non-PSD structure enter?

Sigma2(lam) = d/dlam [ D^4 C_ab(t,t) ] / D^4, with C_ab from corr_op sampled on
n_nodes lambda points and each of the 9 entries splined INDEPENDENTLY.  Three
candidates: (a) the sampled C itself is not smooth, (b) the spline derivative
is a poor estimate of it, (c) the within/cross block assembly.
"""
import sys, numpy as np
sys.path.insert(0, "..")
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)

def eig6(builder, cosg, lam):
    return np.linalg.eigvalsh(core._assemble_6x6(builder, cosg, float(lam)))

print("=== 1. is the non-PSD in the WITHIN block, the CROSS block, or A-C? ===")
b = ds.Sigma2Builder(background=bg, apply_c0=False)
for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    print(f"  gamma={gam}'")
    for lam in (1295.6, 1303.2, 1307.0, 1308.9, 1312.7, 1320.0):
        A = b.matrix(1.0, lam); C = b.matrix(cosg, lam)
        eA = np.linalg.eigvalsh(A)
        eAmC = np.linalg.eigvalsh(A - C); eApC = np.linalg.eigvalsh(A + C)
        e6 = eig6(b, cosg, lam)
        print(f"    lam={lam:7.1f}  min eig A={eA.min():+.3e}  "
              f"A-C={eAmC.min():+.3e}  A+C={eApC.min():+.3e}  6x6={e6.min():+.3e}")

print()
print("=== 2. does it move with n_nodes?  (spline artefact <=> yes) ===")
cosg = float(np.cos(np.deg2rad(17.3 / 60.0)))
lam = np.linspace(406.0, 2313.0, 1000)
for n in (80, 120, 160, 240, 320, 480, 640):
    bb = ds.Sigma2Builder(background=bg, n_nodes=n, apply_c0=False)
    mn = np.array([eig6(bb, cosg, l).min() for l in lam])
    bad = mn < 0
    print(f"    n_nodes={n:4d} (spacing {1922/(n-1):5.2f} Mpc): "
          f"{bad.sum():3d}/1000 nodes non-PSD, most negative {mn.min():+.3e}, "
          f"min/median-positive = {mn.min()/np.median(mn[~bad]):+.4f}")
