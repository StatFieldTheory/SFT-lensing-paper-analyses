import sys
import numpy as np
sys.path.insert(0, "..")
import fk_expect_exact as ex

g = ex.build_grid(n_lambda=1000)
cosg = float(np.cos(np.deg2rad(1.0/60.0)))
lam = np.linspace(1290.0, 1330.0, 81)
s = np.array([ex._CORE._assemble_6x6(g.builder, cosg, float(l))[0,0] for l in lam])
print("fine scan lam 1290..1330 (0.5 Mpc), Sigma2_00:")
for i in range(0, 81, 2):
    print(f"   lam={lam[i]:8.2f}  S00={s[i]: .5e}")
print()
lam2 = np.linspace(1690.0, 1730.0, 41)
s2 = np.array([ex._CORE._assemble_6x6(g.builder, cosg, float(l))[0,0] for l in lam2])
print("fine scan lam 1690..1730 (1 Mpc), Sigma2_00:")
for i in range(0, 41, 2):
    print(f"   lam={lam2[i]:8.2f}  S00={s2[i]: .5e}")
print()
lamc = np.linspace(406.0, 2313.0, 200)
sc = np.array([ex._CORE._assemble_6x6(g.builder, cosg, float(l))[0,0] for l in lamc])
print("coarse profile over the whole ray, Sigma2_00 (every 10th of 200):")
for i in range(0, 200, 10):
    print(f"   lam={lamc[i]:8.1f}  S00={sc[i]: .5e}")
print("number of sign changes over the 200-point coarse scan:",
      int(np.sum(np.diff(np.sign(sc)) != 0)))
