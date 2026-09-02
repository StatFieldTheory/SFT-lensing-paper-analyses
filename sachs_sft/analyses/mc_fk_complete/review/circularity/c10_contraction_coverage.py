"""How much of the F-K index contraction does the kappa-kappa test exercise?

FK_kk contracts F_{0bc} zeta_{bc3} + transpose.  F row 0 is -delta_bc, so the
sum runs over c = 0 (scalar-scalar-scalar) and c = 1,2 (the two spin-2 legs).
Anything in F's spin-2 rows (F_101 = F_202 = -2) never appears.
"""
import _env, numpy as np
import fk_expect_exact as ex

N = 1000
grid = ex.build_grid(n_lambda=N)
kern = grid.dlam * grid.Wd * grid.Hd
print(f"{'gamma':>8}{'c=0 (TTT)':>14}{'c=1 (RePsi)':>14}{'c=2 (ImPsi)':>14}"
      f"{'spin share':>12}")
for gam in (0.5, 1.0, 5.0, 17.3):
    cg = float(np.cos(np.deg2rad(gam / 60.0)))
    z = np.array([ex._DS.zeta6(cg, float(l)) for l in grid.lam])
    parts = []
    for c in range(3):
        v = -np.sum(kern * z[:, c, c, 3]) * 2.0   # F_{0cc} = -1, plus transpose
        parts.append(v)
    tot = sum(parts)
    print(f"{gam:>7.1f}'{parts[0]:>14.4e}{parts[1]:>14.4e}{parts[2]:>14.4e}"
          f"{(parts[1]+parts[2])/tot:>12.4f}")
print("\nF components reachable by the (0,0) entry:",
      sorted({(0, b, c) for b in range(3) for c in range(3) if ex.F[0, b, c]}))
print("F components NEVER reached:",
      sorted({(a, b, c) for a in range(1, 3) for b in range(3)
              for c in range(3) if ex.F[a, b, c]}))
