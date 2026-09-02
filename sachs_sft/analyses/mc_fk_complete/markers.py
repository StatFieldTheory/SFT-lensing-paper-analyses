"""Marker values for the paper's validation figure, on its own GAMMA_MC grid."""
import numpy as np
d = np.load("_papergrid.npz", allow_pickle=True)
print(f"{'gamma':>7} {'MC s=8':>12} {'MC s=4':>12} {'sig->0':>12} {'+/-':>10} "
      f"{'analytic':>12} {'ratio':>7} {'+/-':>7} {'exact s->0/ana':>14}")
rows = []
for gam in sorted(set(d["gamma"])):
    i8 = int(np.where((d["gamma"] == gam) & (d["sigma"] == 8.0))[0][0])
    i4 = int(np.where((d["gamma"] == gam) & (d["sigma"] == 4.0))[0][0])
    s8 = np.asarray(d["seeds"][i8], float); s4 = np.asarray(d["seeds"][i4], float)
    ext = 2.0 * s4 - s8                      # linear in sigma, paired seeds
    ana = float(d["analytic"][i8])
    ex0 = 2.0 * d["exact"][i4] - d["exact"][i8]
    mu, se = ext.mean(), ext.std(ddof=1) / np.sqrt(ext.size)
    rows.append((gam, mu, se, ana))
    print(f"{gam:>7.2f} {s8.mean():>12.4e} {s4.mean():>12.4e} {mu:>12.4e} "
          f"{se:>10.2e} {ana:>12.4e} {mu/ana:>7.4f} {se/ana:>7.4f} "
          f"{ex0/ana:>14.4f}")
r = np.array(rows)
print(f"\nweighted mean ratio over 1'-17.3': "
      f"{np.average(r[:,1]/r[:,3], weights=(r[:,3]/r[:,2])**2):.4f}")
np.savez("_markers_papergrid.npz", gamma=r[:, 0], fk=r[:, 1], err=r[:, 2],
         analytic=r[:, 3])
print("-> _markers_papergrid.npz")
