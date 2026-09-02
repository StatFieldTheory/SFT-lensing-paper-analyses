"""Pool the four independent seed blocks; error from block-to-block scatter."""
import numpy as np
BLK = ["_blk6010001.npz", "_blk6020002.npz", "_blk6030003.npz", "_blk6040004.npz"]
D = [np.load(f, allow_pickle=True) for f in BLK]
gams = sorted(set(D[0]["gamma"]))
print(f"{'gamma':>7} {'blocks (sigma->0, /fold)':>44} {'mean':>8} {'+/-':>7} "
      f"{'FK value':>12} {'+/-':>10}")
out = []
for gam in gams:
    per, ana = [], None
    for d in D:
        i8 = int(np.where((d["gamma"] == gam) & (d["sigma"] == 8.0))[0][0])
        i4 = int(np.where((d["gamma"] == gam) & (d["sigma"] == 4.0))[0][0])
        ana = float(d["analytic"][i8])
        s8 = np.asarray(d["seeds"][i8], float); s4 = np.asarray(d["seeds"][i4], float)
        per.append((2.0 * s4 - s8).mean())        # linear sigma->0, this block
    per = np.array(per)
    m, se = per.mean(), per.std(ddof=1) / np.sqrt(per.size)   # BLOCK scatter
    out.append((gam, m, se, ana))
    print(f"{gam:>7.2f} {'  '.join(f'{x/ana:.4f}' for x in per):>44} "
          f"{m/ana:>8.4f} {se/ana:>7.4f} {m:>12.4e} {se:>10.2e}")
a = np.array(out)
w = 1.0 / (a[:, 2] / a[:, 3]) ** 2
print(f"\ninverse-variance mean of the ratio over the four gamma: "
      f"{np.average(a[:,1]/a[:,3], weights=w):.4f}")
print(f"chi2 against 1: {np.sum(((a[:,1]/a[:,3] - 1)/(a[:,2]/a[:,3]))**2):.2f} "
      f"for {len(a)} points")
np.savez("_markers_pooled.npz", gamma=a[:, 0], fk=a[:, 1], err=a[:, 2],
         analytic=a[:, 3], n_blocks=len(BLK), seeds_per_block=24)
print("-> _markers_pooled.npz")
