"""Definitive gamma = 1' number, pooling every independent seed block."""
import numpy as np
BLOCKS = [("_gate2_highstat.npz", "base 20260827"),
          ("_gate2_highstat_B.npz", "base 314159001"),
          ("_papergrid.npz", "base 20260828")]
pool = {8.0: [], 4.0: []}
print(f"{'block':>16} {'sigma':>6} {'n':>4} {'MC':>12} {'MC/exact':>9} {'MC/ana':>8}")
for f, lbl in BLOCKS:
    d = np.load(f, allow_pickle=True)
    for sig in (8.0, 4.0):
        m = (d["gamma"] == 1.0) & (d["sigma"] == sig)
        if not m.any():
            continue
        i = int(np.where(m)[0][0])
        s = np.asarray(d["seeds"][i], float)
        pool[sig].append(s)
        print(f"{lbl:>16} {sig:>6.1f} {s.size:>4d} {s.mean():>12.5e} "
              f"{s.mean()/d['exact'][i]:>9.4f} {s.mean()/d['analytic'][i]:>8.4f}")
        ex, ana = float(d["exact"][i]), float(d["analytic"][i])
        if sig == 8.0: ex8, ana8 = ex, ana
        else: ex4 = ex
print()
for sig in (8.0, 4.0):
    s = np.concatenate(pool[sig]); se = s.std(ddof=1) / np.sqrt(s.size)
    ex = ex8 if sig == 8.0 else ex4
    print(f"POOLED sigma={sig:>4.1f}  n={s.size}  MC={s.mean():.5e} +/- {se:.2e}"
          f"   MC/exact={s.mean()/ex:.4f} +/- {se/ex:.4f}"
          f"   MC/analytic={s.mean()/ana8:.4f} +/- {se/ana8:.4f}")
n = min(len(np.concatenate(pool[8.0])), len(np.concatenate(pool[4.0])))
a8 = np.concatenate(pool[8.0])[:n]; a4 = np.concatenate(pool[4.0])[:n]
ext = (2.0 * a4 - a8) / ana8
print(f"\nPOOLED sigma_lambda -> 0 (paired, {n} seeds):"
      f"  MC / analytic FK = {ext.mean():.4f} +/- {ext.std(ddof=1)/np.sqrt(n):.4f}")
print(f"  exact expectation, same extrapolation = {(2*ex4-ex8)/ana8:.4f}")
