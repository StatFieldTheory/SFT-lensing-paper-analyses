"""D1. Are the streams default_rng(base + 991*s) correlated?

Direct test on the raw draws, not on the estimator: for one block of S seeds,
take the first L standard normals of each stream and form the full SxS
correlation matrix, both of x and of x^2-1 (the estimator is quartic in x, so
correlations in the squares are what would matter).  Compare against
SeedSequence-spawned streams and against a null of S independent streams.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators"
                   "/code/mc_fk_complete/inputmatch")

S, L = 192, 60_000
BASES = [20260827, 20260828, 314159001, 8675309, 271828183, 6010001, 11000003]

def offdiag(C):
    iu = np.triu_indices(C.shape[0], 1)
    return C[iu]

def report(tag, X):
    Xc = X - X.mean(axis=1, keepdims=True)
    Xc /= np.linalg.norm(Xc, axis=1, keepdims=True)
    c = offdiag(Xc @ Xc.T)
    Y = X ** 2 - 1.0
    Yc = Y - Y.mean(axis=1, keepdims=True)
    Yc /= np.linalg.norm(Yc, axis=1, keepdims=True)
    c2 = offdiag(Yc @ Yc.T)
    exp = 1.0 / np.sqrt(L)
    print(f"  {tag:>26}  x: rms={c.std():.5f} max|c|={np.abs(c).max():.5f} "
          f"({np.abs(c).max()/exp:.2f} sd)   x^2: rms={c2.std():.5f} "
          f"max|c|={np.abs(c2).max():.5f} ({np.abs(c2).max()/exp:.2f} sd)")
    return c, c2

print(f"S={S} streams, L={L} draws each; expected |corr| sd = {1/np.sqrt(L):.5f}")
for base in BASES:
    X = np.empty((S, L))
    for s in range(S):
        X[s] = np.random.default_rng(base + 991 * s).standard_normal(L)
    report(f"int base={base}", X)
ss = np.random.SeedSequence(12345).spawn(S)
X = np.array([np.random.default_rng(c).standard_normal(L) for c in ss])
report("SeedSequence spawn", X)
X = np.array([np.random.default_rng(int(1e9) + 7919 * s).standard_normal(L) for s in range(S)])
report("int step=7919", X)
# consecutive integers: the harshest case for a counter-based seeder
X = np.array([np.random.default_rng(555_000_000 + s).standard_normal(L) for s in range(S)])
report("int step=1", X)
