"""Fresh cos: MULTI-POINT batch FIRST, then scalar.  The untested ordering."""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op
N1 = np.array([0.0, 0.0, 1.0])
def n2c(c):
    c = float(np.clip(c, -1, 1)); s = float(np.sqrt(max(0.0, 1-c*c)))
    return np.array([s, 0.0, c])
TS = np.linspace(420.0, 2320.0, 180)
PROBE = (700.5, 1234.5, 2100.7)
print(f"{'gamma':>8} {'batch-first scalar':>20} {'clean scalar':>20} {'rel diff':>11}")
for gam in (0.5, 1.0, 5.0, 17.3, 2.0, 10.0, 44.4):
    cg = float(np.cos(np.deg2rad(gam/60.0))); n2 = n2c(cg)
    corr_op.C_fn_batch(N1, TS, n2, TS)                 # batch touches this cos FIRST
    a = np.array([float(corr_op.C_fn(N1, float(t), n2, float(t))[0, 0]) for t in PROBE])
    # clean reference: reload the callable in a subprocess-free way -> use a
    # perturbed-but-equal cos is not possible; instead compare to a fresh import
    print(f"{gam:>8.2f} {a[1]:>20.12e}")
print("\nFresh interpreter, scalar only (the production path):")
