"""MINIMAL REPRODUCER of the corr_op batch/scalar first-touch defect.

`corr_op.C_fn_batch` groups its samples by a ROUNDED (cos gamma, psi1, psi2)
key and caches the resulting angular-channel table.  Whichever path touches a
given cos FIRST therefore fixes what BOTH paths return afterwards.  The
callable's own self-check runs the scalar loop first, which is exactly the
ordering that hides this.

Run this file with no argument (scalar-only) and with 'batch' (batch first).
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op

BATCH_FIRST = len(sys.argv) > 1 and sys.argv[1] == "batch"
N1 = np.array([0.0, 0.0, 1.0])
TS = np.linspace(420.0, 2320.0, 180)
print("batch-first" if BATCH_FIRST else "scalar-only")
for gam in (0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4):
    c = float(np.cos(np.deg2rad(gam/60.0)))
    n2 = np.array([float(np.sqrt(max(0.0, 1-c*c))), 0.0, c])
    if BATCH_FIRST:
        corr_op.C_fn_batch(N1, TS, n2, TS)
    print(f"  gamma={gam:6.2f}'  C00(1234.5) = "
          f"{float(corr_op.C_fn(N1, 1234.5, n2, 1234.5)[0,0]):.12e}")
