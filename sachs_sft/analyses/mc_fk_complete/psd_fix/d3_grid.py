"""Do the pathological lambda coincide with the corr_op table's own 21-point grid?"""
import sys, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op

d = np.load(corr_op.TABLE_PATH, allow_pickle=True)
print("table keys:", [k for k in d.files][:20])
for k in d.files:
    a = d[k]
    if a.dtype.kind in "fi" and a.ndim == 1 and 5 <= a.size <= 60:
        print(f"  {k}: size {a.size}, range [{a.min():.1f}, {a.max():.1f}]")
        if a.size <= 30:
            print(f"    {np.round(a, 1)}")
