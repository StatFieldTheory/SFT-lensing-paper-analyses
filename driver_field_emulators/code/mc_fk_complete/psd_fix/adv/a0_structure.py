"""A0: is the corr_op equal-time diagonal EXACTLY piecewise quadratic in lambda?

If yes, the 'true' density is analytic, not an FD convention, and the attack's
sub-claim ("the true density defined by finite differences is itself arbitrary")
fails on structure alone.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as CO

LAM = np.load(CO.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)
print("table nodes:", LAM)
N1 = ds._N1

def Cdiag(cosg, t):
    n2 = ds._n2_of_cos(cosg)
    return np.asarray(CO.C_fn(N1, float(t), n2, float(t)), float)

for cosg, tag in ((1.0, "cos=1"), (float(np.cos(np.deg2rad(1.0/60))), "gam=1'"),
                  (float(np.cos(np.deg2rad(17.3/60))), "gam=17.3'")):
    worst2 = worst3 = 0.0
    for i in range(len(LAM)-1):
        lo, hi = LAM[i], LAM[i+1]
        xs = np.linspace(lo+1e-6*(hi-lo), hi-1e-6*(hi-lo), 9)
        F = np.array([Cdiag(cosg, x) for x in xs])
        u = (xs-lo)/(hi-lo)
        for a in range(3):
            for b in range(3):
                y = F[:,a,b]
                s = max(np.max(np.abs(y)), 1e-300)
                for deg, store in ((2,'2'),(3,'3')):
                    c = np.polyfit(u, y, deg)
                    r = np.max(np.abs(np.polyval(c,u)-y))/s
                    if deg==2: worst2 = max(worst2, r)
                    else: worst3 = max(worst3, r)
    print(f"{tag}: max rel residual  deg2={worst2:.2e}   deg3={worst3:.2e}")
