"""Independent verification of three load-bearing claims (synthesis agent)."""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
_MC = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt")
sys.path.insert(0, str(_MC))
import driver_stats as ds
import corr_op_C_callable as co
from sigma2_repaired import R1Knots, TABLE_LAM

N1 = ds._N1
n2_1 = ds._n2_of_cos(1.0)

print("TABLE_LAM =", TABLE_LAM)

# ---- V1: is C(t,t) exactly piecewise quadratic in a cell? -------------------
print("\n[V1] deg-2 fit residual of C00(t,t) inside cells (cos=1)")
for (a,b) in [(1150.,1300.),(1300.,1447.),(1700.,1900.),(2160.,2250.)]:
    t = np.linspace(a+1e-6, b-1e-6, 25)
    C = np.array([co.C_fn(N1, float(x), n2_1, float(x))[0,0] for x in t])
    co2 = np.polyfit(t, C, 2); r2 = np.max(np.abs(np.polyval(co2,t)-C))/np.max(np.abs(C))
    co3 = np.polyfit(t, C, 3); r3 = np.max(np.abs(np.polyval(co3,t)-C))/np.max(np.abs(C))
    print(f"  cell [{a:.0f},{b:.0f}] deg2 rel resid {r2:.2e}   deg3 {r3:.2e}")

print("\n[V1b] one-sided dC00/dlam at node 1300 (h -> 0)")
node = 1300.0
for h in (10.,5.,2.,1.,0.5,0.1):
    Lm = (co.C_fn(N1,node,n2_1,node)[0,0]-co.C_fn(N1,node-h,n2_1,node-h)[0,0])/h
    Rp = (co.C_fn(N1,node+h,n2_1,node+h)[0,0]-co.C_fn(N1,node,n2_1,node)[0,0])/h
    print(f"  h={h:6.2f}  left {Lm:+.6e}   right {Rp:+.6e}")
Cm = co.C_fn(N1,node-1e-9,n2_1,node-1e-9)[0,0]; Cp = co.C_fn(N1,node+1e-9,n2_1,node+1e-9)[0,0]
print(f"  value continuity: C(1300-)={Cm:.12e}  C(1300+)={Cp:.12e}  reldiff {abs(Cp-Cm)/Cm:.2e}")

# ---- V2: PSD at the production failure node --------------------------------
print("\n[V2] 6x6 eigenvalues at the failure node, stock vs R1  (apply_c0=False)")
stock = ds.Sigma2Builder(apply_c0=False)
r1 = R1Knots(apply_c0=False)
def six(bld, cg, lam):
    w = bld.matrix(1.0, lam); c = bld.matrix(cg, lam)
    M = np.zeros((6,6)); M[:3,:3]=w; M[3:,3:]=w; M[:3,3:]=c; M[3:,:3]=c.T
    return M
for g_arcmin in (1.0, 17.3):
    cg = math.cos(math.radians(g_arcmin/60.0))
    for lam in (1308.93, 1310.84):
        es = np.linalg.eigvalsh(six(stock,cg,lam)); er = np.linalg.eigvalsh(six(r1,cg,lam))
        print(f"  g={g_arcmin:5.2f}' lam={lam:8.2f}  stock minEig/maxEig {es.min()/es.max():+.4e}"
              f"   R1 {er.min()/er.max():+.4e}")
    print(f"     stock Sigma00(1308.93)={stock.matrix(1.0,1308.93)[0,0]:+.6e}"
          f"   R1 {r1.matrix(1.0,1308.93)[0,0]:+.6e}")
