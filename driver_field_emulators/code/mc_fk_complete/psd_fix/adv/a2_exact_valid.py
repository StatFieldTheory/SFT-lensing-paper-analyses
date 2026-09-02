"""A2: validate the parameter-free Exact density three ways, then reproduce the
brief's d6 numbers so the harness is anchored to the published protocol."""
import numpy as np, sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix/adv")
import advlib as A
ds = A.ds

bg = A.BG
E  = A.Exact(background=bg, apply_c0=False)
R0 = A.R.R0Reference(background=bg, apply_c0=False)

cos1 = float(np.cos(np.deg2rad(1.0/60)))
# (a) Exact vs R0 with a range of FD steps, MID-CELL (far from any node)
mids = [(A.TABLE_LAM[i]+A.TABLE_LAM[i+1])/2 for i in range(11)]
print("(a) Exact vs R0 mid-cell, Sigma00, per fd_h:")
for h in (4.0, 2.0, 1.0, 0.5, 0.2, 0.05, 0.01):
    r0 = A.R.R0Reference(background=bg, apply_c0=False, fd_h=h)
    e = np.array([E.matrix(cos1, m)[0,0] for m in mids])
    v = np.array([r0.matrix(cos1, m)[0,0] for m in mids])
    print(f"   fd_h={h:6.3f}  max|R0/Exact-1| = {np.max(np.abs(v/e-1)):.3e}")

# (b) integral identity: int_410^t D^4 Sigma dlam + D^4 C|_410  ==  D^4 C(t,t)
from numpy.polynomial.legendre import leggauss
def d4C(t):
    C,_ = E.C_and_dC(cos1, t); return float(bg.D(t))**4 * C[0,0]
def integ(t):
    inner = A.TABLE_LAM[(A.TABLE_LAM>406+1e-9)&(A.TABLE_LAM<t-1e-9)]
    ed = np.concatenate([[406.0], inner, [t]]); x,w = leggauss(40); s=0.0
    for lo,hi in zip(ed[:-1],ed[1:]):
        la = 0.5*(hi-lo)*x+0.5*(hi+lo); ww = 0.5*(hi-lo)*w
        s += float(np.sum(ww*np.array([float(bg.D(l))**4*E.matrix(cos1,float(l))[0,0] for l in la])))
    return s
print("(b) integral identity  int D^4 Sigma + D^4C(406) vs D^4C(t):")
for t in (600.,1000.,1300.,1301.,1800.,2200.,2310.):
    lhs = integ(t)+d4C(406.0+1e-9); rhs = d4C(t)
    print(f"   t={t:7.1f}  rel err = {abs(lhs/rhs-1):.3e}")

# (c) reproduce the brief / d6 numbers exactly (GL-256 global rule)
print("(c) d6 protocol (global GL-256), median order0_mc/O0_analysis3:")
b_cur = ds.Sigma2Builder(background=bg, apply_c0=False)
b_r1  = A.R.R1Knots(background=bg, apply_c0=False)
b_r0  = A.R.R0Reference(background=bg, apply_c0=False)
b_r2  = A.R.R2Density(background=bg, apply_c0=False)
for name,b in (("current",b_cur),("R1 knots",b_r1),("R0 reference",b_r0),
               ("R2 density",b_r2),("EXACT",E)):
    m = float(np.median(A.anchor(b, lambda lo,hi: A.rule_gl(lo,hi,256))))
    print(f"   {name:14s} {m:.5f}")
