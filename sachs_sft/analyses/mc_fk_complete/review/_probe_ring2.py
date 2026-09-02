import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import _bootstrap
ds, bg, _ = _bootstrap.wire()
cg = float(np.cos(1.0*np.pi/(180*60)))
N=1000; LMIN=406.0; LF=2313.029
dl=(LF-LMIN)/(N-1); lam=LMIN+dl*np.arange(N)
B = bg.Background(lam_min=120.0, lam_max=LF+5.0)
D=np.asarray(B.D(lam))
w=np.full(N,dl); w[0]=w[-1]=dl/2
# discrete window Wd[j] = sum_{k>=j} w_k (D[j]/D[k])^2
tail=np.cumsum((w/D**2)[::-1])[::-1]; Wd=D**2*tail
# H[j] = D^4 sum_{k>=j} w_k Wd[k]/D[k]^4
t4=np.cumsum((w*Wd/D**4)[::-1])[::-1]; Hd=D**4*t4
S6=np.array([ds.sigma2_6x6(cg,float(l)) for l in lam])
Z6=np.array([ds.zeta6(cg,float(l)) for l in lam])
zc = -(Z6[:,0,0,3]+Z6[:,1,1,3]+Z6[:,2,2,3])
integ = Wd*Hd*zc
print("FK continuum-on-MC-grid = %.6e"%(2*np.sum(w*integ)))
wt = w*integ; wt=wt/wt.sum()
cum=np.cumsum(wt)
print("lambda at 10/25/50/75/90 pct of FK weight:",
      [float(np.interp(p,cum,lam)) for p in (0.1,0.25,0.5,0.75,0.9)])
for sig in (8.0,4.0):
    rho=np.exp(-dl/sig); Sch=dl*(1+rho)/(1-rho)
    V=S6/(2*sig)
    Vc=np.empty_like(V)
    for k in range(N):
        d,U=np.linalg.eigh(V[k]); Vc[k]=(U*np.clip(d,0,None))@U.T
    P=np.empty_like(V); P[0]=Vc[0]
    for k in range(1,N): P[k]=rho*rho*P[k-1]+(1-rho*rho)*Vc[k]
    acc=np.zeros((6,6)); Chat=np.empty_like(V)
    for k in range(N):
        acc = rho*acc + P[k]
        Chat[k] = dl*(acc + P[k]*rho/(1-rho))
    r = Chat[:,0,0]/(Sch*V[:,0,0])
    print("sig=%.1f kernel-weighted <r>=%.4f <r^2>=%.4f   plain med %.4f"
          %(sig, float(np.sum(wt*r)), float(np.sum(wt*r*r)), np.median(r)))
    # where is r far from 1, and how much weight there
    bad = np.abs(r-1)>0.5
    print("        weight in |r-1|>0.5 nodes: %.4f (n=%d)"%(wt[bad].sum(), bad.sum()))
# profile of V00 vs lam in the FK-weighted core
print("\nlam, V00 (sig8), Wd*Hd*zc normalised")
for k in range(0,N,50):
    print("%8.1f  %.4e  %.4e"%(lam[k], S6[k,0,0], wt[k]))
