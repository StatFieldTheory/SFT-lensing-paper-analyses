import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import _bootstrap
ds, bg, _ = _bootstrap.wire()
cg=float(np.cos(1.0*np.pi/(180*60)))
LMIN,LF=406.0,2313.029
B=bg.Background(lam_min=120.0,lam_max=LF+5.0)
def kern_cont(n):
    u=np.linspace(LMIN,LF,n); Du=np.asarray(B.D(u)); du=np.diff(u)
    i2=1/Du**2; c2=np.concatenate(([0],np.cumsum(.5*(i2[1:]+i2[:-1])*du)))
    W=Du**2*(c2[-1]-c2)
    i4=W/Du**4; c4=np.concatenate(([0],np.cumsum(.5*(i4[1:]+i4[:-1])*du)))
    H=Du**4*(c4[-1]-c4)
    return u,W,H
for n in (1001,4001,20001,80001):
    u,W,H=kern_cont(n)
    print("n=%6d  W(1000)=%.6e H(1000)=%.6e  W(2000)=%.6e H(2000)=%.6e"%(
        n,np.interp(1000,u,W),np.interp(1000,u,H),np.interp(2000,u,W),np.interp(2000,u,H)))
# discrete
for N in (1000,2000,4000):
    dl=(LF-LMIN)/(N-1); lam=LMIN+dl*np.arange(N); D=np.asarray(B.D(lam))
    w=np.full(N,dl); w[0]=w[-1]=dl/2
    tail=np.cumsum((w/D**2)[::-1])[::-1]; Wd=D**2*tail
    inv=np.zeros(N); inv[:-1]=dl*Wd[1:]/D[:-1]**4
    Hd=D**4*np.cumsum(inv[::-1])[::-1]
    print("N=%5d  Wd(1000)=%.6e Hd(1000)=%.6e  Wd(2000)=%.6e Hd(2000)=%.6e"%(
        N,np.interp(1000,lam,Wd),np.interp(1000,lam,Hd),np.interp(2000,lam,Wd),np.interp(2000,lam,Hd)))
# zeta smoothness on the fine grid
lz=np.linspace(1900,2100,201)
zz=np.array([-(lambda Z: Z[0,0,3]+Z[1,1,3]+Z[2,2,3])(ds.zeta6(cg,float(x))) for x in lz])
d2=zz[2:]-2*zz[1:-1]+zz[:-2]
print("zeta rel |2nd diff| median: %.3e"%np.median(np.abs(d2)/np.abs(zz[1:-1])))
# full discrete FK vs N
for N in (1000,2000,4000,8000):
    dl=(LF-LMIN)/(N-1); lam=LMIN+dl*np.arange(N); D=np.asarray(B.D(lam))
    w=np.full(N,dl); w[0]=w[-1]=dl/2
    tail=np.cumsum((w/D**2)[::-1])[::-1]; Wd=D**2*tail
    inv=np.zeros(N); inv[:-1]=dl*Wd[1:]/D[:-1]**4
    Hd=D**4*np.cumsum(inv[::-1])[::-1]
    zc=np.array([-(lambda Z: Z[0,0,3]+Z[1,1,3]+Z[2,2,3])(ds.zeta6(cg,float(x))) for x in lam])
    print("N=%5d  FK_discrete=%.6e   (trap) %.6e"%(N,2*np.sum(dl*Wd*Hd*zc),2*np.sum(w*Wd*Hd*zc)))
