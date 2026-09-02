import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import _bootstrap
ds, bg, _ = _bootstrap.wire()
cg = float(np.cos(1.0*np.pi/(180*60)))
N=1000; LMIN=406.0; LF=2313.029
dl=(LF-LMIN)/(N-1); lam=LMIN+dl*np.arange(N)
S6=np.array([ds.sigma2_6x6(cg,float(l)) for l in lam])
Z6=np.array([ds.zeta6(cg,float(l)) for l in lam])
print("Sigma2_00: min %.3e max %.3e"%(S6[:,0,0].min(),S6[:,0,0].max()))
zc = -(Z6[:,0,0,3]+Z6[:,1,1,3]+Z6[:,2,2,3])
print("zeta contraction -sum_c z_cc3: min %.3e max %.3e"%(zc.min(),zc.max()))
# smoothness diagnostics: relative node-to-node variation
def rough(x):
    d2 = x[2:]-2*x[1:-1]+x[:-2]
    return np.median(np.abs(d2))/np.median(np.abs(x))
print("roughness Sigma00 %.3e   zeta_c %.3e"%(rough(S6[:,0,0]), rough(zc)))
for sig in (8.0,4.0):
    rho=np.exp(-dl/sig); Sch=dl*(1+rho)/(1-rho)
    V=S6/(2*sig)
    # exact stationary P of the AR chain (clip V to PSD first)
    Vc=np.empty_like(V)
    for k in range(N):
        d,U=np.linalg.eigh(V[k]); Vc[k]=(U*np.clip(d,0,None))@U.T
    P=np.empty_like(V); P[0]=Vc[0]
    for k in range(1,N): P[k]=rho*rho*P[k-1]+(1-rho*rho)*Vc[k]
    # Chat(k) = dl[ sum_{j<=k} rho^{k-j} P[j] + sum_{j>k} rho^{j-k} P[k] ]
    acc=np.zeros((6,6)); Chat=np.empty_like(V)
    for k in range(N):
        acc = rho*acc + P[k]
        Chat[k] = dl*(acc + P[k]*rho/(1-rho))
    r = Chat[:,0,0]/(Sch*V[:,0,0])
    ok = np.isfinite(r)
    print("sigma=%.1f  rho=%.4f S=%.4f  Chat00/(S*V00): med %.4f  mean %.4f  <r^2> %.4f  p5 %.3f p95 %.3f"
          %(sig,rho,Sch,np.median(r[ok]),np.mean(r[ok]),np.mean(r[ok]**2),np.percentile(r[ok],5),np.percentile(r[ok],95)))
    rp = P[:,0,0]/Vc[:,0,0]
    print("        P00/V00: med %.4f p5 %.3f p95 %.3f"%(np.median(rp),np.percentile(rp,5),np.percentile(rp,95)))
