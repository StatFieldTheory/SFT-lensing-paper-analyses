import numpy as np
d=np.load("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix/_r8.npz")
G=d["GAM"]; S=d["stock_I"]; R=d["R1_I"]
def col(I,nm):
    return {"kk":I[:,0,0],"xi+":I[:,1,1]+I[:,2,2],"xi-":I[:,1,1]-I[:,2,2],"kg":I[:,0,1]}[nm]
print(f"{'gam':>8} " + " ".join(f"{n+' stock':>12} {n+' R1/st':>9}" for n in ("kk","xi+")))
for i in range(len(G)):
    if G[i] < 100: continue
    print(f"{G[i]:>8.2f} " + " ".join(
        f"{col(S,n)[i]:>12.4e} {col(R,n)[i]/col(S,n)[i]:>9.6f}" for n in ("kk","xi+")))
print("\n=== R1/stock restricted to |value| > 1% of its peak (signal region) ===")
for n in ("kk","xi+","xi-","kg"):
    s=col(S,n); r=col(R,n); m=np.abs(s) > 0.01*np.abs(s).max()
    q=r[m]/s[m]
    print(f"  {n:>3}: n={m.sum():>2}/60  gam in [{G[m].min():.2f}',{G[m].max():.2f}']"
          f"  min {q.min():.6f} max {q.max():.6f} spread {(q.max()/q.min()-1)*100:+.4f}%")
print("\n=== and restricted to gamma <= 180' (the regime the paper uses) ===")
for n in ("kk","xi+","xi-","kg"):
    s=col(S,n); r=col(R,n); m=(G<=180.0)
    q=r[m]/s[m]
    print(f"  {n:>3}: min {q.min():.6f} max {q.max():.6f} spread {(q.max()/q.min()-1)*100:+.4f}%")
