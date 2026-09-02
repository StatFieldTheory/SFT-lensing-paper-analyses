import sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cleanroom_mc as cm
m = cm.Model(1.0, 1000, 8.0, q_basis="chat", c_sub="P")
print("fk_discrete = %.6e ; fk_continuum = %.6e ; target %.4e"%(m.fk_discrete(), m.fk_continuum(), 1.598e-5))
rng = np.random.default_rng(7)
t0=time.time(); d = cm.run_batch(m, 4000, rng); print("4000 reals in %.1f s"%(time.time()-t0))
n=4000
def loo(x): return (x.sum()-x)/(n-1)
terms = {"X0*Kz3": (d["X0"]-loo(d["X0"]))*d["Kz3"],
         "Y0*Kg3": (d["Y0"]-loo(d["Y0"]))*d["Kg3"],
         "X3*Kz0": (d["X3"]-loo(d["X3"]))*d["Kz0"],
         "Y3*Kg0": (d["Y3"]-loo(d["Y3"]))*d["Kg0"]}
tot=sum(terms.values())
for k,v in terms.items():
    print("%-8s mean %+.4e  sd %.4e  sem %.3e"%(k, v.mean(), v.std(), v.std()/np.sqrt(n)))
print("TOTAL    mean %+.4e  sd %.4e  sem %.3e  sd/mean=%.2f"%(tot.mean(), tot.std(), tot.std()/np.sqrt(n), tot.std()/abs(tot.mean())))
print("uncentred sd/mean:", end=" ")
tot2 = d["X0"]*d["Kz3"]+d["Y0"]*d["Kg3"]+d["X3"]*d["Kz0"]+d["Y3"]*d["Kg0"]
print("%.2f"%(tot2.std()/abs(tot2.mean())))
# correlations between the ray0 and ray3 halves
h0 = terms["X0*Kz3"]+terms["Y0*Kg3"]; h3=terms["X3*Kz0"]+terms["Y3*Kg0"]
print("corr(h0,h3) = %.3f ; sd h0 %.3e sd h3 %.3e"%(np.corrcoef(h0,h3)[0,1], h0.std(), h3.std()))
print("means: X0 %.3e Y0 %.3e Kz3 %.3e Kg3 %.3e"%(d["X0"].mean(),d["Y0"].mean(),d["Kz3"].mean(),d["Kg3"].mean()))
print("sds  : X0 %.3e Y0 %.3e Kz3 %.3e Kg3 %.3e"%(d["X0"].std(),d["Y0"].std(),d["Kz3"].std(),d["Kg3"].std()))
print("corr(Y0,Kg3)=%.4f corr(X0,Kz3)=%.4f"%(np.corrcoef(d["Y0"],d["Kg3"])[0,1], np.corrcoef(d["X0"],d["Kz3"])[0,1]))
# kurtosis check
print("kurtosis of total: %.1f"%(((tot-tot.mean())**4).mean()/tot.var()**2))
