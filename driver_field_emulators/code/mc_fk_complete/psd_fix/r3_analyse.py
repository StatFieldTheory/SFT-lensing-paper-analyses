import numpy as np
P="/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix/_r2.npz"
d=np.load(P); GAM=d["GAM"]; a3=np.interp(GAM,d["g_a3"],d["o0_a3"])
S,R = {k[6:]:d[k] for k in d if k.startswith("stock_")}, {k[3:]:d[k] for k in d if k.startswith("R1_")}

print("=== MC-realised discrete O0 vs what order0_mc(n_gauss=256) REPORTS (kk) ===")
print(f"{'gam':>7} | {'stock d1k/GL':>12} {'stock d2k/GL':>12} {'stock d4k/GL':>12}"
      f" | {'R1 d1k/GL':>10} {'R1 d2k/GL':>10} {'R1 d4k/GL':>10}")
for i,g in enumerate(GAM[:13]):
    print(f"{g:>7.1f} | " + " ".join(f"{S[k][i,0,0]/S['gl'][i,0,0]:>12.5f}" for k in ("d1000","d2000","d4000"))
          + " | " + " ".join(f"{R[k][i,0,0]/R['gl'][i,0,0]:>10.5f}" for k in ("d1000","d2000","d4000")))

print("\n=== convergence of the discrete rule to the panelised truth (kk) ===")
print(f"{'gam':>7} | {'stock d1k/pan':>13} {'d2k/pan':>9} {'d4k/pan':>9} | {'R1 d1k/pan':>11} {'d2k/pan':>9} {'d4k/pan':>9}")
for i,g in enumerate(GAM[:13]):
    print(f"{g:>7.1f} | " + " ".join(f"{S[k][i,0,0]/S['panel'][i,0,0]:>13.5f}" if k=='d1000' else f"{S[k][i,0,0]/S['panel'][i,0,0]:>9.5f}" for k in ("d1000","d2000","d4000"))
          + " | " + " ".join(f"{R[k][i,0,0]/R['panel'][i,0,0]:>11.5f}" if k=='d1000' else f"{R[k][i,0,0]/R['panel'][i,0,0]:>9.5f}" for k in ("d1000","d2000","d4000")))

print("\n=== ATTACK (a): is the GL-256 shift gamma-INDEPENDENT?  R1/stock, kk ===")
for rule in ("gl","d4000","panel"):
    r = R[rule][:,0,0]/S[rule][:,0,0]
    lo,hi = r[:12].min(), r[:12].max()
    print(f"{rule:>7}: 0.5'-114' min {lo:.5f} max {hi:.5f}  spread {(hi/lo-1)*100:+.3f}%   "
          f"full 0.5'-600' min {r.min():.5f} max {r.max():.5f} spread {(r.max()/r.min()-1)*100:+.3f}%")

print("\n=== ATTACK (b)+(d): per-component R1/stock, CROSS block (panel rule) ===")
lbl = ["kk","k g+","k gx","g+ k","g+g+","g+gx","gx k","gx g+","gxgx"]
print(f"{'gam':>7} " + " ".join(f"{l:>8}" for l in lbl))
for i,g in enumerate(GAM[:13]):
    rr = (R["panel"][i]/S["panel"][i]).ravel()
    print(f"{g:>7.1f} " + " ".join(f"{v:>8.5f}" for v in rr))

print("\n=== same, GL-256 rule (what the code actually uses) ===")
for i,g in enumerate(GAM[:13]):
    rr = (R["gl"][i]/S["gl"][i]).ravel()
    print(f"{g:>7.1f} " + " ".join(f"{v:>8.5f}" for v in rr))

print("\n=== WITHIN-ray block (cos=1) R1/stock ===")
for rule in ("gl_w","d1000_w","panel_w"):
    print(f"{rule:>9}: " + " ".join(f"{v:>9.5f}" for v in (R[rule]/S[rule]).ravel()))

print("\n=== xi channels built from the cross block (panel rule) ===")
print(f"{'gam':>7} {'xi+ R1/st':>10} {'xi- R1/st':>10} {'kg R1/st':>10} {'kk R1/st':>10}"
      f" | {'xi-/xi+ stock':>13} {'xi-/xi+ R1':>11}")
for i,g in enumerate(GAM[:13]):
    for nm,M in (("s",S["panel"][i]),("r",R["panel"][i])):
        xp = M[1,1]+M[2,2]; xm = M[1,1]-M[2,2]; kg = M[0,1]
        if nm=="s": sp,sm,sk = xp,xm,kg
        else: rp,rm,rk = xp,xm,kg
    print(f"{g:>7.1f} {rp/sp:>10.5f} {rm/sm:>10.5f} {rk/sk:>10.5f} "
          f"{R['panel'][i,0,0]/S['panel'][i,0,0]:>10.5f} | {sm/sp:>13.5e} {rm/rp:>11.5e}")
