"""What the STOCK MC actually realises at the indefinite nodes.

At the defect nodes _cholesky_psd replaces M by Lf Lf^T (eigenvalues clipped at
0).  The MC therefore injects  (Lf Lf^T)  there, not M.  The exact discrete
Order-0 the MC realises is  sum_j dlam (Lf Lf^T)[j]_00 Wd[j]^2 / dlam  -- i.e.
the floored covariance, not the builder's Sigma2.  This quantifies the bias the
repair removes, and corrects my own disc-integral comparison.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R

LAM_F = bgm.LAM_SOURCE_BASELINE; LAM_LO = 406.0
bg = bgm.Background(lam_min=120.0, lam_max=LAM_F + 5.0)
for N in (600, 1000, 2000):
    lam = np.linspace(LAM_LO, LAM_F, N); dlam = float(lam[1]-lam[0])
    D = np.asarray(bg.D(lam), float)
    w = np.full(N, dlam); w[0] = w[-1] = 0.5*dlam
    Wd = D**2*np.cumsum((w/D**2)[::-1])[::-1]
    for gam in (1.0, 17.3):
        c = float(np.cos(np.deg2rad(gam/60.0)))
        print(f"\n--- N={N} dlam={dlam:.2f} gamma={gam}' ---")
        for bname, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
            b = cls(background=bg, apply_c0=False)
            raw = np.empty(N); flo = np.empty(N); bad = []
            for k, l in enumerate(lam):
                M = core._assemble_6x6(b, c, float(l))
                Mi = M*dlam
                raw[k] = M[0,0]
                try:
                    np.linalg.cholesky(Mi); flo[k] = M[0,0]
                except np.linalg.LinAlgError:
                    ev, V = np.linalg.eigh(0.5*(Mi+Mi.T))
                    L = V*np.sqrt(np.maximum(ev,0.0))[None,:]
                    flo[k] = float((L@L.T)[0,0])/dlam
                    bad.append((k, float(l), ev.min(), ev.max(), M[0,0], flo[k]))
            I_raw = float(np.sum(dlam*raw*Wd**2)); I_flo = float(np.sum(dlam*flo*Wd**2))
            print(f"  {bname:>5}: n_floored={len(bad):>2}  "
                  f"I(builder Sigma2)={I_raw:.6e}  I(what MC injects)={I_flo:.6e}  "
                  f"MC/builder={I_flo/I_raw:.6f}")
            for k, l, e0, e1, s00, f00 in bad:
                print(f"          node {k:>4} lam={l:8.2f}  eig[min,max]=[{e0:+.3e},{e1:+.3e}]"
                      f"  Sigma2_00={s00:+.4e} -> injected {f00:+.4e}")
