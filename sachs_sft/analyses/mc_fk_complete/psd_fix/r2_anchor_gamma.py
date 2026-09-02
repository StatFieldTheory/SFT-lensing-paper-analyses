"""ATTACK (a)+(b)+(d): is the R1 shift gamma-independent, and does it treat the
cross-ray block differently from the within-ray block?

Three integrals of the SAME density, all noise-free:
  O0_GL256  : exactly what driver_stats.order0_mc / sachs_mc_core.simulate report
              as `o0_analytic` (global 256-pt Gauss-Legendre).
  O0_disc(N): sum_j dlam Sigma2[j] Wd[j]^2 -- the EXACT discrete quantity the MC
              realises on its own uniform grid (build_grid's Wd), no MC noise.
  O0_panel  : Gauss-Legendre panelised at the 21 corr_op table nodes, 64 pts per
              panel -- converged for a density with jumps at those nodes.
Full 3x3 so the xi_+/xi_-/kappa-gamma structure can be checked, and both the
cross-ray (cos=gamma) and within-ray (cos=1) blocks.
"""
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
import sigma2_repaired as R

LAM_F = bgm.LAM_SOURCE_BASELINE
LAM_LO = 406.0
bg = bgm.Background(lam_min=120.0, lam_max=LAM_F + 5.0)

def wd_of(lam):
    D = np.asarray(bg.D(lam), float)
    dlam = float(lam[1] - lam[0])
    w = np.full(lam.size, dlam); w[0] = w[-1] = 0.5 * dlam
    return D**2 * np.cumsum((w / D**2)[::-1])[::-1], dlam

def G_cont(la):
    """continuum window G(la) = int_la^{lam_f} (D(la)/D(t))^2 dt, 64-pt GL."""
    nd, wt = leggauss(64)
    Dla = np.asarray(bg.D(la), float)
    out = np.empty_like(np.asarray(la, float))
    for i, l in enumerate(np.atleast_1d(la)):
        t = 0.5*(LAM_F-l)*nd + 0.5*(LAM_F+l); jw = 0.5*(LAM_F-l)*wt
        out[i] = float(np.sum(jw * (Dla[i]/np.asarray(bg.D(t)))**2))
    return out

def mats(b, cosg, lam):
    return np.array([b.matrix(cosg, float(l)) for l in lam])

def I_disc(b, cosg, N):
    lam = np.linspace(LAM_LO, LAM_F, N)
    Wd, dlam = wd_of(lam)
    S = mats(b, cosg, lam)
    return np.einsum("j,jab->ab", dlam * Wd**2, S)

def I_gl(b, cosg, n_gauss=256):
    nd, wt = leggauss(n_gauss)
    la = 0.5*(LAM_F-LAM_LO)*nd + 0.5*(LAM_F+LAM_LO)
    jw = 0.5*(LAM_F-LAM_LO)*wt
    S = mats(b, cosg, la); G = G_cont(la)
    return np.einsum("j,jab->ab", jw*G**2, S)

BREAKS = R._cells(LAM_LO, LAM_F)          # table nodes clipped to the source
def I_panel(b, cosg, npt=64):
    nd, wt = leggauss(npt); tot = np.zeros((3,3))
    for lo, hi in zip(BREAKS[:-1], BREAKS[1:]):
        la = 0.5*(hi-lo)*nd + 0.5*(hi+lo); jw = 0.5*(hi-lo)*wt
        S = mats(b, cosg, la); G = G_cont(la)
        tot += np.einsum("j,jab->ab", jw*G**2, S)
    return tot

d = np.load("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
            "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz",
            allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)

GAM = np.array([0.5,1.0,2.0,3.0,5.0,8.0,12.0,17.3,25.0,44.4,70.0,114.3,
                180.0,293.9,450.0,600.0])
out = {}
for name, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
    b = cls(background=bg, apply_c0=False)
    rec = {k: [] for k in ("gl","d1000","d2000","d4000","panel",
                           "gl_w","d1000_w","panel_w")}
    for gam in GAM:
        c = float(np.cos(np.deg2rad(gam/60.0)))
        rec["gl"].append(I_gl(b, c)); rec["panel"].append(I_panel(b, c))
        rec["d1000"].append(I_disc(b, c, 1000)); rec["d2000"].append(I_disc(b, c, 2000))
        rec["d4000"].append(I_disc(b, c, 4000))
    # within-ray block (cos=1) once
    rec["gl_w"] = I_gl(b, 1.0); rec["panel_w"] = I_panel(b, 1.0)
    rec["d1000_w"] = I_disc(b, 1.0, 1000)
    out[name] = {k: np.array(v) for k, v in rec.items()}
    print(f"[{name}] done", flush=True)

np.savez("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/"
         "mc_fk_complete/psd_fix/_r2.npz", GAM=GAM, g_a3=g_a3, o0_a3=o0_a3,
         **{f"{n}_{k}": v for n, r in out.items() for k, v in r.items()})

a3 = np.interp(GAM, g_a3, o0_a3)
print("\n=== kappa-kappa (cross block, [0,0]) ===")
print(f"{'gam':>7} {'O0_a3':>12} | {'stock GL256':>11} {'R1 GL256':>10} {'R1/st':>8}"
      f" | {'stock disc1k':>12} {'R1 disc1k':>11} {'R1/st':>8}"
      f" | {'stock panel':>11} {'R1 panel':>10} {'R1/st':>8}")
for i, gam in enumerate(GAM):
    s_gl, r_gl = out["stock"]["gl"][i,0,0], out["R1"]["gl"][i,0,0]
    s_d,  r_d  = out["stock"]["d1000"][i,0,0], out["R1"]["d1000"][i,0,0]
    s_p,  r_p  = out["stock"]["panel"][i,0,0], out["R1"]["panel"][i,0,0]
    print(f"{gam:>7.1f} {a3[i]:>12.4e} | {s_gl:>11.4e} {r_gl:>10.4e} {r_gl/s_gl:>8.5f}"
          f" | {s_d:>12.4e} {r_d:>11.4e} {r_d/s_d:>8.5f}"
          f" | {s_p:>11.4e} {r_p:>10.4e} {r_p/s_p:>8.5f}")
print("\n=== ratio to O0_analysis3 (the ANCHOR), kk ===")
for rule in ("gl","d1000","panel"):
    s = out["stock"][rule][:,0,0]/a3; r = out["R1"][rule][:,0,0]/a3
    print(f"{rule:>7}: stock median {np.median(s[:12]):.5f}  R1 median {np.median(r[:12]):.5f}"
          f"  shift {np.median(r[:12])/np.median(s[:12])-1:+.4%}")
