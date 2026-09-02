"""Is the +1.61% anchor shift real, or a quadrature artifact of a global
Gauss-Legendre rule applied to a DISCONTINUOUS integrand?"""
import sys, time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, _bg, _ds
import sigma2_repaired as R

bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"),
            allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
GAM = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3, 293.9, 596.9])

B = {"stock": _ds.Sigma2Builder(background=bgd, apply_c0=False),
     "R1Knots": R.R1Knots(background=bgd, apply_c0=False),
     "indep": IndepDensity(background=bgd, apply_c0=False)}

# --- quadrature abscissae + weights, window precomputed once per rule --------
def rule_global(n):
    nod, wt = leggauss(n)
    la = 0.5*(LAMF-406.0)*nod + 0.5*(LAMF+406.0)
    return la, 0.5*(LAMF-406.0)*wt

def rule_panel(npp):
    edges = [406.0] + [float(x) for x in NODES if 406.0 < x < LAMF] + [LAMF]
    nod, wt = leggauss(npp); LA = []; W = []
    for a, c in zip(edges[:-1], edges[1:]):
        LA.append(0.5*(c-a)*nod + 0.5*(c+a)); W.append(0.5*(c-a)*wt)
    return np.concatenate(LA), np.concatenate(W)

def rule_trap(n):
    la = np.linspace(406.0, LAMF, n)
    w = np.full(n, (LAMF-406.0)/(n-1)); w[0] *= 0.5; w[-1] *= 0.5
    return la, w

RULES = {}
for name, (la, w) in (("GL256", rule_global(256)), ("GL1024", rule_global(1024)),
                      ("GL4096", rule_global(4096)), ("panel48", rule_panel(48)),
                      ("trap2001", rule_trap(2001))):
    RULES[name] = (la, w, _ds.order0_window(bgd, la, LAMF, n_gauss=64)**2)
print("abscissae per rule:", {k: len(v[0]) for k, v in RULES.items()})

med = {}
for bn, b in B.items():
    t0 = time.time(); row = {}
    for rn, (la, w, G2) in RULES.items():
        rr = []
        for gam in GAM:
            cg = float(np.cos(np.deg2rad(gam/60.0)))
            s = np.array([b.matrix(cg, float(l))[0, 0] for l in la])
            rr.append(float(np.sum(w*s*G2))/float(np.interp(gam, g_a3, o0_a3)))
        row[rn] = float(np.median(rr))
    med[bn] = row
    print(f"  {bn} done in {time.time()-t0:.1f}s")

print("\nmedian over the 10 paper gammas of  order0 / O0_analysis3 (raw)")
print(f"{'builder':>10} " + " ".join(f"{r:>10}" for r in RULES))
for bn in B:
    print(f"{bn:>10} " + " ".join(f"{med[bn][r]:>10.5f}" for r in RULES))
print("\nshift vs stock, SAME rule:")
print(f"{'builder':>10} " + " ".join(f"{r:>10}" for r in RULES))
for bn in B:
    print(f"{bn:>10} " + " ".join(f"{med[bn][r]/med['stock'][r]-1:>+9.4%} " for r in RULES))
print("\nrule dependence WITHIN each builder (relative to panel48):")
for bn in B:
    print(f"  {bn:>10}: " + " ".join(
        f"{r} {med[bn][r]/med[bn]['panel48']-1:+.4%}" for r in RULES))
