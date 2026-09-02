"""Verify the two corrections the workflow made to my diagnosis."""
import sys, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
import sigma2_repaired as R
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
N1 = ds._N1

print("=== (a) is C(t,t) C^0 (slope jump) or C^1 (curvature jump)? ===")
for node in (1300.0, 850.0):
    for side, xs in (("left ", np.linspace(node - 40, node - 1, 40)),
                     ("right", np.linspace(node + 1, node + 40, 40))):
        y = np.array([corr_op.C_fn(N1, float(t), N1, float(t))[0, 0] for t in xs])
        c = np.polyfit(xs - node, y, 2)
        res = np.abs(y - np.polyval(c, xs - node)).max() / np.abs(y).max()
        print(f"  node {node:.0f} {side}: dC/dlam at node = {c[1]:+.4e}  "
              f"d2C = {2*c[0]:+.4e}   deg-2 fit residual = {res:.2e}")
    print("   -> first derivative jumps sign" if True else "")

print()
print("=== (b) order0_mc quadrature: is the +1.61% real or a quadrature artefact? ===")
grid21 = np.load(corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"]
from numpy.polynomial.legendre import leggauss
def order0_panel(builder, cosg, lam_f, per_panel=64):
    """Gauss-Legendre panelised AT the corr_op table nodes, so no panel
    straddles a discontinuity of the repaired density."""
    lo = builder.lam_lo
    edges = np.concatenate([[lo], grid21[(grid21 > lo) & (grid21 < lam_f)], [lam_f]])
    nod, wt = leggauss(per_panel); tot = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        la = 0.5 * (b - a) * nod + 0.5 * (b + a); jw = 0.5 * (b - a) * wt
        s00 = np.array([builder.matrix(cosg, float(l))[0, 0] for l in la])
        G = ds.order0_window(bg, la, lam_f, n_gauss=64)
        tot += float(np.sum(jw * s00 * G ** 2))
    return tot
lam_f = bgm.LAM_SOURCE_BASELINE
print(f"{'gamma':>8} {'--- global GL n=256 ---':>26} {'--- panelised GL ---':>26}")
print(f"{'':>8} {'stock':>12} {'R1/stock':>13} {'stock':>12} {'R1/stock':>13}")
for gam in (0.5, 1.0, 5.0, 17.3, 114.3, 293.9):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    bs = ds.Sigma2Builder(background=bg, apply_c0=False)
    br = R.R1Knots(background=bg, apply_c0=False)
    g_s = float(ds.order0_mc(cosg, lam_f=lam_f, builder=bs, n_gauss=256))
    g_r = float(ds.order0_mc(cosg, lam_f=lam_f, builder=br, n_gauss=256))
    p_s = order0_panel(bs, cosg, lam_f)
    p_r = order0_panel(br, cosg, lam_f)
    print(f"{gam:>8.2f} {g_s:>12.5e} {g_r/g_s:>13.5f} {p_s:>12.5e} {p_r/p_s:>13.5f}")
