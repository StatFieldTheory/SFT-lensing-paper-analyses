"""What does the Sigma2 defect actually touch?

Two checkable statements:
 (1) The FK estimator is INVARIANT under a global rescaling of Sigma2, because
     Q ~ zeta/Sigma2^2 while each of the two propagator legs is ~ Sigma2.
     If true, the 1.61% normalisation shift cannot move the FK Monte-Carlo.
 (2) The FF channel is ~ Sigma2^2, so it scales.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import fk_expect_exact as ex
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R

gam, sig, N = 1.0, 8.0, 400
cosg = float(np.cos(np.deg2rad(gam / 60.0)))

print("=== (1) is the FK estimator invariant under Sigma2 -> c * Sigma2 ? ===")
class Scaled(ds.Sigma2Builder):
    scale = 1.0
    def matrix(self, cos_gamma, lam):
        return self.scale * super().matrix(cos_gamma, lam)

base = None
for c in (1.0, 1.0161, 1.10, 2.00):
    g = ex.build_grid(n_lambda=N)
    b = Scaled(background=g.bg, apply_c0=False); b.scale = c
    g.builder = b
    V, A, Q, rho = ex.node_stats(g, cosg, sig, calibrate="smeared")
    T1, T2 = ex.expectation(g, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
    o0 = float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE, builder=b, n_gauss=256))
    if base is None: base, o0b = tot, o0
    print(f"  Sigma2 x {c:6.4f}:  FK = {tot:.6e}  (ratio {tot/base:.6f})   "
          f"Order-0 ratio {o0/o0b:.6f}")

print()
print("=== (2) the actual repair: FK exact expectation, stock vs R1 ===")
for name, cls in (("stock", ds.Sigma2Builder), ("R1 knots", R.R1Knots)):
    g = ex.build_grid(n_lambda=N)
    g.builder = cls(background=g.bg, apply_c0=False)
    V, A, Q, rho = ex.node_stats(g, cosg, sig, calibrate="smeared")
    zi = ex.zeta_injected(ex.calib_matrix(g, V, A, rho, sig, "smeared"), Q)
    ref = ex.fk_reference(g, zi)[0, 3]
    T1, T2 = ex.expectation(g, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
    qmax = np.abs(Q).reshape(N, -1).max(axis=1)
    print(f"  {name:>9}: FK exact = {tot:.6e}   local ref = {ref:.6e}   "
          f"FK/ref = {tot/ref:.5f}   max|Q| = {qmax.max():.3e}  "
          f"nodes >10x median = {(qmax > 10*np.median(qmax)).sum()}")

print()
print("=== (3) FF markers vs the analytic FF line in the deployed figure cache ===")
d = np.load("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
            "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz",
            allow_pickle=True)
g_a, ff = np.asarray(d["g"], float), np.asarray(d["ff_mom"], float)
gm, mc, se = (np.asarray(d[k], float) for k in ("g_mc", "ffc_mc", "ffc_se"))
print(f"{'gamma':>8} {'MC FF':>13} {'+/-':>10} {'analytic FF':>13} {'MC/analytic':>12}")
for i in range(min(6, len(gm))):
    a = float(np.interp(gm[i], g_a, ff))
    print(f"{gm[i]:>8.2f} {mc[i]:>13.4e} {se[i]:>10.1e} {a:>13.4e} {mc[i]/a:>12.4f}")
