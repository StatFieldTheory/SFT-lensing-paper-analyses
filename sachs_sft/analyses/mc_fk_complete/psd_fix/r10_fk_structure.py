"""Independent check of the structural reason FK cannot move with the builder.

`fk_reference` is evaluated against `zeta_injected(M, Q)` -- the cumulant the
deformation actually injects.  If  cum3_from_Q(solve_Q(M, Z), M) == sym(Z)
exactly, then the injected cumulant is the SAME tensor whatever covariance M the
deformation was calibrated against, so the local FK reference is blind to the
builder.  Tested here on the real Sigma2/zeta at several nodes with BOTH
builders, and on random M/Z as a negative-free control.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
import sigma2_repaired as R
import fk_expect_exact as ex

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
c = float(np.cos(np.deg2rad(1.0/60.0)))
LAMS = [800.0, 1308.93, 1310.84, 1400.0, 2100.65, 2280.0]
print("=== zeta_injected(M, solve_Q(M, Z)) vs sym(Z), per node, both builders ===")
print(f"{'lam':>10} {'builder':>7} {'Sigma2_00':>13} {'max|inj-Z|/|Z|':>15} "
      f"{'Q scale':>11}")
store = {}
for l in LAMS:
    Z = ds.zeta6(c, float(l))
    for bname, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
        b = cls(background=bg, apply_c0=False)
        M = core._assemble_6x6(b, c, float(l))
        Q = ds.solve_Q(M, Z)
        inj = ds.cum3_from_Q(Q, M)
        Zs = (Z + Z.transpose(0,2,1) + Z.transpose(1,0,2) + Z.transpose(1,2,0)
              + Z.transpose(2,0,1) + Z.transpose(2,1,0)) / 6.0
        err = np.max(np.abs(inj - Zs)) / max(np.max(np.abs(Zs)), 1e-300)
        store[(l, bname)] = inj
        print(f"{l:>10.2f} {bname:>7} {M[0,0]:>13.4e} {err:>15.3e} "
              f"{np.max(np.abs(Q)):>11.3e}")
    d = np.max(np.abs(store[(l,'stock')] - store[(l,'R1')])) / \
        max(np.max(np.abs(store[(l,'stock')])), 1e-300)
    print(f"{'':>10} {'-> injected cumulant stock vs R1, max rel diff:':>60} {d:.3e}")

print("\n=== negative control: random M (not from any builder) ===")
rng = np.random.default_rng(7)
for trial in range(3):
    Aa = rng.standard_normal((6,6)); M = Aa@Aa.T + 6*np.eye(6)
    Zr = rng.standard_normal((6,6,6)); Zr = (Zr + Zr.transpose(0,2,1) +
        Zr.transpose(1,0,2) + Zr.transpose(1,2,0) + Zr.transpose(2,0,1) +
        Zr.transpose(2,1,0))/6.0
    Q = ds.solve_Q(M, Zr); inj = ds.cum3_from_Q(Q, M)
    print(f"  trial {trial}: max|inj - Z|/|Z| = "
          f"{np.max(np.abs(inj-Zr))/np.max(np.abs(Zr)):.3e}")
