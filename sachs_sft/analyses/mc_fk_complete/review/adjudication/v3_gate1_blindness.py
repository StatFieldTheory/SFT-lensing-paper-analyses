"""ADJUDICATOR check 3: is the gate-1 placement identity blind to the physics?

The adversarial_tautology agent claims the identity holds at 4e-16 inside
deliberately WRONG but self-consistent systems (wrong propagator exponent,
D=1, random F, random Q, random covariance).  Verified here directly by
substituting into the same two routines gate1_placement_share.py uses.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import fk_expect_exact as ex
import gate1_placement_share as g1

N = 24; sig = 60.0; gam = 1.0
cosg = float(np.cos(np.deg2rad(gam/60.0)))
rng = np.random.default_rng(20260828)

def remake(grid, D):
    """Rebuild resp/Wd/Hd from a substituted D, exactly as build_grid does."""
    n = grid.lam.size; dlam = grid.dlam
    resp = np.empty(n); resp[0] = 1.0; resp[1:] = (D[:-1]/D[1:])**2
    w = np.full(n, dlam); w[0] = w[-1] = 0.5*dlam
    Wd = D**2 * np.cumsum((w/D**2)[::-1])[::-1]
    gj = dlam*Wd[1:]/D[:-1]**4
    t4 = np.zeros(n); t4[:-1] = np.cumsum(gj[::-1])[::-1]
    Hd = D**4*t4
    return ex.Grid(grid.lam, dlam, D, resp, Wd, Hd, grid.bg, grid.builder)

base = ex.build_grid(n_lambda=N, lam_min=406.0)
V, A, Q, rho = ex.node_stats(base, cosg, sig)

def run(tag, grid, A_, Q_, F6_):
    truth, *_ = g1.discrete_truth(grid, A_, Q_, rho, F6_)
    # expectation() hardwires F6/FS6 at module scope; patch them
    oldF, oldFS = ex.F6, ex.FS6
    ex.F6 = F6_
    ex.FS6 = F6_ + F6_.transpose(0, 2, 1)
    try:
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A_, Q=Q_, rho=rho)
    finally:
        ex.F6, ex.FS6 = oldF, oldFS
    new = (T1+T2) + (T1+T2).T
    err = np.abs(new-truth).max()/np.abs(truth).max()
    share = (T1+T1.T)[0,3]/truth[0,3]
    print(f"{tag:<38} truth_kk={truth[0,3]:>12.4e}  identity_err={err:.3e}  "
          f"T1share={share:>8.4f}")

print(f"toy grid N={N} sigma={sig} gamma={gam}'")
run("TRUE system", base, A, Q, ex.F6)
for p in (1.0, 6.0):
    Dp = base.D**(p/2.0)          # gives propagator (D_j/D_k)^p
    run(f"propagator exponent {p:.2f}", remake(base, Dp), A, Q, ex.F6)
run("D = 1 (background off)", remake(base, np.ones_like(base.D)), A, Q, ex.F6)
Frand = rng.normal(size=(6,6,6))
run("F -> random (6,6,6), no symmetry", base, A, Q, Frand)
Qr = rng.normal(size=Q.shape); Qr = 0.5*(Qr + Qr.transpose(0,1,3,2))
run("Q -> random symmetric", base, A, Qr, ex.F6)
Ar = rng.normal(size=(N,6,6)); Ar = np.einsum("kij,klj->kil", Ar, Ar)
run("A -> random SPD sequence", base, Ar, Q, ex.F6)
run("ALL wrong at once", remake(base, np.ones_like(base.D)), Ar, Qr, Frand)
