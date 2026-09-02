"""Do the teeth survive in the MONTE-CARLO arm, not just in the Wick evaluator?

Each defect is injected into the simulated ODE (mc_inj), and the MC is compared
BOTH to the paper fold (does the comparison catch it?) and to the exact
expectation of the same defective estimator (is the MC faithfully solving the
defective dynamics?).
"""
import _env, time, numpy as np
import fk_expect_exact as ex, inj_exact as ij, mc_inj
from run_gate2 import analytic_fk

GAM, SIG, N, NR, SEEDS = 1.0, 8.0, 1000, 16000, 6
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
ana = float(analytic_fk([GAM])[0])
grid = ex.build_grid(n_lambda=N)
gp1 = ij.build_grid(n_lambda=N, resp_pow=1.0)
nodes = ex.node_stats(grid, cosg, SIG)
V, A, Q, rho = nodes
MAP = {"none": ("none", grid), "F_new_state": ("acausal_Phi", grid),
       "F_no_leg_sym": ("F_no_leg_sym", grid), "resp_pow1": (None, gp1)}
print(f"gamma={GAM}' sigma={SIG} N={N} n_real={NR} seeds={SEEDS}  fold={ana:.6e}")
print(f"{'defect in the ODE':<18}{'MC':>13}{'MC/fold':>10}{'exact/fold':>12}"
      f"{'MC/exact':>10}{'seedSEM':>10}")
for d, (xd, g) in MAP.items():
    if xd is None:
        T1, T2 = ij.expectation(g, cosg, SIG, V, A, Q, rho)
    else:
        T1, T2 = ij.expectation(g, cosg, SIG, V, A, Q, rho, defect=xd)
    exv = ij.fk_kk(T1, T2, xd or "none")
    mcg = grid if d != "resp_pow1" else gp1
    vals = []
    t0 = time.time()
    for s in range(SEEDS):
        vals.append(mc_inj.simulate(cosg, mcg, nodes, SIG, NR, 2000,
                                    20260828 + 7919 * s, defect=d))
    v = np.array(vals); sem = v.std(ddof=1) / np.sqrt(v.size)
    print(f"{d:<18}{v.mean():>13.5e}{v.mean()/ana:>10.4f}{exv/ana:>12.4f}"
          f"{v.mean()/exv:>10.4f}{sem/abs(v.mean()):>9.2%}  [{time.time()-t0:.0f}s]",
          flush=True)
