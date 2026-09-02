"""ADVERSARIAL TEST A3 -- the central experiment.

CLAIM UNDER ATTACK: "the Monte-Carlo and the exact expectation agree only because
they share the author's conventions and therefore share his bugs; the agreement
measures nothing."

The experiment.  fk_complete_core.simulate_fk_complete and fk_expect_exact.expectation
both accept a Grid object.  Every convention that lives in the Grid (the background
D, the response ratios resp, the LOS window Wd, the F-leg response Hd) is therefore
literally SHARED: pass a doctored Grid to both and neither can notice.  So doctor
one convention at a time, consistently, and measure

    MC / exact          -- does the MC-vs-exact agreement notice?    (prediction: no)
    exact / paper fold  -- does the comparison to the INDEPENDENT object notice?

Doctoring is done by replacing the background D itself, so resp, Wd and Hd all
follow consistently and the doctored world is a self-consistent stochastic model
with the wrong propagator -- exactly the "shared bug" the attack posits.

Everything is evaluated AT AN EXACT NODE of the paper's fold table, so that the
linear-in-gamma interpolation of run_gate2.analytic_fk plays no role.
"""
import sys, os, math, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex
import fk_complete_core as fc
import _bootstrap

# ---- fold nodes, no interpolation ------------------------------------------
d = np.load(_bootstrap.FOLD, allow_pickle=True)
m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
gs, vsf = [], np.asarray(d["value"], float)[m]
for x, y in zip(d["x"][m], d["y"][m]):
    x = np.asarray(x, float); y = np.asarray(y, float)
    gs.append(math.degrees(math.acos(float(np.clip(
        np.dot(x/np.linalg.norm(x), y/np.linalg.norm(y)), -1, 1)))) * 60)
gs = np.asarray(gs); o = np.argsort(gs); gs = gs[o]; vsf = vsf[o]
print("fold nodes near 1':", np.array2string(gs[:5], precision=5),
      "values", np.array2string(vsf[:5], precision=6))

GAM = float(gs[1])                      # exact table node ~ 1.01546'
FOLDV = float(vsf[1])
print(f"\nUsing the EXACT fold node gamma = {GAM:.6f}'  fold = {FOLDV:.6e}")
print(f"  (run_gate2.analytic_fk interpolates LINEARLY in gamma; at 1.0' that gives "
      f"{np.interp(1.0, gs, vsf):.6e}, log-log interp gives "
      f"{math.exp(np.interp(math.log(1.0), np.log(gs), np.log(np.abs(vsf)))):.6e}"
      f"  -> a {100*(np.interp(1.0,gs,vsf)/math.exp(np.interp(math.log(1.0),np.log(gs),np.log(np.abs(vsf))))-1):.2f}% "
      f"choice in the denominator of the headline number)")

N = 1000
SIG = 8.0
NREAL = 16_000
SEEDS = (20260901, 20260902)
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))

base = ex.build_grid(n_lambda=N)
t0 = time.time()
nodes = ex.node_stats(base, cosg, SIG, calibrate="smeared")
V, A, Q, rho = nodes
print(f"node_stats: {time.time()-t0:.1f}s   (independent of the doctored D)")
Mcal = ex.calib_matrix(base, V, A, rho, SIG, "smeared")
zi = ex.zeta_injected(Mcal, Q)


def doctor(power):
    """A self-consistent world in which the linear response is (D_j/D_k)^(2*power)."""
    D = base.D ** power
    resp = np.empty(N); resp[0] = 1.0; resp[1:] = (D[:-1]/D[1:]) ** 2
    w = np.full(N, base.dlam); w[0] = w[-1] = 0.5*base.dlam
    tail = np.cumsum((w/D**2)[::-1])[::-1]
    Wd = D**2 * tail
    gj = base.dlam * Wd[1:] / D[:-1]**4
    t4 = np.zeros(N); t4[:-1] = np.cumsum(gj[::-1])[::-1]
    Hd = D**4 * t4
    return ex.Grid(base.lam, base.dlam, D, resp, Wd, Hd, base.bg, base.builder)


worlds = {
    "baseline: response (D_j/D_k)^2": 1.0,
    "shared bug: response ^2 -> ^2.06": 1.03,
    "shared bug: response ^2 -> ^2.4": 1.2,
    "shared bug: response ^2 -> ^4": 2.0,
    "shared bug: response ^2 -> ^1": 0.5,
}
print()
print("=" * 112)
print(f"{'world':<34} {'exact':>12} {'MC':>12} {'MC/exact':>12} {'exact/fold':>11}"
      f" {'ref0/fold':>10} {'T1/(T1+T2)':>11}")
print("=" * 112)
for name, p in worlds.items():
    g = doctor(p)
    T1, T2 = ex.expectation(g, cosg, SIG, V=V, A=A, Q=Q, rho=rho)
    tot = (T1 + T2) + (T1 + T2).T
    exact = float(tot[0, 3])
    t1o = float((T1 + T1.T)[0, 3])
    ref0 = float(ex.fk_reference(g, zi)[0, 3])
    mcs = []
    for s in SEEDS:
        r = fc.simulate_fk_complete(GAM, SIG, n_lambda=N, n_real=NREAL,
                                    batch_size=2000, seed=s, grid=g, nodes=nodes)
        mcs.append(r.kk)
    mc = float(np.mean(mcs)); sd = float(np.std(mcs, ddof=1)/np.sqrt(len(mcs)))
    print(f"{name:<34} {exact:12.5e} {mc:12.5e} {mc/exact:8.4f}+-{sd/abs(exact):.4f}"
          f" {exact/FOLDV:11.4f} {ref0/FOLDV:10.4f} {t1o/exact:11.4f}")
print("=" * 112)
print("READ: MC/exact stays at 1 in every world (the two share the Grid, so neither")
print("      can see the doctoring).  exact/fold moves, i.e. the comparison that")
print("      carries the physics is the one against the INDEPENDENT fold, not MC-vs-exact.")
