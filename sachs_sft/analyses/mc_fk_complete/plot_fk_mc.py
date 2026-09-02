"""Diagnostic figure for this folder (NOT a paper figure)."""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(p):
    d = np.load(p, allow_pickle=True)
    return {k: d[k] for k in d.files}


def row(d, gam, sig):
    i = int(np.where((d["gamma"] == gam) & (d["sigma"] == sig))[0][0])
    return (np.asarray(d["seeds"][i], float), float(d["exact"][i]),
            float(d["analytic"][i]), float(d["vr"][i]))


d3 = load("_gate3_sigma.npz")
g4 = load("_gate4_gamma.npz")
hA, hB = load("_gate2_highstat.npz"), load("_gate2_highstat_B.npz")
hC = load("_papergrid.npz")
rA, rB = load("_reseedA.npz"), load("_reseedB.npz")

fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.3))

# --- left: sigma_lambda dependence at gamma = 1' ----------------------------
o = np.argsort(d3["sigma"]); s = d3["sigma"][o]
mc, err = [], []
for sig in s:
    if sig in (4.0, 8.0):                       # pooled 88-seed points
        a, ex, ana, _ = row(hA, 1.0, sig); b, _, _, _ = row(hB, 1.0, sig)
        c, _, _, _ = row(hC, 1.0, sig)
        v = np.concatenate([a, b, c])
    else:
        v, ex, ana, _ = row(d3, 1.0, sig)
    mc.append(v.mean() / ana); err.append(v.std(ddof=1) / np.sqrt(v.size) / ana)
ax[0].errorbar(s, mc, yerr=err, fmt="o", ms=5, color="C0", zorder=3,
               label="Monte-Carlo")
ax[0].plot(s, (d3["exact"] / d3["analytic"])[o], "s--", ms=4, color="C3",
           label="exact expectation of the estimator")
ax[0].plot(s, (d3["vr"] / d3["analytic"])[o], "^:", ms=4, color="0.45",
           label=r"$T_1$ alone (published structure)")
fit = np.polyfit(s, (d3["exact"] / d3["analytic"])[o], 1)
ss = np.linspace(0, s.max(), 50)
ax[0].plot(ss, np.polyval(fit, ss), "-", lw=0.8, color="C3", alpha=0.5)
ax[0].plot([0], [np.polyval(fit, 0)], "*", ms=11, color="C3", zorder=4)
ax[0].axhline(1.0, color="k", lw=0.8)
ax[0].set_xlabel(r"$\sigma_\lambda$  [Mpc]   (colored-noise regulator)")
ax[0].set_ylabel("ratio to analytic FK")
ax[0].set_title(r"$\gamma=1'$: removing the regulator")
ax[0].set_xlim(-1.2, 34); ax[0].set_ylim(0, 1.25)
ax[0].legend(fontsize=8, loc="center right")

# --- right: gamma sweep at sigma = 8 ---------------------------------------
gams = np.array(sorted(set(g4["gamma"])))
mc, err, ex8, ex0, vr = [], [], [], [], []
for gam in gams:
    parts = []
    e8 = ana = v8 = None
    for src in (g4, rA, rB, hA, hB, hC):
        if ((src["gamma"] == gam) & (src["sigma"] == 8.0)).any():
            a, e8, ana, v8 = row(src, gam, 8.0)
            parts.append(a)
    v = np.concatenate(parts)
    _, e4, _, _ = row(g4, gam, 4.0)
    mc.append(v.mean() / ana); err.append(v.std(ddof=1) / np.sqrt(v.size) / ana)
    ex8.append(e8 / ana); ex0.append((2 * e4 - e8) / ana); vr.append(v8 / ana)
ax[1].errorbar(gams, mc, yerr=err, fmt="o", ms=5, color="C0", zorder=3,
               label=r"Monte-Carlo, $\sigma_\lambda=8$")
ax[1].plot(gams, ex8, "s--", ms=4, color="C3",
           label=r"exact expectation, $\sigma_\lambda=8$")
ax[1].plot(gams, ex0, "-", lw=1.6, color="C2",
           label=r"exact expectation, $\sigma_\lambda\!\to\!0$")
ax[1].plot(gams, vr, "^:", ms=4, color="0.45",
           label=r"$T_1$ alone, $\sigma_\lambda=8$")
ax[1].axhline(1.0, color="k", lw=0.8)
ax[1].axhspan(0.95, 1.05, color="0.88", zorder=0)
ax[1].set_xscale("log")
ax[1].set_xlabel(r"$\gamma$  [arcmin]")
ax[1].set_ylabel("ratio to analytic FK")
ax[1].set_title("separation sweep")
ax[1].set_ylim(-0.15, 1.25)
ax[1].legend(fontsize=8, loc="center left")
for a_ in ax:
    a_.grid(alpha=0.25, lw=0.5)
fig.tight_layout()
fig.savefig("fk_mc_complete.pdf"); fig.savefig("fk_mc_complete.png", dpi=150)
print("wrote fk_mc_complete.pdf / .png")
