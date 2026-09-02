"""How much does the sigma_lambda -> 0 intercept move under alternative models?

The point of this file: the intercept does NOT have to be fitted at all.  On a
fixed lattice the regulator can simply be switched off, so the limit is a
DIRECTLY COMPUTED number.  Every fitted model is then scored against it.
"""
import numpy as np

ANA = 1.59771599e-05   # paper fold, order 2, a=b=0, gamma = 1'

d = np.load("_fixN1000_g1.npz")
o = np.argsort(d["sigma"])[::-1]
s, ex_, ri = d["sigma"][o], d["exact"][o], d["ref_inj"][o]
r_ana = ex_ / ANA
r_inj = ex_ / ri

print("EXACT (noise-free) ladder, FIXED N = 1000 (the production grid), gamma = 1'")
print(f"{'sigma':>7} {'s/dlam':>7} {'exact':>15} {'exact/analytic':>15} {'exact/ref_inj':>14} "
      f"{'deficit/sig':>12}")
dl = float(d["dlam"][0])
for i in range(s.size):
    print(f"{s[i]:>7.2f} {s[i]/dl:>7.2f} {ex_[i]:>15.8e} {r_ana[i]:>15.7f} {r_inj[i]:>14.7f} "
          f"{(1-r_inj[i])/s[i]:>12.6f}")

TRUE = r_ana[s == 0.25][0]
TRUE_INJ = r_inj[s == 0.25][0]
print(f"\nDIRECTLY COMPUTED limit (sigma = 0.25, i.e. sigma/dlam = 0.13, rho = "
      f"{np.exp(-dl/0.25):.2e} -- the regulator is simply OFF):")
print(f"   exact/ref_inj  = {TRUE_INJ:.7f}   (theory says exactly 1)")
print(f"   exact/analytic = {TRUE:.7f}   <-- the number every fit below is trying to guess")


def take(*sig):
    i = [int(np.where(s == v)[0][0]) for v in sig]
    return s[i], r_ana[i]


def lin2(a, b):
    x, y = take(a, b)
    return y[1] + (y[1] - y[0]) * x[1] / (x[0] - x[1])


def polyfit_int(sigs, deg):
    x, y = take(*sigs)
    return np.polyval(np.polyfit(x, y, deg), 0.0)


def freep(sigs):
    """c0 - c1 sigma^p through three points, solved exactly."""
    x, y = take(*sigs)
    (s1, s2, s3), (y1, y2, y3) = x, y
    from scipy.optimize import brentq
    f = lambda p: (y2 - y1) * (s3**p - s2**p) - (y3 - y2) * (s2**p - s1**p)
    lo, hi = 0.2, 4.0
    if f(lo) * f(hi) > 0:
        return np.nan, np.nan
    p = brentq(f, lo, hi)
    c1 = (y2 - y1) / (s1**p - s2**p)
    return y1 + c1 * s1**p, p


print("\n" + "=" * 78)
print("ALTERNATIVE MODELS, scored against the directly computed limit")
print("=" * 78)
print(f"{'model':>46} {'intercept':>11} {'bias vs truth':>14}")
rows = [
    ("linear, 2 pts (8,4)  <-- THE PUBLISHED RECIPE", lin2(8, 4)),
    ("linear, 2 pts (16,8)", lin2(16, 8)),
    ("linear, 2 pts (4,2)", lin2(4, 2)),
    ("linear, 2 pts (2,1)", lin2(2, 1)),
    ("linear, 2 pts (12,6)", lin2(12, 6)),
    ("quadratic, 3 pts (8,4,2)", polyfit_int([8, 4, 2], 2)),
    ("quadratic, 3 pts (16,8,4)", polyfit_int([16, 8, 4], 2)),
    ("cubic, 4 pts (16,8,4,2)", polyfit_int([16, 8, 4, 2], 3)),
    ("linear LSQ, 5 pts (16,12,8,6,4)", polyfit_int([16, 12, 8, 6, 4], 1)),
    ("quadratic LSQ, 6 pts (16..3)", polyfit_int([16, 12, 8, 6, 4, 3], 2)),
]
p8 = freep([8, 4, 2])
rows.append((f"c0 - c1 sigma^p, free p = {p8[1]:.3f}, pts (8,4,2)", p8[0]))
p16 = freep([16, 8, 4])
rows.append((f"c0 - c1 sigma^p, free p = {p16[1]:.3f}, pts (16,8,4)", p16[0]))
for lbl, v in rows:
    print(f"{lbl:>46} {v:>11.6f} {100*(v-TRUE):>13.3f}%")
vals = np.array([v for _, v in rows])
print(f"\n  spread across ALL models          : {100*(vals.max()-vals.min()):.3f}%")
print(f"  every model is HIGH of the truth by : {100*(vals.min()-TRUE):.3f}% to "
      f"{100*(vals.max()-TRUE):.3f}%")
print(f"  quoted Monte-Carlo error bar        : 0.560%")
