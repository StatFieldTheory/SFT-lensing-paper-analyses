"""A5: the charge head-on -- 'R1 just moves the knots until the negative
eigenvalue disappears'.  If that were true, ANY 21 break points would do.
So break the spline at deliberately WRONG places, with identical knot count and
identical sampling density, and see whether PSD is restored and whether the
density matches the parameter-free EXACT."""
import numpy as np, sys
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix/adv")
import advlib as A
from scipy.interpolate import CubicSpline
ds = A.ds; bg = A.BG; N = ds.N_COMP
E = A.Exact(background=bg, apply_c0=False)
TAB = A.TABLE_LAM


class Broken(A.R.R1Knots):
    """R1 recipe but with USER-SUPPLIED break points."""
    def __init__(self, *a, breaks=None, dx=12.0, **kw):
        super().__init__(*a, **kw); self._brk = np.asarray(breaks, float); self.dx = dx
    def _density_splines(self, cos_gamma):
        k = round(float(cos_gamma), 12)
        if k in self._cache: return self._cache[k]
        b = np.unique(np.concatenate([[self.lam_lo], self._brk, [self.lam_hi]]))
        b = b[(b >= self.lam_lo - 1e-9) & (b <= self.lam_hi + 1e-9)]
        pc = []
        for lo, hi in zip(b[:-1], b[1:]):
            n = max(6, int(np.ceil((hi-lo)/self.dx))+1)
            xs = np.linspace(lo, hi, n); f = self._sample_f(cos_gamma, xs)
            pc.append((lo, hi, [[CubicSpline(xs, f[:,a,c]) for c in range(N)] for a in range(N)]))
        self._cache[k] = pc; return pc


cos1 = float(np.cos(np.deg2rad(1.0/60)))
lamq = np.linspace(410, 2310, 601)
ex = np.array([E.matrix(cos1,float(l))[0,0] for l in lamq])
rng = np.random.default_rng(7)
inner = TAB[1:-1]
variants = [
    ("TRUE table nodes (R1)",     inner),
    ("nodes + 30 Mpc",            inner + 30.0),
    ("nodes - 30 Mpc",            inner - 30.0),
    ("cell MIDPOINTS",            0.5*(TAB[:-1]+TAB[1:])),
    ("19 uniform breaks",         np.linspace(406, 2328, 21)[1:-1]),
    ("19 random breaks (seed 7)", np.sort(rng.uniform(420, 2320, 19))),
    ("nodes jittered +-8 Mpc",    inner + rng.uniform(-8, 8, inner.size)),
    ("nodes jittered +-2 Mpc",    inner + rng.uniform(-2, 2, inner.size)),
]
med = lambda b,r: float(np.median(A.anchor(b,r)))
PAN = lambda lo,hi: A.rule_panel(lo,hi,48)
GL  = lambda lo,hi: A.rule_gl(lo,hi,256)
mE_pan = med(E, PAN)
b_cur = ds.Sigma2Builder(background=bg, apply_c0=False)
mC_gl = med(b_cur, GL); mC_pan = med(b_cur, PAN)
print(f"{'break placement':>26} {'#brk':>5} {'nonPSD/601':>11} {'minEig/maxEig':>14} "
      f"{'medErr vs EXACT':>16} {'maxErr':>10} {'GLshift':>9} {'PANshift':>9}")
for nm, brk in variants:
    b = Broken(background=bg, apply_c0=False, breaks=brk)
    v = np.array([b.matrix(cos1,float(l))[0,0] for l in lamq])
    rel = np.abs(v/ex-1)
    bad = 0; worst = 1e9
    for l in lamq:
        W=b.matrix(1.0,float(l)); X=b.matrix(cos1,float(l))
        M=np.zeros((6,6)); M[:3,:3]=W; M[3:,3:]=W; M[:3,3:]=X; M[3:,:3]=X.T
        ev=np.linalg.eigvalsh(M)
        if ev.min()<0: bad+=1
        worst=min(worst, ev.min()/max(abs(ev.max()),1e-300))
    print(f"{nm:>26} {len(np.atleast_1d(brk)):5d} {bad:11d} {worst:14.3e} "
          f"{np.median(rel):16.3e} {np.max(rel):10.3e} "
          f"{med(b,GL)/mC_gl-1:+8.4%} {med(b,PAN)/mC_pan-1:+8.4%}")

print("\nWhere are R1's only large density errors?  (0.25 Mpc scan, gamma=1')")
b_r1 = A.R.R1Knots(background=bg, apply_c0=False)
lam = np.arange(410, 2310.01, 0.25)
exf = np.array([E.matrix(cos1,float(l))[0,0] for l in lam])
v = np.array([b_r1.matrix(cos1,float(l))[0,0] for l in lam])
rel = np.abs(v/exf-1)
big = lam[rel > 0.1]
print(f"   n(>10%) = {big.size};  those lambda = {np.array2string(big, precision=2)}")
print(f"   interior table nodes in range = {np.array2string(TAB[(TAB>410)&(TAB<2310)], precision=1)}")
m = ~np.isin(lam, TAB)
print(f"   excluding the exact node points: median {np.median(rel[m]):.3e}  max {np.max(rel[m]):.3e}")
