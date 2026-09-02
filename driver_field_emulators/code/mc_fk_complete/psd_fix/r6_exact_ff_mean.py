"""Noise-free counterpart of the FF disconnected channel, plus a PSD sweep in
gamma out to 600'.

The MC recursion  s_k = resp_k s_{k-1} + F(s_{k-1}) dlam + dW_k  gives, exactly,
    <kappa>_a = dlam * sum_j Wd[j] * F_abc V[j-1]_bc ,
    V[m] = sum_{i<=m} (D_i/D_m)^4 Sigma2_i dlam ,
with Wd[j] = sum_{k>=j} w_k (D_j/D_k)^2 -- the same Wd build_grid computes.
This is the deterministic mean second-order convergence the MC estimates as
`mean_kappa_on`, and its square is the FF disconnected piece.  It is LINEAR in
Sigma2, so it isolates the repair's effect with zero Monte-Carlo noise.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R

LAM_F = bgm.LAM_SOURCE_BASELINE; LAM_LO = 406.0
bg = bgm.Background(lam_min=120.0, lam_max=LAM_F + 5.0)

def mean_kappa_exact(b, cosg, N):
    lam = np.linspace(LAM_LO, LAM_F, N); dlam = float(lam[1]-lam[0])
    D = np.asarray(bg.D(lam), float)
    w = np.full(N, dlam); w[0] = w[-1] = 0.5*dlam
    Wd = D**2 * np.cumsum((w/D**2)[::-1])[::-1]
    S = np.array([core._assemble_6x6(b, cosg, float(l)) for l in lam])
    V = np.zeros((N, 6, 6)); acc = np.zeros((6,6))
    for m in range(N):
        acc = (D[m-1]/D[m])**4 * acc if m else 0.0*acc
        acc = acc + S[m]*dlam
        V[m] = acc
    out = np.zeros(6)
    for j in range(N):
        Vp = V[j-1] if j else np.zeros((6,6))
        for r in (0,3):
            out[r]   += Wd[j]*dlam*(-(Vp[r,r]+Vp[r+1,r+1]+Vp[r+2,r+2]))
            out[r+1] += Wd[j]*dlam*(-2.0*Vp[r,r+1])
            out[r+2] += Wd[j]*dlam*(-2.0*Vp[r,r+2])
    return out, S

print("=== exact mean second-order convergence <kappa> (FF disconnected root) ===")
print(f"{'gam':>6} {'N':>6} {'stock <k>_0':>14} {'R1 <k>_0':>14} {'R1/stock':>10} "
      f"{'disc ratio':>11}")
for gam in (1.0, 10.0, 60.0, 300.0):
    for N in (600, 1000):
        c = float(np.cos(np.deg2rad(gam/60.0)))
        vals = {}
        for nm, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
            b = cls(background=bg, apply_c0=False)
            vals[nm] = mean_kappa_exact(b, c, N)[0]
        s, r = vals["stock"][0], vals["R1"][0]
        print(f"{gam:>6.1f} {N:>6} {s:>14.7e} {r:>14.7e} {r/s:>10.6f} "
              f"{(r/s)**2:>11.6f}")

print("\n=== PSD of the 6x6 Sigma2 out to gamma = 600', both builders ===")
print(f"{'gam':>7} | {'stock n_bad/N':>13} {'min eig/max':>13} | {'R1 n_bad/N':>11} {'min eig/max':>13}")
NSCAN = 1200
lam_scan = np.linspace(LAM_LO, LAM_F, NSCAN)
for gam in (0.5, 1.0, 5.0, 17.3, 60.0, 114.3, 200.0, 300.0, 450.0, 600.0):
    c = float(np.cos(np.deg2rad(gam/60.0))); row = []
    for nm, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
        b = cls(background=bg, apply_c0=False)
        nbad = 0; worst = np.inf
        for l in lam_scan:
            M = core._assemble_6x6(b, c, float(l))
            wv = np.linalg.eigvalsh(0.5*(M+M.T))
            rat = wv.min()/max(abs(wv).max(), 1e-300)
            worst = min(worst, rat)
            if wv.min() < -1e-14*abs(wv).max(): nbad += 1
        row.append((nbad, worst))
    print(f"{gam:>7.1f} | {row[0][0]:>8}/{NSCAN:<4} {row[0][1]:>13.3e} |"
          f" {row[1][0]:>6}/{NSCAN:<4} {row[1][1]:>13.3e}")
