"""ADVERSARIAL TEST A4 -- the pieces the fold comparison cannot see, tested by sampling.

(1) INJECTION FIDELITY.  Under the "smeared" calibration the field is BUILT to carry
    the tabulated zeta, so zeta_injected == sym(zeta_tab) is a near-tautology.  Measure
    how near, because that number is the ceiling on what any sigma->0 comparison can test.

(2) THE 0.5 IN THE DEFORMATION, solve_Q AND THE C-SUBTRACTION, tested by BRUTE-FORCE
    SAMPLING and no Wick algebra at all.  Draw z ~ N(0, B), form f_h = z + h*0.5*Q:(zz - C),
    and evaluate the O(h) coefficient of the sample third moment by a central difference
    of the ACTUAL products with common random numbers.  Compare to the tabulated zeta.
    A wrong 0.5, a wrong solve_Q convention, or the wrong subtracted constant all fail.

(3) THE A<->B SYMMETRISATION.  Is T1+T2 already symmetric (in which case "+transpose" is
    a doubling and the factor 2 is a free choice) or genuinely not (in which case the
    symmetrisation is forced)?
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex

DS = ex._DS
N, SIG = 1000, 8.0
GAM = 1.01546
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
grid = ex.build_grid(n_lambda=N)
V, A, Q, rho = ex.node_stats(grid, cosg, SIG, calibrate="smeared")
B = ex.calib_matrix(grid, V, A, rho, SIG, "smeared")
zi = ex.zeta_injected(B, Q)
zt = np.array([DS.zeta6(cosg, float(l)) for l in grid.lam])
zs = (zt + zt.transpose(0, 1, 3, 2) + zt.transpose(0, 2, 1, 3)
      + zt.transpose(0, 2, 3, 1) + zt.transpose(0, 3, 1, 2)
      + zt.transpose(0, 3, 2, 1)) / 6.0

print("=" * 92)
print("A4.1  injection fidelity at the production settings (N=1000, sigma=8, gamma=1.015')")
print("=" * 92)
sc = np.max(np.abs(zt), axis=(1, 2, 3))
r_sym = np.max(np.abs(zi - zs), axis=(1, 2, 3)) / sc
r_raw = np.max(np.abs(zi - zt), axis=(1, 2, 3)) / sc
r_asym = np.max(np.abs(zs - zt), axis=(1, 2, 3)) / sc
print(f"  |zeta_inj - sym(zeta_tab)| / |zeta_tab| : median {np.median(r_sym):.3e} "
      f"max {r_sym.max():.3e}")
print(f"  |zeta_inj - zeta_tab|      / |zeta_tab| : median {np.median(r_raw):.3e} "
      f"max {r_raw.max():.3e}")
print(f"  |sym(zeta_tab) - zeta_tab| / |zeta_tab| : median {np.median(r_asym):.3e} "
      f"max {r_asym.max():.3e}    (the table's leg asymmetry)")
kern = grid.dlam * grid.Wd * grid.Hd
for nm, z in (("tabulated (raw legs)", zt), ("symmetrised", zs), ("injected", zi)):
    M = np.einsum("m,Abc,mbcB->AB", kern, ex.F6, z, optimize=True)
    print(f"  fk_reference[{nm:>20}] = {(M+M.T)[0,3]:.6e}")
print("  => under the smeared calibration the sigma->0 answer is FIXED to the")
print("     symmetrised-table fold by construction; a sigma->0 gate can only test the KERNEL.")

print()
print("=" * 92)
print("A4.2  BRUTE FORCE: does f = z + 0.5 Q:(z z - C) actually carry the tabulated zeta?")
print("      (central difference in the deformation amplitude, common random numbers,")
print("       no Wick contraction anywhere in the estimator)")
print("=" * 92)
rng = np.random.default_rng(4242)
for k in (120, 500, 880):
    Bk, Qk = B[k], Q[k]
    d, U = np.linalg.eigh(Bk)
    d = np.clip(d, 0.0, None)
    L = U @ np.diag(np.sqrt(d))
    m = 4_000_000
    xi = rng.standard_normal((6, m))
    z = L @ xi
    ww = np.einsum("im,jm->ijm", z, z, optimize=True) - Bk[:, :, None]
    g = 0.5 * np.einsum("auv,uvm->am", Qk, ww, optimize=True)
    sg = np.sqrt(np.mean(g[0] ** 2)); sz = np.sqrt(np.mean(z[0] ** 2))
    h = 1e-3 * sz / sg
    fp = z + h * g
    fm = z - h * g
    # sample third moments of the ACTUAL fields, then central difference
    def m3(f, idx):
        a, b, c = idx
        return float(np.mean(f[a] * f[b] * f[c]))
    print(f"  node k={k}  lam={grid.lam[k]:.1f}  h={h:.3e}   (m = {m:,} samples)")
    print(f"    {'(a,b,c)':>10} {'brute-force O(h)':>18} {'sym(zeta_tab)':>16} {'ratio':>9}"
          f" {'MC err est':>12}")
    for idx in ((0, 0, 0), (0, 0, 3), (0, 3, 3), (3, 3, 3), (0, 1, 1), (1, 1, 3)):
        a, b, c = idx
        est = (m3(fp, idx) - m3(fm, idx)) / (2 * h)
        # per-sample estimator variance for the error bar
        per = g[a]*z[b]*z[c] + z[a]*g[b]*z[c] + z[a]*z[b]*g[c]
        err = float(np.std(per)) / math.sqrt(m)
        tgt = zs[k][idx]
        print(f"    {str(idx):>10} {est:18.6e} {tgt:16.6e} {est/tgt:9.4f} {err:12.2e}")
    # falsification controls at one component
    idx = (0, 0, 3)
    a, b, c = idx
    est = (m3(fp, idx) - m3(fm, idx)) / (2 * h)
    print(f"    controls at {idx}:  measured/target = {est/zs[k][idx]:.4f}")
    print(f"       if the code's 0.5 were 1.0   -> would read {2*est/zs[k][idx]:.4f}")
    print(f"       if Q solved against Sigma2   -> ratio "
          f"{ (lambda Qn: (0.5*np.einsum('auv,uvm->am',Qn,ww,optimize=True))) and 0:.0f}", end="")
    Qn = DS.solve_Q(V[k] * (2.0 * SIG), np.ascontiguousarray(zt[k]))
    gn = 0.5 * np.einsum("auv,uvm->am", Qn, ww, optimize=True)
    hn = 1e-3 * sz / np.sqrt(np.mean(gn[0] ** 2))
    fpn, fmn = z + hn * gn, z - hn * gn
    estn = (m3(fpn, idx) - m3(fmn, idx)) / (2 * hn)
    print(f"  -> {estn/zs[k][idx]:.4f}   (nominal calibration, same sampling)")
    del z, ww, g, fp, fm, xi

print()
print("=" * 92)
print("A4.3  is T1+T2 already symmetric?  (if it were, '+transpose' would be a bare factor 2)")
print("=" * 92)
T1, T2 = ex.expectation(grid, cosg, SIG, V=V, A=A, Q=Q, rho=rho)
S = T1 + T2
print(f"  (T1+T2)[0,3] = {S[0,3]:.6e}    (T1+T2)[3,0] = {S[3,0]:.6e}"
      f"    ratio = {S[0,3]/S[3,0]:.6f}")
print(f"  antisymmetric part / symmetric part (full 6x6, max) = "
      f"{np.max(np.abs(S-S.T))/np.max(np.abs(S+S.T)):.4f}")
print(f"  T1[0,3] = {T1[0,3]:.6e}  T1[3,0] = {T1[3,0]:.6e}   ratio {T1[0,3]/T1[3,0]:.6f}")
print(f"  T2[0,3] = {T2[0,3]:.6e}  T2[3,0] = {T2[3,0]:.6e}   ratio {T2[0,3]/T2[3,0]:.6f}")
print("  sharp form check  M = T1 + T2^T  (the placement proof's unsymmetrised statement):")
print(f"     (T1+T2^T)[0,3] = {(T1+T2.T)[0,3]:.6e} ; (T1+T2^T)[3,0] = {(T1+T2.T)[3,0]:.6e}")
