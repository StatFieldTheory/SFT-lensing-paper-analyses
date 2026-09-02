import numpy as np, _bootstrap
ds, bgm, core = _bootstrap.wire()
gam=1.0; cosg=float(np.cos(np.deg2rad(gam/60.0))); sig=8.0
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE+5.0)
builder = ds.Sigma2Builder(background=bg, apply_c0=False)
for lam in (500.0, 1200.0, 2000.0):
    S2 = core._assemble_6x6(builder, cosg, lam); V = S2/(2*sig)
    Z6 = ds.zeta6(cosg, lam)
    sym = (Z6 + Z6.transpose(0,2,1) + Z6.transpose(1,0,2) + Z6.transpose(1,2,0)
           + Z6.transpose(2,0,1) + Z6.transpose(2,1,0))/6.0
    asym = np.abs(Z6-sym).max()/np.abs(Z6).max()
    Zt = Z6/(2*sig)**2
    Q = ds.solve_Q(V, Zt); rt = ds.cum3_from_Q(Q, V)
    # the FK contraction actually used: -sum_c zeta_{cc3}
    tru = -sum(Z6[c,c,3] for c in range(3))
    inj = -sum(rt[c,c,3] for c in range(3))*(2*sig)**2
    print(f"lam={lam:7.1f} asym={asym:.3e} rt_relerr={np.abs(rt-Zt).max()/np.abs(Zt).max():.4f} "
          f"contract true={tru:.5e} injected={inj:.5e} ratio={inj/tru:.5f}")
    # rcond sensitivity
    for rc in (1e-6, 1e-10, 1e-14):
        Q2 = ds.solve_Q(V, Zt, rcond=rc); rt2 = ds.cum3_from_Q(Q2, V)
        inj2 = -sum(rt2[c,c,3] for c in range(3))*(2*sig)**2
        print(f"     rcond={rc:.0e} -> ratio={inj2/tru:.5f}  |Q|max={np.abs(Q2).max():.3e}")
