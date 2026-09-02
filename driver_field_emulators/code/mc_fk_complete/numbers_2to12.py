"""FK as a fraction of Order-0 over 2'-12', the range current cosmic-shear
analyses actually use, for all four 2PCFs and all five source planes."""
import math
import numpy as np
P = "driver_field_emulators/products/"

def fold(path, order, a, b):
    d = np.load(P + path, allow_pickle=True)
    m = (d["a"] == a) & (d["b"] == b) & (d["order"] == order)
    g, tf = [], np.asarray(d["t_final"], float)[m] if "t_final" in d.files else None
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x/np.linalg.norm(x), y/np.linalg.norm(y)), -1, 1)))) * 60)
    return np.asarray(g), np.asarray(d["value"], float)[m], tf

GS = np.array([2.0, 5.0, 8.0, 12.0])

print("=== z_s = 5, the main analysis: FK / Order-0 [%] ===")
c = np.load("SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt/"
            "outputs/appendix_mc_curve.npz", allow_pickle=True)
go, o0 = np.asarray(c["g"], float), np.asarray(c["o0"], float)
gk, fk, _ = fold("table_permclosed_cut15360_permfix_xi.npz", 2, 0, 0)
o = np.argsort(gk); gk, fk = gk[o], fk[o]
print(f"{'gamma':>7} {'FK':>12} {'FK/O0 %':>9}")
for t in GS:
    a = float(np.interp(t, gk, fk)); b = float(np.interp(t, go, o0))
    print(f"{t:>7.1f} {a:>12.4e} {100*a/b:>9.2f}")
lo = min(100*float(np.interp(t, gk, fk))/float(np.interp(t, go, o0)) for t in GS)
hi = max(100*float(np.interp(t, gk, fk))/float(np.interp(t, go, o0)) for t in GS)
print(f"  -> kappa-kappa over 2'-12': {lo:.2f}% to {hi:.2f}%")

print("\n=== multi-z, kappa-kappa: FK / Order-0 [%] at each source plane ===")
gfk, vfk, tfk = fold("table_permclosed_cut1000_permfixmultiz_xi.npz", 2, 0, 0)
g00, v00, t00 = fold("order0_multiz_xi.npz", 0, 0, 0)
zs = sorted(set(np.round(tfk, 3)))
print(f"{'t_final':>9} {'2 arcmin':>10} {'5':>8} {'8':>8} {'12':>8}   {'range':>14}")
for t in zs:
    mf = np.abs(tfk - t) < 1e-6; m0 = np.abs(t00 - t) < 1e-6
    if mf.sum() < 4 or m0.sum() < 4: continue
    gf, vf = gfk[mf], vfk[mf]; of = np.argsort(gf); gf, vf = gf[of], vf[of]
    g0, v0 = g00[m0], v00[m0]; o0i = np.argsort(g0); g0, v0 = g0[o0i], v0[o0i]
    r = [100*float(np.interp(x, gf, vf))/float(np.interp(x, g0, v0)) for x in GS]
    print(f"{t:>9.1f} " + " ".join(f"{v:>8.2f}" for v in r) +
          f"   {min(r):.2f}-{max(r):.2f}%")
