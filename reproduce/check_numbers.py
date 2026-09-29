#!/usr/bin/env python
"""Recompute every quantitative claim of the manuscript from the products on disk.

    python reproduce/check_numbers.py            # markdown table on stdout

Writes reproduce/numbers.json beside it. Each row names the claim, what the paper
says, what the products give today, and which files it came from. Run it after any
change to a product: a row that moves is a stop-and-report item for the operator,
not a licence to edit the manuscript."""
import json, math, os, subprocess, sys
from pathlib import Path
os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
from product_paths import active_manifest, resolve_product

REPO = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing")
SSA = REPO / "SFT-lensing-paper-analyses"; SACHS = SSA / "sachs_sft"
A3 = SACHS / "analyses" / "analysis3"; MC = SACHS / "analyses" / "mc_sachs_2pt"
P = SACHS / "sftwick_outputs" / "2PCF" / "cutoff_ladder"
MZ = SACHS / "sftwick_outputs" / "2PCF" / "multiz"
PERM = SACHS / "callables" / "kappa3_vertex" / "equal_time_limber_cut15360_permaware"
RB = SACHS / "callables" / "kappa3_vertex" / "rebuild"
MCF = SACHS / "analyses" / "mc_fk_complete"
ETL = SACHS / "callables" / "kappa3_vertex" / "equal_time_limber"
RUNS = SACHS / "sftwick_outputs" / "2PCF"
O0 = RUNS / "C_corr_op_O0" / "xi_C_corr_op_O0.npz"
FF = RUNS / "C_corr_op_K_limber_FF" / "xi_C_corr_op_K_limber_FF.npz"
FK_DEP = RUNS / "C_corr_op_K_limber_FK" / "xi_C_corr_op_K_limber_FK.npz"
FK_CONV = (RUNS / "C_corr_op_K_limber_FK_cut15360_permfix"
           / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz")
O0 = resolve_product("order0", O0)
FF = resolve_product("ff", FF)
FK_CONV = resolve_product("fk", FK_CONV)
P = resolve_product("cutoff_ladder", P)
MULTIZ = resolve_product("multiz", A3 / "outputs" / "multiz_kappa_2pcf_5z.npz")
MULTIZ_FF = resolve_product("multiz_ff", MZ / "multiz_ff_real_all5.npz")
MC_MARKERS = resolve_product("mc_markers", MCF / "_markers_pooled.npz")
MC_CACHE = resolve_product("mc_cache", MC / "outputs" / "appendix_mc_curve.npz")
KERNEL_CHECK = resolve_product("kernel_crosscheck", RB / "products" / "fk_kernel_crosscheck.npz")
ZETA_SLICES = resolve_product("zeta_slices", ETL / "outputs" / "zeta_bands_cut15360_gmax85.npz")
ARCMIN = 180 * 60 / np.pi
OBS = {"xi_kappa": [((0, 0), 1.0)], "xi_plus": [((1, 1), 1.0), ((2, 2), 1.0)],
       "xi_minus": [((1, 1), 1.0), ((2, 2), -1.0)], "xi_kappa_gamma_t": [((0, 1), -1.0)]}
rel = lambda p: str(Path(p).resolve().relative_to(REPO)) if str(p).startswith(str(REPO)) else str(p)


def load_pairs(path, order, t_final=None):
    d = np.load(path, allow_pickle=True)
    a, b, o = np.asarray(d["a"]), np.asarray(d["b"]), np.asarray(d["order"])
    v = np.asarray(d["value"], float); x, y = d["x"], d["y"]
    t = np.asarray(d["t_final"], float) if "t_final" in d.files else None
    keep = o == order
    if t_final is not None:
        keep &= np.isclose(t, t_final, rtol=1e-9)
    grouped, gamma = {}, None
    for pair in {(int(ai), int(bi)) for ai, bi in zip(a[keep], b[keep])}:
        m = keep & (a == pair[0]) & (b == pair[1]); idx = np.flatnonzero(m)
        g = np.array([math.acos(float(np.clip(np.dot(
            np.asarray(x[i], float) / np.linalg.norm(np.asarray(x[i], float)),
            np.asarray(y[i], float) / np.linalg.norm(np.asarray(y[i], float))), -1, 1))) * ARCMIN
            for i in idx])
        s = np.argsort(g)
        if gamma is None:
            gamma = g[s]
        grouped[pair] = v[m][s]
    return gamma, grouped


def obs(path, order, t_final=None):
    g, gr = load_pairs(path, order, t_final)
    return g, {n: sum(sg * gr[p] for p, sg in combos) for n, combos in OBS.items()
               if all(p in gr for p, _ in combos)}


rows = []
def add(claim, quoted, value, src):
    rows.append(dict(claim=claim, quoted=quoted, recomputed=value, source=src))
    print(f"| {claim} | {quoted} | {value} | {src} |")

lin = lambda g, arr, x: float(np.interp(x, g, arr))
BAND = np.array([2.0, 5.0, 8.0, 12.0])
ACTIVE = active_manifest() is not None


def quoted(current, historical):
    """Keep historical reproduction statements separate from the active revision."""
    return current if ACTIVE else historical


def first_crossing(test, gamma):
    indices = np.flatnonzero(test)
    if not indices.size:
        return "no crossing on the sampled grid"
    index = int(indices[0])
    previous = f" (previous grid point {gamma[index - 1]:.1f}')" if index else " (first sampled point)"
    return f"first grid point: {gamma[index]:.1f}'" + previous


print("| claim | quoted in paper | recomputed today | source products |")
print("|---|---|---|---|")

g, o0 = obs(O0, 0); gff, ff = obs(FF, 2); gfk, fk = obs(FK_CONV, 2); _, fkd = obs(FK_DEP, 2)
if not (np.allclose(g, gff) and np.allclose(g, gfk)):
    raise ValueError("Order-0, FF and FK angular grids differ")
assert abs(g[0] - 0.5) < 1e-6
kk = "xi_kappa"
historical_fk_value = f"{fkd[kk][0]:+.4e} (historical source plane)"
if active_manifest() is None:
    historical_fk_value += f" = {100*fkd[kk][0]/o0[kk][0]:.3f}% of historical O0"
add("Historical FK(0.5') sweep (cut1000, old callable), kk", "historical baseline only, not the active prediction",
    historical_fk_value, rel(FK_DEP))
add("FK(0.5') finite-cutoff fold (ell_max=15360, perm-aware), kk", quoted("1.91% of O0", "2.31% of O0 (insights, cutoff subsec.)"),
    f"{fk[kk][0]:+.4e} = {100*fk[kk][0]/o0[kk][0]:.2f}% of O0", rel(FK_CONV))
add("FF(0.5') vs O0, kk", quoted("2.28e-6, about 0.27% of Order-0", "few x1e-6, about 0.2% of Order-0"),
    f"FF={ff[kk][0]:+.3e}, O0={o0[kk][0]:+.3e}, FF/O0={100*ff[kk][0]/o0[kk][0]:.2f}%", rel(FF))
add("FF overtakes |O0| (kk)", quoted("first sampled crossing near 183'", "beyond ~200' (insights); ~190' (session notes)"),
    first_crossing(ff[kk] > np.abs(o0[kk]), g), rel(FF) + " + " + rel(O0))
if ACTIVE:
    add("FF overtakes |FK| (kk)", "first sampled crossing near 44'",
        first_crossing(ff[kk] > np.abs(fk[kk]), g), rel(FF) + " + " + rel(FK_CONV))
plateau = ff[kk][-8:]
add("FF plateau at large gamma (kk)", "varies by < 1% out to largest separation",
    f"FF over last 8 grid points: {plateau.min():.3e}..{plateau.max():.3e} (spread {100*(plateau.max()/plateau.min()-1):.2f}%)", rel(FF))

for name, label in (("xi_kappa", "xi_kappa"), ("xi_plus", "xi_+"), ("xi_minus", "xi_-"), ("xi_kappa_gamma_t", "xi_kappa_gamma_t")):
    r = [100 * lin(g, fk[name], x) / lin(g, o0[name], x) for x in BAND]
    inband = (g >= 2.0) & (g <= 12.0)
    rg = 100 * fk[name][inband] / o0[name][inband]
    add(f"FK/O0 over 2'-12', {label}",
        quoted({"xi_kappa": "about 1.5-1.6%", "xi_plus": "about 0.45-0.75%",
                "xi_minus": "about 14% at 2', about 1% beyond 5'",
                "xi_kappa_gamma_t": "negative, about -2.4% to -1.2%"}[name],
        {"xi_kappa": "1.3-1.7%", "xi_plus": "1.3-1.7% (nearly equal to xi_kappa)",
         "xi_minus": "~9% at 2' edge, near 1% above 3' (8.60/3.91/2.19% at 2.06/2.61/3.31')",
         "xi_kappa_gamma_t": "1-2%"}[name]),
        f"interp at 2,5,8,12: {r[0]:.2f}/{r[1]:.2f}/{r[2]:.2f}/{r[3]:.2f}%; grid pts in band: " +
        ", ".join(f"{gi:.2f}':{ri:.2f}%" for gi, ri in zip(g[inband], rg)),
        rel(FK_CONV) + " / " + rel(O0))
for name, label in (("xi_kappa", "xi_kappa"), ("xi_plus", "xi_+")):
    r = [lin(g, fk[name], x) / lin(g, ff[name], x) for x in BAND]
    add(f"FK/FF over 2'-12', {label}", quoted("several times the FF term", "several times the FF term (commit 44aade3: kk 6.9->3.5, ++ 11.9->10.9)"),
        f"at 2/5/8/12': {r[0]:.1f}/{r[1]:.1f}/{r[2]:.1f}/{r[3]:.1f}", rel(FK_CONV) + " / " + rel(FF))
for name in ("xi_minus", "xi_kappa_gamma_t"):
    r = [lin(g, fk[name], x) / lin(g, ff[name], x) for x in BAND]
    add(f"Signed FK/FF in {name} over the band", quoted("FK magnitude exceeds FF in all four observables", "FK remains above FF in all four"),
        f"FK/FF at 2/5/8/12': {r[0]:.1f}/{r[1]:.1f}/{r[2]:.1f}/{r[3]:.1f}", rel(FK_CONV) + " / " + rel(FF))

# cutoff ladder
CUTS = (960, 1920, 3840, 7680, 15360)
lad = {}
for model in ("tree", "bihalofit"):
    vals, band_vals = [], []
    for cut in CUTS:
        gg, ob = obs(P / f"table_{model}_cut{cut}_r4_xi.npz", 2)
        vals.append(100 * ob[kk][0] / o0[kk][0])
        band_vals.append([100 * lin(gg, ob[kk], x) / lin(g, o0[kk], x) for x in BAND])
    lad[model] = (np.array(vals), np.array(band_vals))
v = lad["tree"][0]; inc = np.diff(v); ratios = v[1:] / v[:-1]
add("Cutoff ladder, tree, FK/O0 at 0.5'", quoted("measured finite-cutoff values, ell_max 960 to 15360", "0.48% -> 2.31% for ell_max 960 -> 15360"),
    " -> ".join(f"{c}:{x:.2f}%" for c, x in zip(CUTS, v)) + f"; successive ratios {np.round(ratios,2).tolist()}", rel(P) + "/table_tree_cut*_r4_xi.npz")
add("Final doubling adds (0.5')", quoted("measured change from ell_max=7680 to 15360", "14%"), f"{100*(v[-1]/v[-2]-1):.1f}%", "same")
bv = lad["tree"][1]
add("Final doubling adds across 2'-12'", quoted("finite-cutoff sensitivity, not an infinite-cutoff bound", "8-9%"), ", ".join(f"{x:.0f}':{100*(bv[-1][i]/bv[-2][i]-1):.1f}%" for i, x in enumerate(BAND)), "same")
if not ACTIVE:
    r_last = inc[-1] / inc[-2]
    extrap = v[-1] + inc[-1] * r_last / (1 - r_last)
    add("Historical geometric extrapolation (0.5')", "historical estimate, not a verified bound",
        f"last-increment ratio {r_last:.2f}: {extrap:.2f}% (= {100*(extrap/v[-1]-1):.0f}% above the final finite-cutoff value)", "same")
vb = lad["bihalofit"][0]
add("BiHalofit ladder at 0.5'", quoted("measured finite-cutoff sensitivity", "keeps growing, no turnover"), " -> ".join(f"{c}:{x:.2f}%" for c, x in zip(CUTS, vb)), rel(P) + "/table_bihalofit_cut*_r4_xi.npz")

# multi-z
d = np.load(MULTIZ, allow_pickle=True)
ffr = np.load(MULTIZ_FF, allow_pickle=True)
z = np.asarray(d["z"], float); gz = np.asarray(d["gamma"], float)
o0z, fkz = np.asarray(d["o0"], float), np.asarray(d["fk"], float); ffz = np.asarray(ffr["ff"], float)
half = [100 * fkz[i, 0] / o0z[i, 0] for i in range(len(z))]
add("Multi-z FK/O0 at 0.5' (z_s=1,1.7,2.5,3.2,4)", quoted("1.09/1.44/1.66/1.77/1.85%, increasing toward 1.91% at z_s=5", "session note: 1.36/1.77/2.04/2.14/2.24%, monotone toward 2.31% at z=5"),
    "/".join(f"{x:.2f}" for x in half) + "%", rel(MULTIZ))
bz = [[100 * lin(gz, fkz[i], x) / lin(gz, o0z[i], x) for x in BAND] for i in range(len(z))]
add("Multi-z FK/O0 over 2'-12' by z_s", quoted("about 0.86-0.94% at z_s=1, increasing with source redshift", "0.9-1.1% at z_s=1 (intro); grows with source redshift"),
    "; ".join(f"z={z[i]:g}: {min(b):.2f}-{max(b):.2f}%" for i, b in enumerate(bz)), "same")
ffplat = [ffz[i, -5:].mean() for i in range(len(z))]
add("Multi-z FF plateaus (real per-z sweeps)", quoted("large-separation plateau, without an independent mean-square accuracy claim", "plateau = <kappa>^2(z_s) to 0.1-2.2%"),
    "; ".join(f"z={z[i]:g}: {ffplat[i]:.3e}" for i in range(len(z))) + " (mean over last 5 gammas)", rel(MULTIZ_FF))

# harmonic space: B/E, EB, C_kk correction
sys.path.insert(0, str(A3))
import plot_analysis3_cl_decomposition as C  # noqa: E402
C.FK_NPZ = FK_CONV
ELL = C.ELL
g_fk, fkg = C.load_sweep_order(C.FK_NPZ, 2); g_ff, ffg = C.load_sweep_order(C.FF_NPZ, 2); g0, o0g = C.load_sweep_order(C.O0_NPZ, 0)
s22 = C.build_curved_matrix(g_fk, ELL, 2, 2); s2m2 = C.build_curved_matrix(g_fk, ELL, 2, -2); s00 = C.build_curved_matrix(g0, ELL, 0, 0)
def pol(gr):
    pbb = C.forward_curved(C._combine(gr, [((1, 1), 1.0), ((2, 2), 1.0)]), s22, dc_subtract=False)
    mbb = C.forward_curved(C._combine(gr, [((1, 1), 1.0), ((2, 2), -1.0)]), s2m2, dc_subtract=False)
    eb = C.forward_curved(C._combine(gr, [((1, 2), 1.0)]), s2m2, dc_subtract=False)
    return 0.5 * (pbb + mbb), 0.5 * (pbb - mbb), eb
ee_fk, bb_fk, eb_fk = pol(fkg); ee_ff, bb_ff, eb_ff = pol(ffg)
i60, i1500 = int(np.argmin(np.abs(ELL - 60))), int(np.argmin(np.abs(ELL - 1500)))
band = (ELL >= 50) & (ELL <= 1500)
add("FK recovered BB/EE numerical residual at ell=60 and ell=1500", "physical FK BB vanishes at this order, finite-grid recovery leaves a residual",
    f"ell={ELL[i60]}: {bb_fk[i60]/ee_fk[i60]:.3f}; ell={ELL[i1500]}: {bb_fk[i1500]/ee_fk[i1500]:.3f}; median over 50<=ell<=1500: {np.median(bb_fk[band]/ee_fk[band]):.3f}",
    rel(FK_CONV) + " via plot_analysis3_cl_decomposition transform")
add("FF recovered polarization diagnostic", "FF can generate BB, the recovered ratio is not an accuracy-certified physical amplitude",
    f"median EE/BB over band: {np.median(ee_ff[band]/bb_ff[band]):.1f} (BB/EE median {np.median(bb_ff[band]/ee_ff[band]):.4f})", rel(FF))
add("EB parity null", "parity-even inputs give zero EB, recovered FF residual is numerical",
    f"max|EB_FK|={np.max(np.abs(eb_fk)):.2e}; max|EB_FF|/max|EE_FF|={np.max(np.abs(eb_ff))/np.max(np.abs(ee_ff)):.1e}", "same")
o0kk = C._combine(o0g, [((0, 0), 1.0)]); ffkk = C._combine(ffg, [((0, 0), 1.0)]); fkkk = C._combine(fkg, [((0, 0), 1.0)])
cl0 = C.forward_curved(o0kk, s00, dc_subtract=True)
clf = C.forward_curved(o0kk + ffkk + fkkk, s00, dc_subtract=True)
corr = 100 * (clf - cl0) / cl0
add("Recovered C_kk total Order-2 correction over 50<=ell<=1500", quoted("about 1.5-2.0%, finite-angle transform diagnostic", "slowly rising one to two percent"),
    f"{corr[band].min():.2f}% (ell={ELL[band][np.argmin(corr[band])]}) .. {corr[band].max():.2f}% (ell={ELL[band][np.argmax(corr[band])]})", "same")
pbb0 = C.forward_curved(C._combine(o0g, [((1, 1), 1.0), ((2, 2), 1.0)]), s22, dc_subtract=False)
pbbf = C.forward_curved(C._combine(o0g, [((1, 1), 1.0), ((2, 2), 1.0)]) + C._combine(ffg, [((1, 1), 1.0), ((2, 2), 1.0)]) + C._combine(fkg, [((1, 1), 1.0), ((2, 2), 1.0)]), s22, dc_subtract=False)
corr2 = 100 * (pbbf - pbb0) / pbb0
add("Recovered C_EE+BB total Order-2 correction over 50<=ell<=1500", quoted("about 0.4-1.4%, finite-angle transform diagnostic", "one to two percent"), f"{corr2[band].min():.2f}% .. {corr2[band].max():.2f}%", "same")

# validation numbers
mk = np.load(MC_MARKERS)
if active_manifest() is None:
    fold_at = np.exp(np.interp(np.log(mk["gamma"]), np.log(g), np.log(fk[kk])))
else:
    fold_at = np.interp(mk["gamma"], g, fk[kk])
add("FK Monte-Carlo markers / fold", quoted("ratios within about 1%, including finite-sample and discretization effects", "about a percent at 1' (CHANGES: 1.012+/-0.005, 1.012+/-0.008, 1.008+/-0.014, 0.999+/-0.033)"),
    "; ".join(f"{gm:g}': {f/a:.3f}+/-{e/a:.3f}" for gm, f, e, a in zip(mk["gamma"], mk["fk"], mk["err"], mk["analytic"])) +
    f" (markers' stored analytic vs today's fold interp: {np.max(np.abs(mk['analytic']/fold_at-1)):.1e})",
    rel(MC_MARKERS))
if active_manifest() is None:
    hl = subprocess.run([sys.executable, "headline.py"], cwd=MCF, capture_output=True, text=True).stdout
    pooled_headline = hl.strip().splitlines()[-2].strip()
else:
    first = int(np.argmin(np.abs(mk["gamma"] - 1.0)))
    pooled_headline = f"{mk['fk'][first] / mk['analytic'][first]:.4f} +/- {mk['err'][first] / mk['analytic'][first]:.4f}"
add("Pooled sigma_lambda->0 MC/analytic FK at 1'", "reproduces to about a percent", pooled_headline, rel(MC_MARKERS))
cache = np.load(MC_CACHE, allow_pickle=True)
ffm = np.interp(cache["g_mc"], cache["g_ff"] if "g_ff" in cache.files else cache["g"], cache["ff_mom"])
pull = np.abs(cache["ffc_mc"] - ffm) / cache["ffc_se"]
add("Local-covariance FF Monte-Carlo diagnostic", quoted("subpercent agreement, maximum pull about 2.4 standard errors at N=4000", "within 1.3 standard errors at every separation"),
    f"max pull {pull.max():.2f} at {cache['g_mc'][np.argmax(pull)]:g}' against the interpolated dense reference", rel(MC_CACHE))
if KERNEL_CHECK.suffix == ".json":
    kernel_record = json.loads(KERNEL_CHECK.read_text())
    finest = int(np.argmax(kernel_record["n_lambda"]))
    kernel_values = np.asarray(kernel_record["local_fk"], float)[finest]
    kc = {"gamma": kernel_record["gamma"],
          "ratio": kernel_values / np.asarray(kernel_record["fold_values"], float)}
else:
    kc = np.load(KERNEL_CHECK)
add("Kernel quadrature cross-check / fold", quoted("within 1% at the four sampled angles (1', 2.6', 6.7', 17.3')", "about a percent at arcminutes, within 7% out to 17'"),
    ", ".join(f"{gm:g}':{r:.3f}" for gm, r in zip(kc["gamma"], kc["ratio"])), rel(KERNEL_CHECK))

# zeta slices sign change
zb = np.load(ZETA_SLICES, allow_pickle=True)
gz2 = np.asarray(zb["gamma_arcmin"], float); zs = np.asarray(zb["z_shells"], float)
out = []
for ch in ("TTT", "TTP", "Bmod", "Dmod"):
    arr = np.asarray(zb[f"full_{ch}"], float)   # (n_gamma, n_shell)?
    if arr.shape[0] != gz2.size:
        arr = arr.T
    flips = []
    for j in range(arr.shape[1]):
        s = np.sign(arr[:, j]); idx = np.flatnonzero(s[1:] * s[:-1] < 0)
        if idx.size:
            index = int(idx[0])
            flips.append(f"z={zs[j]:.1f} between {gz2[index]:.0f}' and {gz2[index + 1]:.0f}'")
    out.append(f"{ch}: " + (", ".join(flips) if flips else "no sign change"))
add("zeta slices: sign-change brackets in the driving-field shells", "channel-dependent sign changes, shell redshift rather than source redshift", " | ".join(out), rel(ZETA_SLICES))

# linear n_eff crossing
pk = np.loadtxt("/Users/zzhang/projects/angular_statistics/canoes/examples/data/PCAMBz0.txt")
k, pkv = pk[:, 0], pk[:, 1]
neff = np.gradient(np.log(pkv), np.log(k))
kc2 = k[np.flatnonzero(neff < -2)[0]]
add("linear P(k): n_eff crosses -2 at", "k=0.21 h/Mpc (nonlinear one-halo delays it to 8.4 h/Mpc)", f"{kc2:.3f} h/Mpc (finite-difference n_eff on PCAMBz0.txt); nonlinear value not re-derived here (note/sec_uv.tex)", "canoes/examples/data/PCAMBz0.txt")
add("Fig 12 truncation angle", "83 deg", f"gamma_max = {g[-1]/60:.2f} deg", rel(O0))

Path(__file__).resolve().parent.joinpath("numbers.json").write_text(json.dumps(rows, indent=1))
print(f"\n[{len(rows)} rows written to reproduce/numbers.json]")
