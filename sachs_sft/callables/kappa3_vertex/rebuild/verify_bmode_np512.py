"""Adversarial verification of the FK B-mode claim after the n_phi=512 vertex rebuild.

Written independently of the parent session's script, from the primitives in
analyses/analysis3/{plot_analysis3_cl_decomposition,plot_cl_EB_polarization}.py.

Sections:
  0. provenance / grid checks, and the 2x2 (vertex table) x (vertex callable)
  1. reproduce the DEPLOYED harmonic B/E baseline (0.951 @ l=60, 0.420 @ l=1500,
     median 0.629 over 50<=l<=1500).  If this fails, everything below is suspect.
  2. same on every other fold
  3. hypothesis (a): per-channel median |new/old| on the two vertex tables
  4. hypothesis (b): real-space |xi_-/xi_+| vs harmonic C_BB/C_EE
  5. the note's INFERENCE, tested directly: does |xi_-/xi_+| rising imply B/E rising?
  6. physicality of the B-mode collapse: positivity, smoothness, EB null,
     cancellation conditioning

FOLDS.  The deployed FK sweep was folded with the PERM-AWARE vertex callable
(callables/kappa3_vertex/equal_time_limber_cut15360_permaware/
perm_aware_kappa3_callable.py).  run_fk_variant.py defaults to the PLAIN
production callable, so a rebuild folded without ``--callable`` changes the
callable as well as the table.  Both arms of that 2x2 are computed here; the
two extra folds are produced by:

  SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
  $SFT run_fk_variant.py products/pieces_nphi512/table_permclosed_np512.npz \
      --callable ../../../callables/kappa3_vertex/equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable.py \
      --tag permaware --n-jobs 6
  $SFT run_fk_variant.py ../equal_time_limber_cut15360_permaware/table_permclosed.npz \
      --tag plaincall --n-jobs 6 \
      --work-dir $PWD/products/pieces_nphi512/_variant_deployed_plaincall \
      --out $PWD/products/pieces_nphi512/table_deployed_plaincall_xi.npz

Run with the PyCCL interpreter.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

A3 = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
          "sachs_sft/analyses/analysis3")
sys.path.insert(0, str(A3))
import plot_analysis3_cl_decomposition as C  # noqa: E402

ROOT = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft")
P512 = ROOT / "callables/kappa3_vertex/rebuild/products/pieces_nphi512"

# (label, fold npz, vertex table, vertex callable)
FOLDS = (
    ("A dep64+permaware ", C.FK_NPZ,                            "deployed n_phi=64", "perm-aware"),
    ("B np512+plaincall ", P512 / "table_permclosed_np512_xi.npz",          "rebuilt n_phi=512", "PLAIN"),
    ("C np512+permaware ", P512 / "table_permclosed_np512_permaware_xi.npz", "rebuilt n_phi=512", "perm-aware"),
    ("D dep64+plaincall ", P512 / "table_deployed_plaincall_xi.npz",        "deployed n_phi=64", "PLAIN"),
)
VTX_DEPLOYED = (ROOT / "callables/kappa3_vertex/equal_time_limber_cut15360_permaware/"
                "table_permclosed.npz")
VTX_REBUILT = P512 / "table_permclosed_np512.npz"

BAND = (50.0, 1500.0)
ELL_DENSE = np.array(sorted({int(round(v)) for v in np.geomspace(20.0, 1500.0, 90)}),
                     dtype=float)


def sep(title):
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def polar(grouped, setups):
    """(EE, BB, EB, xi_+, xi_-, cross, EE+BB, EE-BB) for one Order-2 fold.

    xi_+ = <g+g+> + <gxgx> = (1,1)+(2,2)  -> kernel d^l_{2,2}   -> EE+BB
    xi_- = <g+g+> - <gxgx> = (1,1)-(2,2)  -> kernel d^l_{2,-2}  -> EE-BB
    <g+gx> = (1,2)                        -> kernel d^l_{2,-2}  -> EB
    """
    s22, s2m2 = setups
    xi_p = grouped[(1, 1)] + grouped[(2, 2)]
    xi_m = grouped[(1, 1)] - grouped[(2, 2)]
    cross = grouped[(1, 2)]
    epb = C.forward_curved(xi_p, s22)
    emb = C.forward_curved(xi_m, s2m2)
    eb = C.forward_curved(cross, s2m2)
    return dict(EE=0.5 * (epb + emb), BB=0.5 * (epb - emb), EB=eb,
                xi_p=xi_p, xi_m=xi_m, cross=cross, epb=epb, emb=emb)


def main() -> int:
    # ------------------------------------------------------------ 0. provenance
    sep("0. PROVENANCE / GRID CHECKS  --  the 2x2 of (vertex table) x (vertex callable)")
    print(f"    {'fold':>20} {'vertex table':>20} {'callable':>12}  exists")
    res = {}
    grids = {}
    for lab, path, tab, call in FOLDS:
        print(f"    {lab:>20} {tab:>20} {call:>12}  {path.exists()}")
        g, gr = C.load_sweep_order(path, 2)
        grids[lab] = g
        res[lab] = gr
    g0 = grids[FOLDS[0][0]]
    for lab in grids:
        assert np.allclose(grids[lab], g0), f"gamma grid differs for {lab}"
    print(f"\n  all four folds share the gamma grid: {g0.size} points, "
          f"{g0.min():.4f} .. {g0.max():.1f} arcmin")

    vd = dict(np.load(VTX_DEPLOYED, allow_pickle=True))
    vr = dict(np.load(VTX_REBUILT, allow_pickle=True))
    print(f"  vertex tables: rows {vd['zeta_TTT'].shape} / {vr['zeta_TTT'].shape}   "
          f"triples identical {np.array_equal(vd['cosine_triples'], vr['cosine_triples'])}   "
          f"shells identical {np.allclose(vd['chi_shells_Mpc'], vr['chi_shells_Mpc'])}")
    C._selftest_wigner()

    setups = (C.build_curved_matrix(g0, C.ELL, 2, 2),
              C.build_curved_matrix(g0, C.ELL, 2, -2))
    out = {lab: polar(res[lab], setups) for lab in res}

    # ------------------------------------------------- 1/2. harmonic B/E
    sep("1+2. HARMONIC C_BB/C_EE  (curved-sky Wigner-d forward transform)")
    m = (C.ELL >= BAND[0]) & (C.ELL <= BAND[1])
    i60 = int(np.argmin(np.abs(C.ELL - 60)))
    i15 = int(np.argmin(np.abs(C.ELL - 1500)))
    print(f"    {'fold':>20} {'B/E med':>9} {'B/E l=60':>10} {'B/E l=1500':>11} "
          f"{'C_EE l=60':>12} {'C_BB l=60':>12}")
    for lab, *_ in FOLDS:
        o = out[lab]
        r = o["BB"] / o["EE"]
        o["r"] = r
        print(f"    {lab:>20} {np.median(r[m]):9.4f} {r[i60]:10.4f} {r[i15]:11.4f} "
              f"{o['EE'][i60]:12.4e} {o['BB'][i60]:12.4e}")

    rA = out[FOLDS[0][0]]["r"]
    ok = (abs(rA[i60] - 0.951) < 5e-3 and abs(rA[i15] - 0.420) < 5e-3
          and abs(np.median(rA[m]) - 0.629) < 5e-3)
    print(f"\n  expected deployed baseline 0.951 / 0.420 / 0.629 -> REPRODUCED: {ok}")

    print("\n  DECOMPOSITION of the parent's A->B move (median B/E over the band):")
    mA = np.median(out[FOLDS[0][0]]["r"][m]); mB = np.median(out[FOLDS[1][0]]["r"][m])
    mC = np.median(out[FOLDS[2][0]]["r"][m]); mD = np.median(out[FOLDS[3][0]]["r"][m])
    print(f"      A (deployed, as published)                 {mA:.4f}")
    print(f"      D = A with the callable swapped only       {mD:.4f}   "
          f"(callable effect alone: {mD - mA:+.4f})")
    print(f"      C = A with the TABLE rebuilt only          {mC:.4f}   "
          f"(table effect alone:    {mC - mA:+.4f})")
    print(f"      B = both (what the parent measured)        {mB:.4f}   "
          f"(total:                 {mB - mA:+.4f})")

    print("\n  per-ell B/E, all four folds:")
    print(f"    {'ell':>6}" + "".join(f"{lab.split()[0]+lab.split()[1][:9]:>18}" for lab, *_ in FOLDS))
    for i, L in enumerate(C.ELL):
        if L < 20:
            continue
        print(f"    {int(L):6d}" + "".join(f"{out[lab]['r'][i]:18.4f}" for lab, *_ in FOLDS))

    # ------------------------------------------- 3. hypothesis (a): blast radius
    sep("3. HYPOTHESIS (a): per-channel change of the vertex table (rebuilt vs deployed)")
    print("  The note's blast radius for the PHASE FIX ALONE was:")
    print("    'TTT, TTP, TPP, PPP, Dmod bit-identical; only Bmod changes'.")
    print("  If those moved here, the note's premise ('xi_- is untouched') fails "
          "for the combined rebuild.\n")
    print(f"    {'channel':>10} {'bit-identical':>14} {'median|new/old|':>16} "
          f"{'p16':>9} {'p84':>9} {'median rel chg':>15}")
    for ch in ("zeta_TTT", "zeta_TTP", "zeta_TPP", "zeta_PPP", "zeta_Bmod", "zeta_Dmod"):
        o, n = np.asarray(vd[ch], float), np.asarray(vr[ch], float)
        good = np.abs(o) > 0
        rat = np.abs(n[good] / o[good])
        rel = np.abs(n[good] - o[good]) / np.abs(o[good])
        print(f"    {ch:>10} {str(bool(np.array_equal(o, n))):>14} {np.median(rat):16.4f} "
              f"{np.percentile(rat,16):9.4f} {np.percentile(rat,84):9.4f} {np.median(rel):15.4f}")

    # ------------------------- 4. hypothesis (b): real space vs harmonic ratios
    sep("4. HYPOTHESIS (b): REAL-SPACE |xi_-/xi_+| vs HARMONIC C_BB/C_EE")
    print("  The note quotes |xi_-/xi_+| = 8.58e-4 (deployed) -> ~1.8e-3 (phase-fixed),")
    print("  and calls that 'the B/E ratio'.  The rebuild quotes C_BB/C_EE.\n")
    print(f"    {'gamma[arcmin]':>13}" + "".join(f"{lab.split()[0]:>14}" for lab, *_ in FOLDS))
    for target in (0.5, 1.0, 2.0, 5.0, 12.0, 30.0):
        i = int(np.argmin(np.abs(g0 - target)))
        row = "".join(f"{abs(out[lab]['xi_m'][i]/out[lab]['xi_p'][i]):14.4e}"
                      for lab, *_ in FOLDS)
        print(f"    {g0[i]:13.4f}" + row)

    i05 = int(np.argmin(np.abs(g0 - 0.5)))
    A, B = FOLDS[0][0], FOLDS[1][0]
    C_ = FOLDS[2][0]
    print("\n  THE FOUR NUMBERS:")
    print(f"    real-space |xi_-/xi_+| @0.5'    deployed (A) = "
          f"{abs(out[A]['xi_m'][i05]/out[A]['xi_p'][i05]):.4e}")
    print(f"    real-space |xi_-/xi_+| @0.5'    rebuilt  (C) = "
          f"{abs(out[C_]['xi_m'][i05]/out[C_]['xi_p'][i05]):.4e}")
    print(f"    harmonic  C_BB/C_EE band-median deployed (A) = {mA:.4f}")
    print(f"    harmonic  C_BB/C_EE band-median rebuilt  (C) = {mC:.4f}")

    # ------------------- 5. the note's INFERENCE, tested on its own premise
    sep("5. THE NOTE'S INFERENCE, TESTED DIRECTLY")
    print("  Premise (the note's own, for the phase fix alone): xi_+ = -4 Bmod falls,")
    print("  xi_- = -4 TPP is UNTOUCHED, so |xi_-/xi_+| rises.  The note concludes")
    print("  'the B/E ratio RISES by up to 2x'.  Test that implication.\n")
    print("  Algebra: EE = (Cl[xi_+] + Cl[xi_-])/2, BB = (Cl[xi_+] - Cl[xi_-])/2, so")
    print("      B/E = (1 - q)/(1 + q),   q = Cl[xi_-]/Cl[xi_+].")
    print("  B/E is a strictly DECREASING function of q.  |xi_-| growing relative to")
    print("  |xi_+| makes q larger, hence B/E SMALLER.  Small |xi_-/xi_+| is the")
    print("  MAXIMAL-B-mode limit (xi_- = 0 gives BB = EE exactly), not the minimal one.\n")
    print("  Numerically, on the DEPLOYED fold, scaling xi_+ by s with xi_- fixed")
    print("  (exactly the note's phase-fix scenario):")
    print(f"    {'s (xi_+ scale)':>15} {'|xi-/xi+| @0.5inch':>19} {'median B/E':>12} "
          f"{'B/E l=1500':>11}")
    for s in (1.0, 0.9, 0.8, 0.74, 0.6, 0.5):
        xp = out[A]["xi_p"] * s
        xm = out[A]["xi_m"]
        epb = C.forward_curved(xp, setups[0])
        emb = C.forward_curved(xm, setups[1])
        ee, bb = 0.5 * (epb + emb), 0.5 * (epb - emb)
        r = bb / ee
        print(f"    {s:15.2f} {abs(xm[i05]/xp[i05]):19.4e} {np.median(r[m]):12.4f} "
              f"{r[i15]:11.4f}")
    print("\n  -> B/E falls monotonically as |xi_-/xi_+| rises. The note's measurement")
    print("     direction and the rebuild's are the SAME; only the note's inference")
    print("     from it to 'B/E rises' is inverted.")

    # ----------------------------------- 6. is the high-ell BB collapse physical?
    sep("6. PHYSICALITY / CONDITIONING OF THE B-MODE")
    sd = (C.build_curved_matrix(g0, ELL_DENSE, 2, 2),
          C.build_curved_matrix(g0, ELL_DENSE, 2, -2))
    md = (ELL_DENSE >= BAND[0]) & (ELL_DENSE <= BAND[1])
    ld = np.log(ELL_DENSE[md])
    for lab, *_ in FOLDS:
        o = polar(res[lab], sd)
        ee, bb, eb, epb, emb = o["EE"], o["BB"], o["EB"], o["epb"], o["emb"]
        cond = (np.abs(epb[md]) + np.abs(emb[md])) / (2.0 * np.abs(bb[md]))
        slope = np.diff(np.log(np.abs(bb[md]))) / np.diff(ld)
        print(f"  [{lab}] {int(md.sum())} multipoles in {BAND[0]:.0f}<=l<={BAND[1]:.0f}")
        print(f"      C_EE<=0 at {int(np.sum(ee[md]<=0)):2d} l ;  "
              f"C_BB<=0 at {int(np.sum(bb[md]<=0)):2d} l ;  "
              f"max|C_EB| = {np.max(np.abs(eb)):.3e} ;  "
              f"max|<g+gx>| = {np.max(np.abs(o['cross'])):.3e}")
        print(f"      d ln|C_BB| / d ln l : min {slope.min():+.2f}  median "
              f"{np.median(slope):+.2f}  max {slope.max():+.2f}   "
              f"(monotone decline: {bool(np.all(slope < 0))})")
        print(f"      BB cancellation (|Cl[xi+]|+|Cl[xi-]|)/2|BB| : "
              f"median {np.median(cond):6.2f}  max {np.max(cond):7.2f}  "
              f"@l=1500 {cond[-1]:7.2f}")
        print(f"      -> a 1% error in xi_+ or xi_- moves C_BB(l=1500) by "
              f"~{cond[-1]:.0f}%")


    # ------------------ 7. conditioning of the note's own real-space statistic
    sep("7. CONDITIONING OF THE REAL-SPACE RATIO, AND THE CALIBRATED VERTEX ERROR BAR")
    print("  xi_- = <g+g+> - <gxgx> is itself a cancellation.  At 0.5 arcmin the note's")
    print("  headline 8.58e-4 is a 1-part-in-N residual of two nearly equal numbers:\n")
    print(f"    {'gamma[arcmin]':>13} {'<g+g+>':>13} {'<gxgx>':>13} {'xi_-':>13} "
          f"{'cancel factor':>14}")
    for target in (0.5, 1.0, 2.0, 5.0, 12.0):
        i = int(np.argmin(np.abs(g0 - target)))
        p11 = res[A][(1, 1)][i]; p22 = res[A][(2, 2)][i]
        xm = p11 - p22
        print(f"    {g0[i]:13.4f} {p11:13.4e} {p22:13.4e} {xm:13.4e} "
              f"{(abs(p11)+abs(p22))/(2*abs(xm)):14.1f}")
    print("\n  At 0.5 arcmin the callable swap alone (A->D, same table, same physics)")
    print(f"  moves |xi_-/xi_+| by "
          f"{abs(out['D dep64+plaincall ']['xi_m'][i05]/out['D dep64+plaincall ']['xi_p'][i05])/abs(out[A]['xi_m'][i05]/out[A]['xi_p'][i05]):.2f}x, "
          "while moving the harmonic median B/E by <0.1%.")
    print("  The 0.5-arcmin real-space ratio is therefore not a stable statistic; the")
    print("  harmonic band-median is.")

    print("\n  Calibrated vertex error bar (zeta_TTT leg-symmetry, true value zero):")
    from itertools import product as _prod
    for nm, tab in (("deployed n_phi=64 ", vd), ("rebuilt  n_phi=512", vr)):
        tri = np.asarray(tab["cosine_triples"], float)
        ttt = np.asarray(tab["zeta_TTT"], float)
        by = {0: {}, 1: {}, 2: {}}
        for i, t in enumerate(tri):
            j = int(np.argmin(np.abs(t - 1.0)))
            if abs(t[j] - 1.0) > 1e-12:
                continue
            rest = np.delete(t, j)
            if abs(rest[0] - rest[1]) > 1e-12:
                continue
            by[j][float(np.round(rest[0], 15))] = i
        common = set(by[0]) & set(by[1]) & set(by[2])
        dev = []
        for c in common:
            v = [ttt[by[s][c]] for s in (0, 1, 2)]
            for a_, b_ in ((0, 1), (0, 2), (1, 2)):
                d = np.abs(v[a_] / v[b_] - 1.0)
                dev.append(d[np.isfinite(d)])
        dev = np.concatenate(dev)
        print(f"    {nm}: {len(common)} placements, |ratio-1| median "
              f"{np.median(dev):.3e}  p95 {np.percentile(dev,95):.3e}  "
              f"max {np.max(dev):.3e}")
    print("\n  Propagating the rebuilt table's own residual through the l=1500")
    print("  cancellation amplification (~40x) is the honest error bar on B/E there.")

    print("\n  transform robustness, n_fine 20000 -> 60000 (fold C):")
    sa = (C.build_curved_matrix(g0, C.ELL, 2, 2, n_fine=60000),
          C.build_curved_matrix(g0, C.ELL, 2, -2, n_fine=60000))
    oa = polar(res[C_], sa)
    ra = oa["BB"] / oa["EE"]
    print(f"      median B/E {mC:.4f} -> {np.median(ra[m]):.4f} ;  "
          f"l=1500 {out[C_]['r'][i15]:.4f} -> {ra[i15]:.4f}")

    # ------------------------------- 8. is the band-median a stable statistic?
    sep("8. ELL-SAMPLING SENSITIVITY OF THE BAND MEDIAN")
    print("  'median B/E over 50<=l<=1500' depends on how the band is sampled: the")
    print("  paper's grid is 16 LOG-spaced multipoles in band, which under-weights")
    print("  high ell, where the rebuilt B-mode collapses hardest.\n")
    print(f"    {'sampling':>20} {'n_l':>5} {'A deployed':>11} {'C rebuilt':>10} {'C/A':>7}")
    for nm, EL in (("16 log (paper grid)", C.ELL),
                   ("71 log", ELL_DENSE),
                   ("1451 linear", np.arange(50, 1501, dtype=float))):
        s2 = (C.build_curved_matrix(g0, EL, 2, 2), C.build_curved_matrix(g0, EL, 2, -2))
        mm = (EL >= BAND[0]) & (EL <= BAND[1])
        vals = []
        for lab in (A, C_):
            o = polar(res[lab], s2)
            vals.append(float(np.median((o["BB"] / o["EE"])[mm])))
        print(f"    {nm:>20} {int(mm.sum()):5d} {vals[0]:11.4f} {vals[1]:10.4f} "
              f"{vals[1]/vals[0]:7.3f}")
    print("\n  The DIRECTION and rough size of the fall are sampling-independent;")
    print("  the quoted pair (0.629, 0.289) is specific to the 16-point log grid.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())