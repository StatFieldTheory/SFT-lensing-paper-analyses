"""Emit the ℓ-band decomposition of ζ_TTT as a booktabs LaTeX table.

Replaces the (hard-to-read) ℓ-band figure: at the source shell z=5.7, tabulates
the signed ζ_TTT contribution from three ℓ groups (ℓ≤60 exact Wigner-3j,
ℓ∈(60,250] and ℓ∈(250,1000] flat-sky Born-Limber) vs angular scale γ, with the
total and the high-ℓ (ℓ>60) share. Writes a self-contained `table` environment.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_NPZ = _HERE / "ell_band_decomp_results.npz"
_OUT = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
            "_draft_ell_band_table.tex")


def _sci(x: float) -> str:
    """LaTeX scientific notation, e.g. -2.45\\times10^{-11}."""
    if x == 0:
        return r"$0$"
    s = f"{x:.2e}"
    m, e = s.split("e")
    return rf"${m}\times10^{{{int(e)}}}$"


def main() -> None:
    d = np.load(_NPZ, allow_pickle=False)
    g = d["gamma_arcmin"]; z = d["z_shells"]
    iz = len(z) - 1
    B = np.stack([d[f"band{j}_TTT"] for j in range(6)], 0)[:, :, iz]  # (6, ng)
    full = d["full_TTT"][:, iz]
    g_lo = B[0] + B[1]                  # ℓ ≤ 60
    g_mid = B[2] + B[3]                 # ℓ ∈ (60, 250]
    g_hi = B[4] + B[5]                  # ℓ ∈ (250, 1000]
    sh = g_mid + g_hi                   # all high-ℓ
    frac = np.abs(sh) / (np.abs(sh) + np.abs(g_lo)) * 100.0

    sel_gamma = [0.5, 5.45, 21.31, 42.16, 117.31, 326.43]
    rows = [int(np.argmin(np.abs(g - gv))) for gv in sel_gamma]

    lines = [
        r"\begin{table}[t]",
        r"  \centering",
        r"  \caption{$\ell$-band decomposition of the squeezed driving-field "
        r"cumulant $\zeta_{TTT}(\gamma,\lambda)$ at the source shell "
        rf"$z={z[iz]:.1f}$. The exact Wigner-3j "
        r"branch ($\ell\le 60$) is $\gamma$-flat and positive; the flat-sky "
        r"Born--Limber branch ($\ell>60$) is negative and decays with $\gamma$. "
        r"At arcminute scales the high-$\ell$ (small-scale) modes dominate "
        r"$\zeta_{TTT}$ by $\sim\!2$ orders of magnitude and set its sign; the "
        r"low-$\ell$ modes take over only beyond the crossover near "
        r"$\gamma\approx40'$. ``High-$\ell$ share'' is "
        r"$|\zeta^{\ell>60}|/(|\zeta^{\ell>60}|+|\zeta^{\ell\le60}|)$.}",
        r"  \label{tab: zeta ell bands}",
        r"  \begin{tabular}{r r r r r c}",
        r"    \toprule",
        r"    $\gamma$ [arcmin] & $\ell\le 60$ & $\ell\in(60,250]$ & "
        r"$\ell\in(250,1000]$ & total & high-$\ell$ share \\",
        r"    \midrule",
    ]
    for i in rows:
        lines.append(
            f"    {g[i]:.2f} & {_sci(g_lo[i])} & {_sci(g_mid[i])} & "
            f"{_sci(g_hi[i])} & {_sci(full[i])} & {frac[i]:.0f}\\% \\\\"
        )
    lines += [
        r"    \bottomrule",
        r"  \end{tabular}",
        r"\end{table}",
        "",
    ]
    _OUT.write_text("\n".join(lines))
    print(f"[table] wrote {_OUT}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
