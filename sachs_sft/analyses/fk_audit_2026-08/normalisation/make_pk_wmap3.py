"""Generate the linear P(k, z=0) of the Sato et al. 2009 ray-tracing cosmology.

The published cosmic-shear convergence bispectrum of Takahashi et al. 2020
Figure 13 was measured on the Sato et al. 2009 / Kayo et al. 2013 maps, whose
cosmology is WMAP3-like: Omega_m = 0.238, Omega_b = 0.042, n_s = 0.958,
sigma8 = 0.76, h = 0.732 (Kayo et al. 2013, arXiv:1207.6322, Sec. 2).

Runs in the PyCCL conda environment (CAMB), writes a two-column table
``k [h/Mpc]  P [(Mpc/h)^3]`` that the canoes venv reads back.

  /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python \
      normalisation/make_pk_wmap3.py
"""

from __future__ import annotations

from pathlib import Path

import camb
import numpy as np

OMEGA_M = 0.238
OMEGA_B = 0.042
H = 0.732
N_S = 0.958
SIGMA8 = 0.76
OUT = Path(__file__).resolve().parent / "data" / "PCAMB_wmap3_sato2009_z0.txt"


def main() -> None:
    h2 = H * H
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=100.0 * H,
        ombh2=OMEGA_B * h2,
        omch2=(OMEGA_M - OMEGA_B) * h2,
        mnu=0.0,
        omk=0.0,
        tau=0.06,
    )
    pars.InitPower.set_params(ns=N_S, As=2.1e-9)
    pars.set_matter_power(redshifts=[0.0], kmax=600.0)
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)
    s8_raw = float(results.get_sigma8_0())

    k, _z, pk = results.get_matter_power_spectrum(
        minkh=1e-5, maxkh=500.0, npoints=1200
    )
    p = pk[0] * (SIGMA8 / s8_raw) ** 2

    OUT.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "linear matter P(k, z=0), Sato et al. 2009 / Kayo et al. 2013 cosmology\n"
        f"Omega_m={OMEGA_M} Omega_b={OMEGA_B} h={H} n_s={N_S} "
        f"sigma8={SIGMA8} (CAMB raw sigma8={s8_raw:.6f}, rescaled)\n"
        "k [h/Mpc]    P [(Mpc/h)^3]"
    )
    np.savetxt(OUT, np.column_stack([k, p]), header=header)
    print(f"wrote {OUT}  ({k.size} rows, k = {k[0]:.3e} .. {k[-1]:.3e} h/Mpc)")
    print(f"CAMB sigma8 = {s8_raw:.6f} -> rescaled to {SIGMA8}")


if __name__ == "__main__":
    main()
