# Superseded as the production vertex (2026-09-02)

`equal_time_limber_kappa3_callable.py` and `equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz`
are the cut1000, corner-frozen vertex of the June 2026 draft. They are kept byte-identical because
SFT-WL-B vendored them (2026-08-31) and because `build_equal_time_limber_table.py`,
`ell_band_decomp.py` and `plot_zeta_figures.py` here are still the live build and figure-7 code.

The production vertex of the manuscript is `../equal_time_limber_cut15360_permaware/`
(permutation-aware lookup, ell_max = 15360). The callable here sorts the query triple and rounds
the mesh key to 8 decimals, which is wrong for the spin channels (23% on zeta_TPP, 70% on
zeta_Dmod at the collapsed configurations) and cannot represent separations below 0.344 arcmin;
see `../equal_time_limber_cut15360_permaware/README.md`.

The "+3.086e-5" baseline quoted in this folder's README and docstrings is a pre-2026-06-09 stamp;
see `sftwick_outputs/2PCF/C_corr_op_K_limber_FK/README.md` for the history of that number.
