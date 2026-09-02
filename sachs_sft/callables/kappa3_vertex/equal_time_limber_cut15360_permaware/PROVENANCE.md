# PROVENANCE

* `perm_aware_kappa3_callable.py`, `compare_callables.py`, `README.md`, `tests/test_perm_aware.py`
  came from `driver_field_emulators/code/callable_fixed/` (2026-08-26, moved here 2026-09-02 with
  history). The callable and the README are byte-identical to the copies SFT-WL-B vendored on
  2026-08-31; `tests/test_perm_aware.py` and `compare_callables.py` had their path arithmetic
  updated for the new location (old contents at tag `dfe-imported-2026-09`).
* `table_permclosed.npz` is `driver_field_emulators/products/table_permclosed_cut15360.npz`,
  byte-identical (sha256 prefix `b400dc5b11e1033b`), renamed to the name the callable resolves.
  Built 2026-08-26 09:46 by `rebuild/assemble.py` from `rebuild/products/pieces/low_permclosed_r4.npz`
  plus the `band_tree_permclosed_*` and `band_tree_extra_*` octave windows (tree-level bispectrum,
  ell_max = 15360, permutation-closed r4 grid plus extra collapsed rows); logs in `rebuild/logs/`.
* This folder is the production vertex for every FK result in the manuscript since the 2026-08
  revision; the superseded cut1000 folder is `../equal_time_limber/` (see its DEPRECATED.md).
