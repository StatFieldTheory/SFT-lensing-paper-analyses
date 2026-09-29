# Polarization interpretation required by the FK correction

## Scope

This is a manuscript revision boundary, not a new higher-order calculation.
The current scalar, unperturbed-path example uses the linearized map from
Sachs fields to observable shear. Its linear shear has no B component.
The FK topology correlates one linear shear with one shear containing a
single nonlinear Sachs insertion. Consequently FK cannot contribute BB.
Statistical parity additionally removes EB. FF can contribute BB at the
same vertex order.

The original interpretation of FK as comparable E and B power must therefore
be removed, even if the corrected finite-grid transform retains a nonzero
BB residual. A residual is not a physical prediction without a numerical
leakage assessment. No higher-order non-Gaussian BB spectrum is computed in
this revision.

## Existing evidence inspected

- The manuscript itself states that the scalar linear reference has no B-mode,
  in `sections/conclusion.tex` immediately before the polarization discussion.
- `sections/path_int.tex` gives the FK topology and propagators.
- `/Users/zzhang/projects/SFT-WL-B/analysis/wpa/fk_pure_e.py` explains the
  linear external leg and tests the spin-leg placement of the corrected table.
- `/Users/zzhang/projects/SFT-WL-B/docs/leading_bmode_postfix.md` distinguishes
  the BB selection rule from uncomputed higher-order spectra.
- `/Users/zzhang/projects/SFT-WL-B/docs/plan_review_postfix.md` identifies the
  old FK B/E claim as a leg-placement artifact and leaves numerical residual
  and independent-input limitations open.

The earlier project's fold used a different pinned engine and historical
source endpoint. Its numerical values are not adopted here. The new figures
must use this revision's corrected inputs, current engine and actual z_s=5.

## Minimal manuscript consequence

Replace the paragraph interpreting small xi_minus as comparable E and B
power with the linear-external-leg argument. Retain the FF discussion, with
its updated numerical ratio and stated approximation. Label any displayed
FK BB as a numerical residual and do not quote it as signal. The abstract
is explicitly excluded from editing by the author, so any conflict with
its existing polarization statement must be reported separately.
