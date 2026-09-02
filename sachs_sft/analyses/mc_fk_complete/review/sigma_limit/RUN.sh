#!/bin/sh
# Reproduce every number in review/sigma_limit.  ~12 min wall, 2 processes max.
set -e
cd "$(dirname "$0")/../.."
export PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
export OMP_NUM_THREADS=2
S=review/sigma_limit

# --- symbolic: the O(sigma) coefficients, exactly -------------------------
wolframscript -file $S/kink_moments.wl
wolframscript -file $S/kink_toy_numeric.wl     # kink -> O(sigma); smooth -> 0

# --- the load-bearing Wick identity, re-derived from SPEC alone -----------
$PY $S/identity_check.py

# --- the requested fixed sigma/dlam = 4.19 scans --------------------------
$PY $S/scan_sigma.py 1.0 "251:32,500:16,1000:8,1999:4,3997:2" $S/_fixratio_g1.npz
$PY $S/scan_sigma.py 5.0 "251:32,500:16,1000:8,1999:4,3997:2" $S/_fixratio_g5.npz

# --- controls -------------------------------------------------------------
# (a) fixed sigma, varying N: isolates the lattice contribution
$PY $S/scan_sigma.py 1.0 "251:32,500:32,1000:32,500:8,1000:8,2000:8,4000:8" $S/_dlctl_g1.npz
$PY $S/scan_sigma.py 5.0 "251:32,500:32,1000:32,500:8,1000:8,2000:8,4000:8" $S/_dlctl_g5.npz
# (b) production grid N=1000 + a second fixed-ratio series at sigma/dlam=8.39
$PY $S/scan_sigma.py 1.0 "1000:16,1000:8,1000:4,1000:2,2000:8,4000:4,8000:2" $S/_prod_g1.npz
$PY $S/scan_sigma.py 5.0 "1000:16,1000:8,1000:4,1000:2,2000:8,4000:4,8000:2" $S/_prod_g5.npz
# (c) switch the regulator OFF on a fixed lattice: is anything left at sigma=0?
$PY $S/scan_sigma.py 1.0 "400:20,400:10,400:5,400:2,400:1,400:0.5,400:0.25,400:0.125" $S/_white_g1.npz
# (d) the PUBLISHED (nominal) calibration, same scan: linearity is NOT generic
$PY $S/scan_sigma.py 1.0 "251:32,500:16,1000:8,1999:4,3997:2" $S/_nominal_g1.npz nominal

# --- premises and sigma-independent offsets -------------------------------
NL=2000 $PY $S/calibration_limit.py
$PY $S/reference_offsets.py

# --- the verdict ----------------------------------------------------------
$PY $S/verdict.py | tee $S/VERDICT.log
