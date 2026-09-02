"""Shared wiring for the circularity review.  Only reads the production tree."""
import os, sys
from pathlib import Path
os.environ.setdefault("OMP_NUM_THREADS", "2")
_MC = Path(__file__).resolve().parents[2]          # .../mc_fk_complete
if str(_MC) not in sys.path:
    sys.path.insert(0, str(_MC))
CANOES = "/Users/zzhang/projects/angular_statistics/canoes/src"
if CANOES not in sys.path:
    sys.path.insert(0, CANOES)
