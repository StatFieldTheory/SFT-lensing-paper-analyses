#!/usr/bin/env python
"""Compare regenerated figures (each generator's outputs/) with the deployed
``figures/*.pdf`` of the manuscript.

    python reproduce/check_figures.py            # all 17
    python reproduce/check_figures.py --figure 2 3

Verdicts: IDENTICAL(raw) = same bytes; IDENTICAL(content) = same bytes apart from
the PDF creation date; EQUIVALENT = raster identical at 60 dpi but the tight bounding
box drifted by a fraction of a point (font-metric drift between environments);
DIFFERS = content differs.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pdfcmp import compare  # noqa: E402
from regen_figure import ALIASES, FIGURES, PAPER, PKG  # noqa: E402


def verdict(r: dict) -> str:
    if "error" in r:
        return "MISSING " + r["error"]
    if r.get("raw_identical"):
        return "IDENTICAL(raw)"
    if r.get("stripped_identical"):
        return "IDENTICAL(content)"
    if r.get("raster_frac_diff") == 0.0:
        return "EQUIVALENT(bbox drift)"
    return "DIFFERS"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figure", type=int, nargs="+")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args()
    todo = sorted({ALIASES.get(n, n) for n in a.figure}) if a.figure else sorted(FIGURES)
    rows, worst = [], 0
    for n in todo:
        _, _, _, _, outs, papers = FIGURES[n]
        for out, paper in zip(outs, papers):
            r = compare(PAPER / "figures" / paper, PKG / out)
            v = verdict(r)
            r.update(figure=n, paper=paper, verdict=v)
            rows.append(r)
            print(f"[fig {n:2d}] {paper:38s} {v:24s} raster diff frac {r.get('raster_frac_diff', '-')}")
            worst = max(worst, 0 if v.startswith("IDENTICAL") else 1 if v.startswith("EQUIVALENT") else 2)
    if a.json:
        a.json.write_text(json.dumps(rows, indent=1))
    return worst if worst == 2 else 0


if __name__ == "__main__":
    raise SystemExit(main())
