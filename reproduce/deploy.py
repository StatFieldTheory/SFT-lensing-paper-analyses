#!/usr/bin/env python
"""Copy regenerated figures into the manuscript's ``figures/`` under the paper's
file names, recording md5 sums in ``reproduce/deployed_md5.txt``.

    python reproduce/deploy.py --figure 8 10 --dry-run
    python reproduce/deploy.py --all

This is the only step that writes into the paper repository. It refuses to copy
a file whose generator output is missing, and prints old and new md5 so the
change is visible in the git history of the paper repository.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from regen_figure import ALIASES, FIGURES, PAPER, PKG  # noqa: E402


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest() if p.exists() else "-"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figure", type=int, nargs="+")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    todo = sorted(FIGURES) if a.all else sorted({ALIASES.get(n, n) for n in (a.figure or [])})
    if not todo:
        ap.error("give --figure N ... or --all")
    record = PKG / "reproduce" / "deployed_md5.txt"
    lines = {}
    if record.exists():
        for line in record.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                h, name = line.split(None, 1)
                lines[name.strip()] = h
    for n in todo:
        _, _, _, _, outs, papers = FIGURES[n]
        for out, paper in zip(outs, papers):
            src, dst = PKG / out, PAPER / "figures" / paper
            if not src.exists():
                print(f"[fig {n:2d}] MISSING output {out}; run regen_figure.py first")
                continue
            print(f"[fig {n:2d}] {paper}: {md5(dst)[:12]} -> {md5(src)[:12]}{'  (dry run)' if a.dry_run else ''}")
            if not a.dry_run:
                shutil.copy2(src, dst)
                lines[paper] = md5(dst)
    if not a.dry_run:
        body = "# md5 of the deployed figures/*.pdf, written by reproduce/deploy.py\n"
        body += "".join(f"{h}  {name}\n" for name, h in sorted(lines.items()))
        record.write_text(body)
        print(f"recorded {len(lines)} md5 sums in {record.relative_to(PKG)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
