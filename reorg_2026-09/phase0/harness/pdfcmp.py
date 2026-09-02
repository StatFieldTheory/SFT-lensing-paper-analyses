#!/usr/bin/env python
"""Compare two PDFs: raw md5, metadata-date-stripped md5, and a raster diff.

Matplotlib PDFs differ only in /CreationDate when regenerated from identical
inputs (the date string has fixed length so xref offsets are unchanged), so the
stripped md5 is an exact content test. The raster diff catches everything else.
"""
import hashlib, json, re, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
from PIL import Image

DATE_RE = re.compile(rb"/(CreationDate|ModDate)\s*\((?:[^()\\]|\\.)*\)")

def md5(b): return hashlib.md5(b).hexdigest()
def stripped(b): return DATE_RE.sub(lambda m: b"/" + m.group(1) + b"(STRIPPED)", b)

def raster(pdf, dpi):
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "page"
        subprocess.run(["pdftoppm", "-r", str(dpi), "-png", "-singlefile", "-f", "1", "-l", "1",
                        str(pdf), str(out)], check=True, capture_output=True)
        return np.asarray(Image.open(str(out) + ".png").convert("RGB"), dtype=np.int16)

def compare(a, b, dpi=60):
    a, b = Path(a), Path(b)
    res = {"deployed": str(a), "candidate": str(b)}
    if not a.exists() or not b.exists():
        res["error"] = f"missing: {'deployed' if not a.exists() else 'candidate'}"
        return res
    ba, bb = a.read_bytes(), b.read_bytes()
    res.update(md5_deployed=md5(ba), md5_candidate=md5(bb),
               raw_identical=md5(ba) == md5(bb),
               stripped_identical=md5(stripped(ba)) == md5(stripped(bb)))
    ra, rb = raster(a, dpi), raster(b, dpi)
    if ra.shape != rb.shape:
        res.update(raster="SHAPE DIFFERS", shape_deployed=list(ra.shape), shape_candidate=list(rb.shape))
    else:
        d = np.abs(ra - rb).max(axis=2)
        res.update(raster_max_diff=int(d.max()), raster_frac_diff=round(float((d > 16).mean()), 6))
    return res

if __name__ == "__main__":
    r = compare(sys.argv[1], sys.argv[2])
    tag = sys.argv[3] if len(sys.argv) > 3 else ""
    verdict = ("IDENTICAL(raw)" if r.get("raw_identical") else
               "IDENTICAL(content)" if r.get("stripped_identical") else
               "DIFFERS" if "error" not in r else "ERROR")
    print(f"[{tag}] {verdict} raster_frac_diff={r.get('raster_frac_diff')} max={r.get('raster_max_diff')} {r.get('error','')}")
    print(json.dumps(r))
