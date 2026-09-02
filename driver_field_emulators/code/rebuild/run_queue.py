"""Run build jobs with a hard cap on concurrency and free memory.

Written after an out-of-memory reboot. The HIGH branch holds a working
array proportional to (rows x n_ell x n_ell x n_phi), so an unchunked build
on the full triple set is several GB, and eighteen of them at once exhausted
a 103 GB machine. Two guards, because either alone is insufficient:

* each job is chunked over rows inside the build script, which bounds the
  per-process footprint
* this runner starts a job only when measured free memory exceeds a floor,
  and never runs more than `--max-jobs` at once

The free-memory floor matters because this machine may be running other
work. The runner re-checks before every launch rather than assuming the
footprint it measured at the start still holds.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path


def free_gb() -> float:
    """Free plus inactive pages, in GB. Inactive pages are reclaimable."""
    out = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    page = 16384
    values = {}
    for line in out.splitlines():
        if ":" not in line:
            continue
        key, _, value = line.partition(":")
        digits = value.strip().rstrip(".")
        if digits.isdigit():
            values[key.strip()] = int(digits)
    free = values.get("Pages free", 0) + values.get("Pages inactive", 0)
    return free * page / 1e9


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs-json", type=Path, required=True,
                    help="JSON list of {name, argv, log}")
    ap.add_argument("--max-jobs", type=int, default=3)
    ap.add_argument("--min-free-gb", type=float, default=25.0)
    ap.add_argument("--poll-seconds", type=float, default=10.0)
    args = ap.parse_args()

    jobs = json.loads(args.jobs_json.read_text())
    pending = [j for j in jobs if not Path(j["out"]).exists()]
    skipped = len(jobs) - len(pending)
    if skipped:
        print(f"[queue] {skipped} job(s) already have output, skipping")
    running: list[tuple[dict, subprocess.Popen]] = []
    done, failed = [], []

    while pending or running:
        for job, process in list(running):
            if process.poll() is None:
                continue
            running.remove((job, process))
            (done if process.returncode == 0 else failed).append(job["name"])
            state = "ok" if process.returncode == 0 else f"FAILED {process.returncode}"
            print(f"[queue] {job['name']}: {state} "
                  f"({len(done)} done, {len(failed)} failed, "
                  f"{len(pending)} queued)", flush=True)

        while pending and len(running) < args.max_jobs:
            available = free_gb()
            if available < args.min_free_gb and running:
                break
            job = pending.pop(0)
            log = open(job["log"], "w")
            process = subprocess.Popen(job["argv"], stdout=log, stderr=log)
            running.append((job, process))
            print(f"[queue] start {job['name']} "
                  f"(free {available:.0f} GB, {len(running)} running)", flush=True)
            time.sleep(2.0)

        if running:
            time.sleep(args.poll_seconds)

    print(f"[queue] finished: {len(done)} ok, {len(failed)} failed")
    if failed:
        print("[queue] failed: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
