"""Locate superseded inputs that were moved out of the working tree.

The files listed in archived_inputs.json are superseded arrays and tables. The
referee-revision run records pin many of them by path and SHA-256, so the
recomputation drivers still need their bytes. ``locate`` returns a file at its
recorded path when it is present there, and otherwise its archived copy. Either
way a file that the index lists is checked against its recorded hash. A record
therefore keeps its original path, and a guard that compares a record with the
bytes on disk keeps its strength.

Nothing in the figure or number readers uses this module. Those take their
inputs from active_products.json alone.
"""

from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
INDEX = HERE / "archived_inputs.json"


@lru_cache(maxsize=1)
def archive_index() -> dict:
    """Return the tracked index, or an empty one when nothing was archived."""
    if not INDEX.exists():
        return {"files": {}}
    index = json.loads(INDEX.read_text())
    if index.get("schema_version") != 1 or "files" not in index:
        raise ValueError(f"Invalid archived-input index: {INDEX}")
    return index


def _digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _index_key(path: Path) -> str | None:
    """Package-relative name of a recorded path, without following links.

    A link placed at an archived path keeps that path's key, so the bytes
    behind it are still checked against the index.
    """
    for candidate in (Path(os.path.abspath(path)), Path(path).resolve()):
        try:
            return candidate.relative_to(PACKAGE).as_posix()
        except ValueError:
            continue
    return None


def _check(stored: Path, key: str) -> Path:
    record = archive_index()["files"][key]
    if stored.stat().st_size != record["bytes"] or _digest(stored) != record["sha256"]:
        raise ValueError(f"Archived input differs from its recorded hash: {stored}")
    return stored


@lru_cache(maxsize=None)
def _archived_copy(key: str) -> Path:
    index = archive_index()
    stored = (HERE / index["archive_root"] / key).resolve()
    if not stored.is_file():
        raise FileNotFoundError(
            f"Archived input is not on this machine: {key}. Recover it "
            f"from commit {index['recover_from_commit']} of this repository, "
            f"or restore {index['archive_root']} relative to {HERE}."
        )
    return _check(stored, key)


def is_archived(path: Path) -> bool:
    """Whether the index lists this recorded path."""
    return _index_key(path) in archive_index()["files"]


def locate(path: Path) -> Path:
    """Return the readable location of a recorded input.

    A file the index does not list is returned when it exists. A listed file
    is returned from its recorded path when it was restored there, otherwise
    from the archive, and in both cases only after its size and SHA-256 match
    the index, so different bytes written under an archived path are refused.
    Anything else is an error.
    """
    recorded = Path(os.path.abspath(path))
    key = _index_key(recorded)
    listed = key in archive_index()["files"]
    if recorded.is_file():
        return _check(recorded.resolve(), key) if listed else recorded.resolve()
    if not listed:
        raise FileNotFoundError(f"Required input is missing: {recorded}")
    return _archived_copy(key)
