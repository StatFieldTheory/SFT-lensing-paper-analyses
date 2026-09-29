"""Resolve reviewed numerical products without replacing historical arrays.

An explicit active_products.json selects a complete revision bundle. With no
active bundle, historical generators retain their recorded inputs. A partial
or modified active bundle fails rather than falling back to historical data.
"""

from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ACTIVE = HERE / "active_products.json"


@lru_cache(maxsize=1)
def active_manifest() -> dict | None:
    if not ACTIVE.exists():
        return None
    manifest = json.loads(ACTIVE.read_text())
    if manifest.get("schema_version") != 1 or not manifest.get("products"):
        raise ValueError(f"Invalid numerical product manifest: {ACTIVE}")
    return manifest


def _digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


@lru_cache(maxsize=None)
def resolve_product(key: str, historical: Path) -> Path:
    """Return the selected file or directory after checking recorded hashes."""
    manifest = active_manifest()
    if manifest is None:
        return historical
    record = manifest["products"][key]
    path = Path(record["path"])
    if not path.is_absolute():
        path = HERE / path
    path = path.resolve()
    if "sha256" in record:
        if not path.is_file() or _digest(path) != record["sha256"]:
            raise ValueError(f"Numerical product changed or is missing: {key}: {path}")
    elif record.get("files"):
        if not path.is_dir():
            raise ValueError(f"Numerical product directory is missing: {key}: {path}")
        for name, expected in record["files"].items():
            child = (path / name).resolve()
            if not child.is_relative_to(path):
                raise ValueError(f"Product member leaves its recorded directory: {name}")
            if not child.is_file() or _digest(child) != expected:
                raise ValueError(f"Numerical product member changed or is missing: {child}")
    else:
        raise ValueError(f"Product has no recorded content hashes: {key}")
    return path
