"""Resolve the reviewed numerical products selected by active_products.json.

The manifest is required. The historical arrays that readers once fell back to
were archived on 2026-10-02 (reproduce/archived_inputs.json), so a missing,
partial or modified bundle fails instead of selecting older data.
"""

from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path, PureWindowsPath

HERE = Path(__file__).resolve().parent
ACTIVE = HERE / "active_products.json"


@lru_cache(maxsize=1)
def active_manifest() -> dict:
    if not ACTIVE.exists():
        raise FileNotFoundError(
            f"Numerical product manifest is missing: {ACTIVE}. There is no "
            "historical fallback: those inputs were archived on 2026-10-02. "
            "Restore the manifest from Git."
        )
    manifest = json.loads(ACTIVE.read_text())
    if manifest.get("schema_version") != 1 or not manifest.get("products"):
        raise ValueError(f"Invalid numerical product manifest: {ACTIVE}")
    return manifest


def _digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _contained_path(selector: str, *, base: Path, boundary: Path, label: str) -> Path:
    """Resolve a relative selector without allowing host-specific or escaping paths."""
    if not isinstance(selector, str) or not selector:
        raise ValueError(f"Invalid relative numerical product path: {label}")
    relative = Path(selector)
    if relative.is_absolute() or PureWindowsPath(selector).anchor:
        raise ValueError(f"Numerical product path must be relative: {label}")
    resolved = (base / relative).resolve()
    if not resolved.is_relative_to(boundary.resolve()):
        raise ValueError(f"Numerical product path leaves its recorded directory: {label}")
    return resolved


@lru_cache(maxsize=None)
def resolve_product(key: str) -> Path:
    """Return a package-local selected product after checking recorded hashes."""
    record = active_manifest()["products"][key]
    path = _contained_path(
        record["path"], base=HERE, boundary=HERE.parent, label=key
    )
    if "sha256" in record:
        if not path.is_file() or _digest(path) != record["sha256"]:
            raise ValueError(f"Numerical product changed or is missing: {key}: {path}")
    elif record.get("files"):
        if not path.is_dir():
            raise ValueError(f"Numerical product directory is missing: {key}: {path}")
        for name, expected in record["files"].items():
            child = _contained_path(name, base=path, boundary=path, label=name)
            if not child.is_file() or _digest(child) != expected:
                raise ValueError(f"Numerical product member changed or is missing: {child}")
    else:
        raise ValueError(f"Product has no recorded content hashes: {key}")
    return path


def cutoff_series() -> list[dict]:
    """Select explicit model/leg-order series, preserving the recorded legacy layout."""
    record = active_manifest()['products']['cutoff_ladder']
    declared = record.get('series')
    if declared is None:
        return [
            {'key': 'tree', 'model': 'tree', 'leg_order': None,
             'directory': '', 'label': 'tree-level bispectrum'},
            {'key': 'bihalofit', 'model': 'bihalofit', 'leg_order': None,
             'directory': '', 'label': 'BiHalofit (uncontrolled)'},
        ]
    expected = {('tree', 'auto'), ('bihalofit', 'auto'), ('bihalofit', 'given')}
    if not isinstance(declared, list) or len(declared) != len(expected):
        raise ValueError('Complete explicit tree and both BiHalofit leg-order series required')
    result, seen = [], set()
    for item in declared:
        if not isinstance(item, dict) or set(item) != {'model', 'leg_order', 'directory'}:
            raise ValueError('Exact cutoff series declaration required')
        pair = (item['model'], item['leg_order'])
        if pair not in expected or pair in seen or item['directory'] != '_'.join(pair):
            raise ValueError('Duplicate, unsupported or ambiguous cutoff series')
        seen.add(pair)
        key = '_'.join(pair)
        label = ('tree-level bispectrum' if pair[0] == 'tree'
                 else 'BiHalofit (' + pair[1] + ', uncontrolled)')
        result.append({**item, 'key': key, 'label': label})
    return result


def resolve_cutoff_member(series: dict, cutoff: int, *, root: Path | None = None) -> Path:
    """Require the exact selected series and a hash-declared finite-cutoff member."""
    if series not in cutoff_series() or type(cutoff) is not int or cutoff not in (960, 1920, 3840, 7680, 15360):
        raise ValueError('Exact selected finite-cutoff scope required')
    selected = resolve_product('cutoff_ladder')
    if root is not None and Path(root).resolve() != selected:
        raise ValueError('Cutoff renderer root must match the active manifest')
    member = Path(series['directory']) / f"table_{series['model']}_cut{cutoff}_r4_xi.npz"
    if member.as_posix() not in active_manifest()['products']['cutoff_ladder'].get('files', {}):
        raise ValueError('Cutoff member is missing from the hash-declared selection')
    return selected / member
