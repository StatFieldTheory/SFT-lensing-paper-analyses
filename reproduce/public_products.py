"""Export accepted NPZ products without private session metadata.

All retained NPY members are copied as bytes, including their headers. Only
explicitly classified scalar JSON metadata may be replaced. The returned
descriptor contains hashes and array structure, never source paths or metadata
contents; the accepted source archive remains unchanged.
"""

from __future__ import annotations

from collections.abc import Mapping
import hashlib
import io
import json
import os
from pathlib import Path
import re
import stat
import tempfile
from typing import Any
import zipfile

import numpy as np


_METADATA_MEMBERS = frozenset({
    "metadata_json.npy", "cl_table_metadata_json.npy", "cosmo_meta.npy",
})


class PublicProductError(ValueError):
    """An archive does not meet the accepted public-export contract."""


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check_source(source: Path, expected_sha256: str) -> None:
    if not isinstance(expected_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise PublicProductError("Expected source SHA-256 must be 64 lowercase hexadecimal characters")
    if _file_sha256(source) != expected_sha256:
        raise PublicProductError("Accepted source SHA-256 does not match")


def _archive_members(archive: zipfile.ZipFile) -> list[zipfile.ZipInfo]:
    members = archive.infolist()
    names = [member.filename for member in members]
    if len(names) != len(set(names)):
        raise PublicProductError("Duplicate NPZ member names")
    for member in members:
        name = member.filename
        unsafe = (
            name != member.orig_filename
            or "/" in name
            or "\\" in name
            or ":" in name
            or name.startswith(".")
            or not name.endswith(".npy")
            or name == ".npy"
            or stat.S_ISLNK(member.external_attr >> 16)
            or member.is_dir()
        )
        if unsafe:
            raise PublicProductError("Unsafe or unsupported NPZ member name")
    if not members:
        raise PublicProductError("An NPZ product must contain at least one array")
    return sorted(members, key=lambda member: member.filename)


def _array_descriptor(raw: bytes) -> dict[str, Any]:
    try:
        array = np.load(io.BytesIO(raw), allow_pickle=False)
    except (ValueError, OSError, EOFError) as error:
        raise PublicProductError("NPZ member is not a valid non-object NPY array") from error
    if not isinstance(array, np.ndarray) or array.dtype.hasobject:
        if isinstance(array, np.lib.npyio.NpzFile):
            array.close()
        raise PublicProductError("NPZ member is not a valid non-object NPY array")
    return {
        "sha256": hashlib.sha256(raw).hexdigest(),
        "nbytes": len(raw),
        "dtype": array.dtype.str,
        "shape": list(array.shape),
    }


def _json_npy(value: Mapping[str, Any], *, byte_vector: bool = False) -> bytes:
    if not isinstance(value, Mapping):
        raise PublicProductError("Metadata replacement must be a JSON object")
    try:
        encoded = json.dumps(dict(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as error:
        raise PublicProductError("Metadata replacement must be a finite JSON object") from error
    stream = io.BytesIO()
    array = (np.frombuffer(encoded.encode("utf-8"), dtype=np.uint8)
             if byte_vector else np.asarray(encoded))
    np.save(stream, array, allow_pickle=False)
    return stream.getvalue()


def _check_byte_vector_json(raw: bytes, description: dict[str, Any]) -> None:
    if len(description["shape"]) != 1 or description["dtype"] != "|u1":
        raise PublicProductError("Byte-vector metadata must be a one-dimensional uint8 array")
    array = np.load(io.BytesIO(raw), allow_pickle=False)
    try:
        value = json.loads(array.tobytes().decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as error:
        raise PublicProductError("Byte-vector metadata must contain a UTF-8 JSON object") from error
    if not isinstance(value, dict):
        raise PublicProductError("Byte-vector metadata must contain a UTF-8 JSON object")


def _expected_members(
    source: Path,
    expected_sha256: str,
    metadata_replacements: Mapping[str, Mapping[str, Any]] | None,
) -> tuple[dict[str, bytes], dict[str, Any], dict[str, Any]]:
    _check_source(source, expected_sha256)
    replacements = {} if metadata_replacements is None else dict(metadata_replacements)
    if set(replacements) - _METADATA_MEMBERS:
        raise PublicProductError("Replacement requested for an unclassified NPZ member")
    retained, preserved, replaced = {}, {}, {}
    with zipfile.ZipFile(source) as archive:
        members = _archive_members(archive)
        names = {member.filename for member in members}
        if set(replacements) - names:
            raise PublicProductError("Metadata replacement member is absent from the source")
        if "metadata_json.npy" in names:
            replacements.setdefault("metadata_json.npy", {
                "schema_version": 1,
                "original_sha256": expected_sha256,
                "export": "scientific NPY members preserved byte-for-byte",
            })
        for member in members:
            raw = archive.read(member)
            description = _array_descriptor(raw)
            if member.filename in replacements:
                byte_vector = member.filename == "cosmo_meta.npy"
                if byte_vector:
                    _check_byte_vector_json(raw, description)
                elif description["shape"] or description["dtype"][1] not in {"U", "S"}:
                    raise PublicProductError("Classified metadata must be a scalar string array")
                replacement = _json_npy(replacements[member.filename], byte_vector=byte_vector)
                retained[member.filename] = replacement
                replaced[member.filename] = {
                    "original_member_sha256": description["sha256"],
                    "public_member": _array_descriptor(replacement),
                }
            else:
                retained[member.filename] = raw
                preserved[member.filename] = description
    # A source being edited while it is read is not an accepted export input.
    _check_source(source, expected_sha256)
    return retained, preserved, replaced


def verify_npz(
    source: str | os.PathLike[str],
    public: str | os.PathLike[str],
    expected_sha256: str,
    *,
    metadata_replacements: Mapping[str, Mapping[str, Any]] | None = None,
    expected_public_sha256: str | None = None,
) -> dict[str, Any]:
    """Verify every public NPY member against the accepted source or replacement.

    Replacement keys are archive member names, including ``.npy``. Only
    ``metadata_json.npy`` and ``cl_table_metadata_json.npy`` use scalar strings;
    ``cosmo_meta.npy`` uses a one-dimensional uint8 UTF-8 JSON object. An existing
    ``metadata_json.npy`` receives compact provenance by default. Other metadata
    requires explicit replacement; ``config_json.npy`` is always retained.
    """
    source_path, public_path = Path(source), Path(public)
    expected, preserved, replaced = _expected_members(
        source_path, expected_sha256, metadata_replacements
    )
    public_sha256 = _file_sha256(public_path)
    if expected_public_sha256 is not None and public_sha256 != expected_public_sha256:
        raise PublicProductError("Public product SHA-256 does not match")
    with zipfile.ZipFile(public_path) as archive:
        members = _archive_members(archive)
        if {member.filename for member in members} != set(expected):
            raise PublicProductError("Public NPZ members differ from the accepted source")
        for member in members:
            if archive.read(member) != expected[member.filename]:
                raise PublicProductError("Public NPY member differs from its accepted bytes")
    if _file_sha256(public_path) != public_sha256:
        raise PublicProductError("Public product changed during verification")
    return {
        "schema_version": 1,
        "original_sha256": expected_sha256,
        "public_sha256": public_sha256,
        "preserved_members": preserved,
        "replaced_members": replaced,
    }


def export_npz(
    source: str | os.PathLike[str],
    target: str | os.PathLike[str],
    expected_sha256: str,
    *,
    metadata_replacements: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Create a deterministic public NPZ, refusing to overwrite other content.

    Repeating the export to an identical target is allowed. A different existing
    target, a symlink target, or any alias of the accepted source is rejected.
    """
    source_path, target_path = Path(source), Path(target)
    if source_path.resolve() == target_path.resolve() or (
        target_path.exists() and os.path.samefile(source_path, target_path)
    ):
        raise PublicProductError("The public target aliases the accepted source")
    if target_path.is_symlink():
        raise PublicProductError("A public target must not be a symlink")
    members, _, _ = _expected_members(source_path, expected_sha256, metadata_replacements)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=target_path.parent, prefix=".public-npz-", delete=False) as stream:
        temporary = Path(stream.name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for name, raw in members.items():
                member = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                member.compress_type = zipfile.ZIP_DEFLATED
                member.external_attr = 0o600 << 16
                archive.writestr(member, raw, compresslevel=9)
        descriptor = verify_npz(
            source_path, temporary, expected_sha256,
            metadata_replacements=metadata_replacements,
        )
        try:
            # Linking within one directory publishes only a complete archive and
            # atomically refuses an existing target rather than replacing it.
            os.link(temporary, target_path)
        except FileExistsError:
            if target_path.is_symlink() or _file_sha256(target_path) != descriptor["public_sha256"]:
                raise PublicProductError("Public target already contains different content") from None
        return descriptor
    finally:
        temporary.unlink(missing_ok=True)
