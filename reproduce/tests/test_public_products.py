"""Public exports preserve scientific bytes while excluding private metadata."""

from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import sys
import warnings
import zipfile

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from public_products import PublicProductError, export_npz, verify_npz


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def npy(array: np.ndarray) -> bytes:
    stream = io.BytesIO()
    np.save(stream, array)
    return stream.getvalue()


def write_archive(path: Path, members: list[tuple[str, bytes]]) -> str:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        with zipfile.ZipFile(path, "w") as archive:
            for name, raw in members:
                archive.writestr(name, raw)
    return digest(path)


@pytest.fixture
def accepted(tmp_path: Path) -> tuple[Path, str]:
    source = tmp_path / "accepted.npz"
    # Two distinct NaN payloads and signed zero are deliberately represented as
    # raw IEEE-754 words; a numerical equality test would lose this distinction.
    words = np.array([0x7FF8000000000001, 0x7FF8000000000002, 0x8000000000000000, 0], dtype="<u8")
    private = {"session": "private conversation", "path": "/private/author/run"}
    config = {"grid": [1, 2, 3], "retained": "loader configuration"}
    expected = write_archive(source, [
        ("xi.npy", npy(words.view("<f8"))),
        ("grid.npy", npy(np.array([[1, 2], [3, 4]], dtype=">i4", order="F"))),
        ("metadata_json.npy", npy(np.asarray(json.dumps(private)))),
        ("cl_table_metadata_json.npy", npy(np.asarray(json.dumps({"convention": "spin-2"})))),
        ("config_json.npy", npy(np.asarray(json.dumps(config)))),
    ])
    return source, expected


def test_raw_members_and_required_config_are_preserved(accepted, tmp_path):
    source, expected = accepted
    before = source.read_bytes()
    target = tmp_path / "public.npz"
    descriptor = export_npz(source, target, expected)
    with zipfile.ZipFile(source) as original, zipfile.ZipFile(target) as public:
        assert set(original.namelist()) == set(public.namelist())
        for name in original.namelist():
            if name != "metadata_json.npy":
                assert public.read(name) == original.read(name)
        metadata = str(np.load(io.BytesIO(public.read("metadata_json.npy")), allow_pickle=False))
    assert "private conversation" not in metadata
    assert "/private/author/run" not in metadata
    assert expected in metadata
    assert source.read_bytes() == before
    assert descriptor["preserved_members"]["xi.npy"]["shape"] == [4]
    assert descriptor["preserved_members"]["grid.npy"]["dtype"] == ">i4"
    assert set(descriptor["replaced_members"]) == {"metadata_json.npy"}
    assert str(tmp_path) not in json.dumps(descriptor)
    assert descriptor == verify_npz(source, target, expected, expected_public_sha256=digest(target))
    assert export_npz(source, target, expected) == descriptor


def test_explicit_table_metadata_replacement(accepted, tmp_path):
    source, expected = accepted
    replacement = {"convention": "spin-2", "rebuild_manifest": {"pk_path": "power.npz", "sha256": "a" * 64}}
    choices = {"cl_table_metadata_json.npy": replacement}
    target = tmp_path / "public.npz"
    descriptor = export_npz(source, target, expected, metadata_replacements=choices)
    assert set(descriptor["replaced_members"]) == {"metadata_json.npy", "cl_table_metadata_json.npy"}
    with zipfile.ZipFile(source) as original, zipfile.ZipFile(target) as public:
        assert public.read("config_json.npy") == original.read("config_json.npy")
        assert public.read("xi.npy") == original.read("xi.npy")
        text = np.load(io.BytesIO(public.read("cl_table_metadata_json.npy")), allow_pickle=False).item()
        assert json.loads(text) == replacement
    assert descriptor == verify_npz(source, target, expected, metadata_replacements=choices)


def test_byte_vector_metadata_requires_explicit_replacement(tmp_path):
    source = tmp_path / "vertex.npz"
    private = json.dumps({"path": "/private/author/table.npz", "ell_max": 15360}).encode("utf-8")
    expected = write_archive(source, [
        ("cosmo_meta.npy", npy(np.frombuffer(private, dtype=np.uint8))),
        ("vertex.npy", npy(np.arange(3, dtype=np.float64))),
    ])
    retained_target = tmp_path / "retained.npz"
    retained = export_npz(source, retained_target, expected)
    assert not retained["replaced_members"]
    with zipfile.ZipFile(source) as original, zipfile.ZipFile(retained_target) as public:
        assert original.read("cosmo_meta.npy") == public.read("cosmo_meta.npy")
    target = tmp_path / "public.npz"
    choices = {"cosmo_meta.npy": {"path": "table.npz", "ell_max": 15360}}
    descriptor = export_npz(source, target, expected, metadata_replacements=choices)
    assert set(descriptor["replaced_members"]) == {"cosmo_meta.npy"}
    with zipfile.ZipFile(source) as original, zipfile.ZipFile(target) as public:
        assert original.read("vertex.npy") == public.read("vertex.npy")
        array = np.load(io.BytesIO(public.read("cosmo_meta.npy")), allow_pickle=False)
        assert array.ndim == 1 and array.dtype == np.uint8
        assert json.loads(array.tobytes().decode("utf-8")) == choices["cosmo_meta.npy"]
    assert descriptor == verify_npz(source, target, expected, metadata_replacements=choices)


@pytest.mark.parametrize("array", [
    np.array([1, 2], dtype=np.int32),
    np.frombuffer(b"invalid JSON", dtype=np.uint8),
    np.frombuffer(b"[1, 2]", dtype=np.uint8),
    np.frombuffer(b"\xff", dtype=np.uint8),
])
def test_invalid_byte_vector_metadata_is_rejected(tmp_path, array):
    source = tmp_path / "vertex.npz"
    expected = write_archive(source, [("cosmo_meta.npy", npy(array))])
    with pytest.raises(PublicProductError, match="Byte-vector metadata"):
        export_npz(source, tmp_path / "public.npz", expected, metadata_replacements={"cosmo_meta.npy": {}})


def test_wrong_source_hash_does_not_create_target(accepted, tmp_path):
    source, _ = accepted
    target = tmp_path / "not-created" / "public.npz"
    with pytest.raises(PublicProductError, match="source SHA-256"):
        export_npz(source, target, "0" * 64)
    assert not target.parent.exists()


@pytest.mark.parametrize("name", ["../xi.npy", "/xi.npy", "nested/xi.npy", "nested\\xi.npy", "C:xi.npy", "xi.txt", ".npy"])
def test_unsafe_members_are_rejected(tmp_path, name):
    source = tmp_path / "source.npz"
    expected = write_archive(source, [(name, npy(np.arange(2)))])
    with pytest.raises(PublicProductError, match="Unsafe or unsupported"):
        export_npz(source, tmp_path / "public.npz", expected)


def test_duplicate_members_are_rejected(tmp_path):
    source = tmp_path / "source.npz"
    expected = write_archive(source, [("xi.npy", npy(np.arange(2)))] * 2)
    with pytest.raises(PublicProductError, match="Duplicate"):
        export_npz(source, tmp_path / "public.npz", expected)


@pytest.mark.parametrize("name", ["xi.npy", "metadata_json.npy"])
def test_object_arrays_are_rejected(tmp_path, name):
    source = tmp_path / "source.npz"
    expected = write_archive(source, [(name, npy(np.array([{"private": "session"}], dtype=object)))])
    with pytest.raises(PublicProductError, match="non-object"):
        export_npz(source, tmp_path / "public.npz", expected)


@pytest.mark.parametrize("alias", ["same", "symlink", "hardlink"])
def test_source_alias_target_is_rejected(accepted, tmp_path, alias):
    source, expected = accepted
    target = source if alias == "same" else tmp_path / "alias.npz"
    if alias == "symlink":
        target.symlink_to(source)
    elif alias == "hardlink":
        target.hardlink_to(source)
    before = source.read_bytes()
    with pytest.raises(PublicProductError, match="aliases"):
        export_npz(source, target, expected)
    assert source.read_bytes() == before


def test_other_existing_target_is_not_overwritten(accepted, tmp_path):
    source, expected = accepted
    target = tmp_path / "public.npz"
    target.write_bytes(b"preserve this other artifact")
    with pytest.raises(PublicProductError, match="different content"):
        export_npz(source, target, expected)
    assert target.read_bytes() == b"preserve this other artifact"


@pytest.mark.parametrize("name", ["xi.npy", "config_json.npy", "unknown.npy"])
def test_unclassified_replacement_is_rejected(accepted, tmp_path, name):
    source, expected = accepted
    with pytest.raises(PublicProductError, match="unclassified"):
        export_npz(source, tmp_path / "public.npz", expected, metadata_replacements={name: {}})


def test_classified_numeric_member_cannot_be_replaced(tmp_path):
    source = tmp_path / "source.npz"
    expected = write_archive(source, [("metadata_json.npy", npy(np.arange(2)))])
    with pytest.raises(PublicProductError, match="scalar string"):
        export_npz(source, tmp_path / "public.npz", expected)


@pytest.mark.parametrize("mutation", ["modified", "missing", "added"])
def test_public_member_tampering_is_rejected(accepted, tmp_path, mutation):
    source, expected = accepted
    target = tmp_path / "public.npz"
    export_npz(source, target, expected)
    with zipfile.ZipFile(target) as archive:
        members = [(name, archive.read(name)) for name in archive.namelist()]
    if mutation == "modified":
        members = [(name, npy(np.zeros(4)) if name == "xi.npy" else raw) for name, raw in members]
    elif mutation == "missing":
        members = [(name, raw) for name, raw in members if name != "xi.npy"]
    else:
        members.append(("extra.npy", npy(np.zeros(1))))
    tampered = tmp_path / "tampered.npz"
    write_archive(tampered, members)
    with pytest.raises(PublicProductError, match="differs|differ"):
        verify_npz(source, tampered, expected)


def test_public_archive_hash_is_checked(accepted, tmp_path):
    source, expected = accepted
    target = tmp_path / "public.npz"
    export_npz(source, target, expected)
    with pytest.raises(PublicProductError, match="Public product SHA-256"):
        verify_npz(source, target, expected, expected_public_sha256="0" * 64)
