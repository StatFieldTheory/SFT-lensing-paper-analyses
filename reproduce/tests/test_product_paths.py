"""Hash and containment checks for a relocated numerical product selection."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import pytest


SOURCE = Path(__file__).resolve().parents[1] / "product_paths.py"
SPEC = importlib.util.spec_from_file_location("test_product_paths_module", SOURCE)
assert SPEC is not None and SPEC.loader is not None
products = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(products)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def bundle(tmp_path, monkeypatch):
    root = tmp_path / "package"
    here = root / "reproduce"
    here.mkdir(parents=True)
    data = root / "products"
    data.mkdir()
    member = data / "sample.npz"
    member.write_bytes(b"unchanged numerical payload")
    manifest = {"schema_version": 1, "products": {
        "sample": {"path": "../products/sample.npz", "sha256": sha256(member)},
        "directory": {"path": "../products", "files": {"sample.npz": sha256(member)}},
    }}
    active = here / "active_products.json"

    def write():
        active.write_text(json.dumps(manifest))
        products.active_manifest.cache_clear()
        products.resolve_product.cache_clear()

    monkeypatch.setattr(products, "HERE", here)
    monkeypatch.setattr(products, "ACTIVE", active)
    write()
    yield root, manifest, write
    products.active_manifest.cache_clear()
    products.resolve_product.cache_clear()


def test_relative_selection_survives_relocation(bundle, tmp_path, monkeypatch):
    root, _, _ = bundle
    relocated = tmp_path / "relocated"
    shutil.copytree(root, relocated)
    monkeypatch.setattr(products, "HERE", relocated / "reproduce")
    monkeypatch.setattr(products, "ACTIVE", relocated / "reproduce/active_products.json")
    products.active_manifest.cache_clear()
    products.resolve_product.cache_clear()
    assert products.resolve_product("sample") == relocated / "products/sample.npz"
    assert products.resolve_product("directory") == relocated / "products"


def test_missing_manifest_has_no_fallback(bundle):
    root, _, _ = bundle
    (root / "reproduce/active_products.json").unlink()
    with pytest.raises(FileNotFoundError, match="There is no historical fallback"):
        products.resolve_product("sample")


@pytest.mark.parametrize("manifest", [{}, {"schema_version": 2, "products": {"x": {}}}])
def test_invalid_manifest_rejected(bundle, manifest):
    root, _, _ = bundle
    (root / "reproduce/active_products.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Invalid numerical product manifest"):
        products.resolve_product("sample")


@pytest.mark.parametrize("operation", ["missing", "changed"])
def test_missing_or_tampered_payload_rejected(bundle, operation):
    root, _, _ = bundle
    member = root / "products/sample.npz"
    if operation == "missing":
        member.unlink()
    else:
        member.write_bytes(b"changed payload")
    for key in ("sample", "directory"):
        with pytest.raises(ValueError, match="changed or is missing"):
            products.resolve_product(key)


@pytest.mark.parametrize("selector", ["../../outside.npz", "C:\\outside.npz", "\\\\host\\share\\outside.npz"])
def test_escaping_or_windows_absolute_selector_rejected(bundle, selector):
    _, manifest, write = bundle
    manifest["products"]["sample"]["path"] = selector
    write()
    with pytest.raises(ValueError, match="must be relative|leaves its recorded directory"):
        products.resolve_product("sample")


def test_absolute_path_is_rejected_even_inside_package(bundle):
    root, manifest, write = bundle
    manifest["products"]["sample"]["path"] = str(root / "products/sample.npz")
    write()
    with pytest.raises(ValueError, match="must be relative"):
        products.resolve_product("sample")


def test_product_symlink_escape_rejected(bundle, tmp_path):
    root, manifest, write = bundle
    outside = tmp_path / "outside.npz"
    outside.write_bytes(b"outside")
    (root / "products/linked.npz").symlink_to(outside)
    manifest["products"]["sample"] = {"path": "../products/linked.npz", "sha256": sha256(outside)}
    write()
    with pytest.raises(ValueError, match="leaves its recorded directory"):
        products.resolve_product("sample")


def test_directory_symlink_escape_rejected(bundle, tmp_path):
    root, manifest, write = bundle
    outside = tmp_path / "outside"
    outside.mkdir()
    member = outside / "sample.npz"
    member.write_bytes(b"outside")
    (root / "linked").symlink_to(outside, target_is_directory=True)
    manifest["products"]["directory"] = {
        "path": "../linked", "files": {"sample.npz": sha256(member)}}
    write()
    with pytest.raises(ValueError, match="leaves its recorded directory"):
        products.resolve_product("directory")


def test_unhashed_selection_rejected(bundle):
    _, manifest, write = bundle
    del manifest["products"]["sample"]["sha256"]
    write()
    with pytest.raises(ValueError, match="no recorded content hashes"):
        products.resolve_product("sample")


@pytest.mark.parametrize("selector", ["../reproduce/extra.npz", "/absolute.npz", "C:\\outside.npz"])
def test_directory_member_escape_rejected(bundle, selector):
    root, manifest, write = bundle
    extra = root / "reproduce/extra.npz"
    extra.write_bytes(b"inside package, outside selected directory")
    manifest["products"]["directory"]["files"] = {selector: sha256(extra)}
    write()
    with pytest.raises(ValueError, match="must be relative|leaves its recorded directory"):
        products.resolve_product("directory")


def test_directory_member_symlink_escape_rejected(bundle):
    root, manifest, write = bundle
    extra = root / "reproduce/extra.npz"
    extra.write_bytes(b"outside selected directory")
    (root / "products/linked.npz").symlink_to(extra)
    manifest["products"]["directory"]["files"] = {"linked.npz": sha256(extra)}
    write()
    with pytest.raises(ValueError, match="leaves its recorded directory"):
        products.resolve_product("directory")


def add_cutoff_selection(bundle, *, explicit=True):
    root, manifest, write = bundle
    directory = root / "cutoffs"
    directory.mkdir()
    pairs = [("tree", "auto"), ("bihalofit", "auto"), ("bihalofit", "given")]
    declared = [{"model": model, "leg_order": order, "directory": f"{model}_{order}"}
                for model, order in pairs]
    member = directory / "tree_auto/table_tree_cut960_r4_xi.npz"
    member.parent.mkdir()
    member.write_bytes(b"finite cutoff payload")
    record = {"path": "../cutoffs", "files": {str(member.relative_to(directory)): sha256(member)}}
    if explicit:
        record["series"] = declared
    manifest["products"]["cutoff_ladder"] = record
    write()
    return declared


def test_cutoff_series_and_member_preserve_explicit_scope(bundle):
    root, _, _ = bundle
    add_cutoff_selection(bundle)
    selected = products.cutoff_series()
    assert [(item["model"], item["leg_order"]) for item in selected] == [
        ("tree", "auto"), ("bihalofit", "auto"), ("bihalofit", "given")]
    assert products.resolve_cutoff_member(selected[0], 960) == (
        root / "cutoffs/tree_auto/table_tree_cut960_r4_xi.npz")
    with pytest.raises(ValueError, match="missing from the hash-declared selection"):
        products.resolve_cutoff_member(selected[0], 1920)
    with pytest.raises(ValueError, match="root must match"):
        products.resolve_cutoff_member(selected[0], 960, root=root)
    with pytest.raises(ValueError, match="Exact selected finite-cutoff scope"):
        products.resolve_cutoff_member({**selected[0], "directory": "../elsewhere"}, 960)
    for invalid in (True, 1000, 960.0):
        with pytest.raises(ValueError, match="Exact selected finite-cutoff scope"):
            products.resolve_cutoff_member(selected[0], invalid)


def test_legacy_cutoff_series_preserved(bundle):
    add_cutoff_selection(bundle, explicit=False)
    assert [(item["model"], item["leg_order"], item["directory"])
            for item in products.cutoff_series()] == [("tree", None, ""), ("bihalofit", None, "")]


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "directory", "extra_field", "unsupported"])
def test_invalid_cutoff_series_rejected(bundle, mutation):
    _, manifest, write = bundle
    declared = add_cutoff_selection(bundle)
    if mutation == "missing":
        declared.pop()
    elif mutation == "duplicate":
        declared[-1] = declared[0].copy()
    elif mutation == "directory":
        declared[0]["directory"] = "../elsewhere"
    elif mutation == "extra_field":
        declared[0]["fallback"] = True
    else:
        declared[0]["model"] = "unsupported"
    manifest["products"]["cutoff_ladder"]["series"] = declared
    write()
    with pytest.raises(ValueError, match="series"):
        products.cutoff_series()
