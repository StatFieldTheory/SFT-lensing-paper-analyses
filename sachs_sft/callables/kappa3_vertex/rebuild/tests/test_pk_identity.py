"""Exercise P(k) provenance guards through the cached-piece CLI entry points.

The small output stub isolates assembly metadata from the external numerical
backend. These tests validate provenance and array routing, not lensing physics.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import sys
import types

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _common
import assemble
import assemble_flip_odd_high as flip
import merge_bands


SHA_A, SHA_B = "a" * 64, "b" * 64


@dataclass
class OutputStub:
    cosine_triples: np.ndarray
    lambda_shells: np.ndarray
    chi_shells: np.ndarray
    z_shells: np.ndarray
    zeta_TTT: np.ndarray
    zeta_TTP: np.ndarray
    zeta_TPP: np.ndarray
    zeta_PPP: np.ndarray
    ell_max: int
    cosmo_meta: dict

    def save_npz(self, path):
        values = dict(vars(self))
        values["cosmo_meta"] = _common.dump_meta(self.cosmo_meta)
        np.savez(path, **values)


@pytest.fixture
def stub_backend(monkeypatch):
    def combine(low, high):
        channels = {f"zeta_{channel}": getattr(low, f"zeta_{channel}") + getattr(high, f"zeta_{channel}")
                    for channel in _common.CHANNELS}
        return replace(low, **channels, ell_max=high.ell_max, cosmo_meta=dict(high.cosmo_meta))

    backend, sachs, kappa3 = (types.ModuleType(name) for name in ("canoes", "canoes.sachs", "canoes.sachs.kappa3"))
    sachs.kappa3_combine_low_high = combine
    kappa3.Kappa3Output = OutputStub
    monkeypatch.setitem(sys.modules, "canoes", backend)
    monkeypatch.setitem(sys.modules, "canoes.sachs", sachs)
    monkeypatch.setitem(sys.modules, "canoes.sachs.kappa3", kappa3)
    monkeypatch.setattr(assemble, "bootstrap_paths", lambda: None)
    monkeypatch.setattr(flip, "bootstrap_paths", lambda: None)


def build(sha, **extra):
    record = dict(lo=60, hi=120, n_ell=96, n_phi=64, b_model="tree",
                  radial_measure="lambda", pk_table="fiducial.txt", leg_order="auto")
    if sha is not None:
        record["pk_table_sha256"] = sha
    return {**record, **extra}


def write_piece(path, record, triples=None):
    triples = np.array([[0.1, 0.2, 0.3]]) if triples is None else triples
    shells = np.array([1.0, 2.0])
    shape = (len(triples), len(shells))
    np.savez(path, cosine_triples=triples, lambda_shells=shells,
             chi_shells=shells + 1, z_shells=shells / 10,
             ell_max=np.int64(record.get("hi", 60)),
             fingerprint=_common.dump_meta(_common.grid_fingerprint(triples, shells)),
             build=_common.dump_meta(record),
             cosmo_meta=_common.dump_meta({"potential_convention": flip.NEW_MARK}),
             **{f"zeta_{channel}": np.full(shape, index + 1.0)
                for index, channel in enumerate(_common.CHANNELS)},
             bmod=np.full(shape, 5.0), dmod=np.full(shape, 6.0))


def assembly_inputs(tmp_path, first_sha, second_sha, low_sha=None):
    low, first, second = (tmp_path / name for name in ("low.npz", "first.npz", "second.npz"))
    low_build = dict(ell_cut=60, radial_measure="lambda", pk_table="historical.txt")
    if low_sha is not None:
        low_build["pk_table_sha256"] = low_sha
    write_piece(low, low_build)
    write_piece(first, build(first_sha))
    write_piece(second, build(second_sha, lo=120, hi=240))
    return low, first, second


def run_assembly(module, inputs, out, monkeypatch):
    low, first, second = inputs
    monkeypatch.setattr(sys, "argv", [module.__name__, "--low", str(low), "--bands", str(first), str(second),
                                     "--cutoff", "240", "--out", str(out)])
    return module.main()


def test_identity_exposes_explicit_legacy_sentinel():
    assert _common.build_pk_identity({}) == (_common.PK_TABLE_DEFAULT, None)
    assert _common.build_pk_identity({"pk_table": "fiducial.txt", "pk_table_sha256": None}) == ("fiducial.txt", None)
    assert _common.build_pk_identity(build(SHA_A)) == ("fiducial.txt", SHA_A)


@pytest.mark.parametrize("table", ["", ".", "..", "../table.txt", "/table.txt", "a/b.txt", "a\\b.txt", "C:table.txt", None, 123])
def test_identity_rejects_non_basename(table):
    with pytest.raises(ValueError, match="basename"):
        _common.build_pk_identity({"pk_table": table})


@pytest.mark.parametrize("sha", ["", "A" * 64, "a" * 63, "z" * 64, 123, True])
def test_identity_rejects_invalid_hash(sha):
    with pytest.raises(ValueError, match="SHA-256"):
        _common.build_pk_identity({"pk_table_sha256": sha})


@pytest.mark.parametrize("module", [assemble, flip])
@pytest.mark.parametrize("second_sha", [SHA_B, None])
def test_assembly_rejects_same_basename_with_different_identity(module, second_sha, tmp_path, monkeypatch, stub_backend):
    inputs = assembly_inputs(tmp_path, SHA_A, second_sha)
    out = tmp_path / "out.npz"
    with pytest.raises(SystemExit, match=r"P\(k\) table identities"):
        run_assembly(module, inputs, out, monkeypatch)
    assert not out.exists()


@pytest.mark.parametrize("module", [assemble, flip])
@pytest.mark.parametrize("high_sha,low_sha", [(None, None), (SHA_A, None), (SHA_A, SHA_B)])
def test_assembly_retains_common_hash_and_allows_separate_low_identity(module, high_sha, low_sha, tmp_path, monkeypatch, stub_backend):
    inputs = assembly_inputs(tmp_path, high_sha, high_sha, low_sha)
    original = {path: path.read_bytes() for path in inputs}
    out = tmp_path / "out.npz"
    assert run_assembly(module, inputs, out, monkeypatch) == 0
    with np.load(out, allow_pickle=False) as product:
        metadata = _common.load_meta(product["cosmo_meta"])
        assert metadata["pk_table_high"] == "fiducial.txt"
        assert metadata["pk_table_low"] == "historical.txt"
        assert metadata.get("pk_table_high_sha256") == high_sha
        assert metadata.get("pk_table_low_sha256") == low_sha
        if high_sha is None:
            assert "pk_table_high_sha256" not in metadata
        if low_sha is None:
            assert "pk_table_low_sha256" not in metadata
        np.testing.assert_array_equal(product["zeta_TTT"], np.full((1, 2), 3.0))
    assert original == {path: path.read_bytes() for path in inputs}


def run_merge(tmp_path, first_sha, second_sha, monkeypatch):
    first, second, triples, out = (tmp_path / name for name in ("first.npz", "second.npz", "triples.npz", "merged.npz"))
    rows = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    write_piece(first, build(first_sha), rows[:1])
    write_piece(second, build(second_sha), rows[1:])
    np.savez(triples, cosine_triples=rows)
    monkeypatch.setattr(sys, "argv", ["merge_bands", "--parts", str(first), str(second),
                                     "--triples", str(triples), "--out", str(out)])
    return out


@pytest.mark.parametrize("second_sha", [SHA_B, None])
def test_merge_rejects_same_basename_with_different_identity(second_sha, tmp_path, monkeypatch):
    out = run_merge(tmp_path, SHA_A, second_sha, monkeypatch)
    with pytest.raises(SystemExit, match=r"P\(k\) table identities"):
        merge_bands.main()
    assert not out.exists()


@pytest.mark.parametrize("sha", [None, SHA_A])
def test_merge_records_common_identity_without_changing_arrays(sha, tmp_path, monkeypatch):
    out = run_merge(tmp_path, sha, sha, monkeypatch)
    assert merge_bands.main() == 0
    with np.load(out, allow_pickle=False) as product:
        record = _common.load_meta(product["build"])
        assert record["pk_table"] == "fiducial.txt"
        assert record.get("pk_table_sha256") == sha
        if sha is None:
            assert "pk_table_sha256" not in record
        np.testing.assert_array_equal(product["zeta_TTT"], np.ones((2, 2)))
