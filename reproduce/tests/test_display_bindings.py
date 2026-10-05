"""Current publication bindings must reject missing or inconsistent inputs."""
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

PACKAGE = Path(__file__).resolve().parents[2]


def load_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, PACKAGE / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_current_operator_never_falls_back_when_both_bindings_missing():
    module = load_module('operator_inputs_test', 'figures/operator_inputs.py')
    manifest = {'status': 'ACCEPTED_FINITE_PUBLIC_SAVED_PRODUCT_SELECTION', 'products': {}}

    def forbidden():
        pytest.fail('Current selection called historical fallback')

    with pytest.raises(ValueError, match='requires its C callable and table'):
        module.selected_c_callable(lambda: manifest, forbidden, forbidden)


def test_recorded_r1_operator_policy_still_works():
    module = load_module('operator_inputs_legacy_test', 'figures/operator_inputs.py')
    legacy = SimpleNamespace(C_fn_batch=lambda: None, TABLE_PATH=Path('legacy.npz'),
                             __file__=str(PACKAGE / 'legacy.py'))
    result = module.selected_c_callable(
        lambda: {'status': 'accepted_for_revised_figures', 'products': {}},
        lambda: pytest.fail('Legacy selection resolved new product'), lambda: legacy)
    assert result[0] is legacy.C_fn_batch


def test_display_manifest_digest_binding(tmp_path, monkeypatch):
    module = load_module('publication_displays_test', 'reproduce/publication_displays.py')
    active = tmp_path / 'active_products.json'
    active.write_text('{"products": {}}')
    data = {'schema_version': 1, 'policy': 'AUTHOR_RETAINED_APPROXIMATE_DISPLAY',
            'retained_pdf_figures': {str(n): {} for n in (4, 5, 6, 17)},
            'computational_manifest_sha256': hashlib.sha256(active.read_bytes()).hexdigest()}
    manifest = tmp_path / 'publication_displays.json'
    manifest.write_text(json.dumps(data))
    monkeypatch.setattr(module, 'HERE', tmp_path)
    monkeypatch.setattr(module, 'MANIFEST', manifest)
    assert module._manifest() == data
    active.write_text('{"products": {"changed": true}}')
    with pytest.raises(ValueError, match='computational manifest'):
        module._manifest()
