"""Bind operator figures to explicitly selected callable and table source bytes."""
import hashlib
import importlib.util
from pathlib import Path


def selected_c_callable(active_manifest, resolve_product, legacy):
    manifest = active_manifest()
    records = manifest['products']
    keys = ('corr_op_callable', 'corr_op_table')
    present = [key in records for key in keys]
    if not any(present):
        if (manifest.get('status') != 'accepted_for_revised_figures'
                or 'public_lineage' in manifest):
            raise ValueError('Current operator selection requires its C callable and table')
        module = legacy()
        return module.C_fn_batch, Path(module.TABLE_PATH), {
            'scope': 'RECORDED_LEGACY_C21_OPERATOR_FIGURE_INPUT',
            'callable': str(Path(module.__file__).resolve()),
            'table': str(Path(module.TABLE_PATH).resolve()),
        }
    if not all(present):
        raise ValueError('Selected operator callable and table must both be declared')
    source = resolve_product(keys[0])
    table = resolve_product(keys[1])
    with table.open('rb') as stream:
        if hashlib.file_digest(stream, 'sha256').hexdigest() != records[keys[1]].get('sha256'):
            raise ValueError('Selected operator table changed before import')
    raw = source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != records[keys[0]].get('sha256'):
        raise ValueError('Selected operator callable source changed before import')
    name = 'selected_operator_C_' + records[keys[0]]['sha256'][:16]
    spec = importlib.util.spec_from_file_location(name, source)
    module = importlib.util.module_from_spec(spec)
    exec(compile(raw, str(source), 'exec'), module.__dict__)
    if (Path(module.TABLE_PATH).resolve() != table
            or module.TABLE_SHA256 != records[keys[1]].get('sha256')):
        raise ValueError('Selected callable does not bind the exact selected operator table')
    if not callable(module.C_fn_batch):
        raise ValueError('Selected callable lacks its original vectorized-C interface')
    with table.open('rb') as stream:
        if hashlib.file_digest(stream, 'sha256').hexdigest() != records[keys[1]]['sha256']:
            raise ValueError('Selected operator table changed during import')
    final_raw = source.read_bytes()
    if (final_raw != raw
            or hashlib.sha256(final_raw).hexdigest() != records[keys[0]]['sha256']):
        raise ValueError('Selected operator callable source changed during import')
    return module.C_fn_batch, table, {
        'scope': 'SELECTED_C23_OPERATOR_FIGURE_INPUT_NO_GLOBAL_ACCURACY_GRANT',
        'callable': {'path': str(source), 'sha256': records[keys[0]]['sha256']},
        'table': {'path': str(table), 'sha256': records[keys[1]]['sha256']},
        'source_import': 'CHECKED_RAW_BYTES_NO_PYC',
    }
