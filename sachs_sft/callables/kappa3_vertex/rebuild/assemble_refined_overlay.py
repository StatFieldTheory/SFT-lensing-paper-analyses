"""Assemble a new, three-row refinement overlay without numerical integration.

Parent input contract (confirm before real-data assembly)
-------------------------------------------------------
--high accepts eight native build_band NPZs (bmod/dmod, lambda_shells and
chi_shells aliases accepted) OR one summed HIGH NPZ with canonical channels.
Every HIGH input has exactly the three requested cosine triples and the full
16-shell physical-Mpc grid. Row order may differ; shell order may not.

--provenance is a JSON object with schema="refined_overlay_inputs_v1",
base_sha256, low_sha256, pk_table_high, pk_table_high_sha256, source_sha256,
potential_convention, spin2_sign_convention, base_high_source_sha256, and
high_inputs. The base source attestation is tied to base_sha256 and must equal
source_sha256. high_inputs maps
EACH HIGH-file SHA256 to {source_sha256, canoes_commit}. The commit must match
native build.canoes_commit when present. The hash binding is necessary because
historical native build files record commits, not source-content hashes.
Optional source_path/pk_path are rehashed against their declared identities.

Native HIGH build: kind="band", lo,hi,n_ell=384,n_phi=512,b_model="tree",
pk_table,pk_table_sha256,leg_order="auto",radial_measure="lambda".
Summed HIGH build: kind="refined_high", the same settings, band_windows and
band_provenance (eight records with lo,hi,n_ell,n_phi,sha256,source_sha256,
pk_table_sha256). Every source and PK identity must match the manifest.
cosmo_meta must record physical/project/lambda conventions, cosmology, and the
same potential_convention as the base. LOW may retain its historical PK; its
identity is bound independently and must not be relabelled as the HIGH PK.

Output is LOW[key] + sum(HIGH[key]) for ONLY the three requested rows. All
other arrays/values are retained; row-by-band quadrature and provenance are
explicit. Output construction is not numerical/scientific acceptance.

Example (no implicit globs):
 python assemble_refined_overlay.py --base BASE.npz --low LOW.npz \
   --high BAND1.npz BAND2.npz BAND3.npz BAND4.npz BAND5.npz BAND6.npz \
          BAND7.npz BAND8.npz --provenance identity.json --out NEW.npz
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np

CHANNELS = ('zeta_TTT', 'zeta_TTP', 'zeta_TPP', 'zeta_PPP', 'zeta_Bmod', 'zeta_Dmod')
WINDOWS = ((60, 120), (120, 240), (240, 480), (480, 960), (960, 1920),
           (1920, 3840), (3840, 7680), (7680, 15360))
DEFAULT_ROWS = (1671, 1828, 1985)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def metadata(array):
    value = json.loads(np.asarray(array, dtype=np.uint8).tobytes().decode())
    require(isinstance(value, dict), 'metadata must be a JSON object')
    return value


def encode(value):
    return np.frombuffer(json.dumps(value, sort_keys=True).encode(), dtype=np.uint8)


def sha(value):
    require(isinstance(value, str) and len(value) == 64 and
            all(c in '0123456789abcdef' for c in value), 'invalid SHA256 identity')
    return value


def load(path):
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key].copy() for key in data.files}


def field(data, canonical, alias=None):
    if canonical in data:
        value = data[canonical]
        if alias in data:
            require(np.array_equal(value, data[alias]), 'conflicting aliases: ' + canonical)
        return value
    require(alias is not None and alias in data, 'missing field: ' + canonical)
    return data[alias]


def geometry(data):
    triples = np.asarray(field(data, 'cosine_triples'))
    require(triples.ndim == 2 and triples.shape[1] == 3 and
            np.isfinite(triples).all(), 'invalid cosine triples shape/values')
    require(len({tuple(row) for row in triples}) == len(triples), 'duplicate cosine rows')
    require(np.all(abs(triples) <= 1), 'cosine outside [-1,1]')
    shell = {key: np.asarray(field(data, key, alias)) for key, alias in
             [('lambda_shells_Mpc', 'lambda_shells'), ('chi_shells_Mpc', 'chi_shells'),
              ('z_shells', None)]}
    for key, value in shell.items():
        require(value.ndim == 1 and np.isfinite(value).all(), 'invalid shell array: ' + key)
    size = len(shell['lambda_shells_Mpc'])
    require(size == 16 and all(len(x) == size for x in shell.values()), 'expected 16 shells')
    require(np.all(np.diff(shell['lambda_shells_Mpc']) > 0), 'lambda shells not increasing')
    return triples, shell


def row_map(triples, wanted):
    lookup = {tuple(row): i for i, row in enumerate(triples)}
    require(all(tuple(row) in lookup for row in wanted), 'missing exact requested cosine row')
    return np.array([lookup[tuple(row)] for row in wanted], dtype=int)


def channels(data, rows):
    result = {}
    for key in CHANNELS:
        alias = {'zeta_Bmod': 'bmod', 'zeta_Dmod': 'dmod'}.get(key)
        value = np.asarray(field(data, key, alias))
        require(value.shape == (rows, 16), 'channel shape mismatch: ' + key)
        require(value.dtype == np.dtype('float64'), 'channel must be float64: ' + key)
        require(np.isfinite(value).all(), 'non-finite channel: ' + key)
        result[key] = value
    return result


def conventions(meta, base_meta):
    for key, expected in [('units', 'physical'), ('lambda_convention', 'project'),
                          ('radial_measure', 'lambda')]:
        require(meta.get(key) == expected and base_meta.get(key) == expected,
                'convention mismatch: ' + key)
    for key in ('Omega_m', 'h'):
        require(key in meta and meta[key] == base_meta.get(key), 'cosmology mismatch: ' + key)


def assemble(base_path, low_path, high_paths, provenance_path, out_path,
             expected_rows=DEFAULT_ROWS):
    """Validate all inputs, then exclusively create one new output NPZ."""
    paths = [Path(p).resolve() for p in [base_path, low_path, *high_paths, provenance_path]]
    require(len(high_paths) in (1, 8), 'provide one summed HIGH or eight native bands')
    require(len(set(paths)) == len(paths), 'duplicate input paths')
    out = Path(out_path)
    require(out.suffix == '.npz', 'output must be an NPZ')
    require(out.resolve() not in paths and not out.exists() and not out.is_symlink(),
            'output exists or would overwrite an input')
    initial = {p: digest(p) for p in paths}
    base = load(paths[0]); low = load(paths[1])
    manifest = json.loads(paths[-1].read_text())
    require(manifest.get('schema') == 'refined_overlay_inputs_v1', 'unsupported provenance schema')
    require(sha(manifest.get('base_sha256')) == initial[paths[0]], 'base hash mismatch')
    require(sha(manifest.get('low_sha256')) == initial[paths[1]], 'LOW hash mismatch')
    source_hash = sha(manifest.get('source_sha256'))
    pk_hash = sha(manifest.get('pk_table_high_sha256'))
    require(sha(manifest.get('base_high_source_sha256')) == source_hash,
            'base/refined HIGH source mismatch')
    for key, expected in [('source_path', source_hash), ('pk_path', pk_hash)]:
        if key in manifest:
            p = Path(manifest[key]).resolve()
            require(out.resolve() != p, 'output would overwrite provenance source')
            require(digest(p) == expected, 'provenance path hash mismatch: ' + key)
            initial[p] = expected
    bm = metadata(base['cosmo_meta'])
    require(manifest.get('pk_table_high') == bm.get('pk_table_high'), 'base HIGH PK mismatch')
    if 'pk_table_high_sha256' in bm:
        require(bm['pk_table_high_sha256'] == pk_hash, 'base HIGH PK hash mismatch')
    for key in ('potential_convention', 'spin2_sign_convention'):
        require(isinstance(manifest.get(key), str) and manifest[key] == bm.get(key),
                'base sign/source mismatch: ' + key)
    bt, bs = geometry(base); bc = channels(base, len(bt))
    rows = np.asarray(expected_rows)
    require(rows.shape == (3,) and np.issubdtype(rows.dtype, np.integer) and
            len(set(rows.tolist())) == 3 and np.all((rows >= 0) & (rows < len(bt))),
            'expected exactly three distinct valid row indices')
    wanted = bt[rows]
    c = wanted[0, 1]
    require(0 < c < 1 and np.array_equal(wanted,
            np.array([[1., c, c], [c, 1., c], [c, c, 1.]])),
            'expected collapsed slots 0,1,2 at one separation')
    gamma = float(np.degrees(np.arccos(c)) * 60)
    require(abs(gamma - .5) < 1e-8, 'expected nearest 0.5 arcmin family')
    require(bm.get('band_windows') == [list(w) for w in WINDOWS] and
            bm.get('band_quadrature') == [[96, 512]], 'base quadrature is not uniform96/512')
    require(int(np.asarray(base['ell_max'])) == 15360, 'base cutoff mismatch')
    lt, ls = geometry(low); lc = channels(low, len(lt)); lm = metadata(low['cosmo_meta'])
    conventions(lm, bm)
    for key in bs:
        require(np.array_equal(ls[key], bs[key]), 'LOW shell mismatch: ' + key)
    low_rows = row_map(lt, wanted)
    low_build = metadata(low['build']) if 'build' in low else {}
    require(int(np.asarray(low['ell_max'])) == 60 and
            low_build.get('ell_cut', 60) == 60, 'LOW cutoff mismatch')
    if 'pk_table' in low_build:
        require(low_build['pk_table'] == bm.get('pk_table_low'), 'LOW PK mismatch')
    total = {key: np.zeros((3, 16), dtype=np.float64) for key in CHANNELS}
    provenance = []
    for path in paths[2:-1]:
        data = load(path); triples, shells = geometry(data)
        require(len(triples) == 3, 'HIGH must contain exactly three rows')
        selected = row_map(triples, wanted)
        values = channels(data, 3); meta = metadata(data['cosmo_meta'])
        conventions(meta, bm)
        for key, expected in [('n_ell_quad', 384), ('n_phi_quad', 512)]:
            if key in meta:
                require(meta[key] == expected, 'HIGH metadata resolution mismatch: ' + key)
        require(meta.get('potential_convention') == manifest['potential_convention'],
                'HIGH sign/source mismatch')
        if 'spin2_sign_convention' in meta:
            require(meta['spin2_sign_convention'] == manifest['spin2_sign_convention'],
                    'HIGH spin sign mismatch')
        for key in bs:
            require(np.array_equal(shells[key], bs[key]), 'HIGH shell mismatch: ' + key)
        build = metadata(data['build'])
        for key, expected in [('n_ell', 384), ('n_phi', 512), ('leg_order', 'auto'),
                              ('radial_measure', 'lambda'), ('b_model', bm.get('b_delta_model')),
                              ('pk_table', manifest['pk_table_high']),
                              ('pk_table_sha256', pk_hash)]:
            require(build.get(key) == expected, 'HIGH build mismatch: ' + key)
        record = manifest.get('high_inputs', {}).get(initial[path])
        require(isinstance(record, dict) and record.get('source_sha256') == source_hash,
                'HIGH source identity missing/mismatch')
        if 'canoes_commit' in build:
            require(record.get('canoes_commit') == build['canoes_commit'], 'HIGH commit mismatch')
        if 'source_sha256' in build:
            require(build['source_sha256'] == source_hash, 'HIGH source hash mismatch')
        if build.get('kind') == 'band':
            require(len(high_paths) == 8, 'one native band is not a summed HIGH')
            window = (build.get('lo'), build.get('hi'))
            require(window in WINDOWS and int(np.asarray(data['ell_max'])) == window[1],
                    'invalid native band window/cutoff')
            entries = [dict(lo=window[0], hi=window[1], n_ell=384, n_phi=512,
                            sha256=initial[path], source_sha256=source_hash,
                            pk_table_sha256=pk_hash)]
        else:
            require(len(high_paths) == 1 and build.get('kind') == 'refined_high',
                    'invalid summed HIGH kind')
            require(build.get('band_windows') == [list(w) for w in WINDOWS] and
                    int(np.asarray(data['ell_max'])) == 15360, 'summed HIGH windows/cutoff mismatch')
            entries = build.get('band_provenance')
            require(isinstance(entries, list) and len(entries) == 8,
                    'summed HIGH needs eight band provenance records')
            for entry in entries:
                sha(entry.get('sha256'))
                for key, expected in [('n_ell', 384), ('n_phi', 512),
                                      ('source_sha256', source_hash), ('pk_table_sha256', pk_hash)]:
                    require(entry.get(key) == expected, 'summed band identity mismatch: ' + key)
        provenance.extend(dict(entry, container_sha256=initial[path], path=str(path),
                               leg_order='auto', units='physical', radial_measure='lambda',
                               lambda_convention='project',
                               potential_convention=manifest['potential_convention'],
                               spin2_sign_convention=manifest['spin2_sign_convention'])
                          for entry in entries)
        for key in CHANNELS:
            total[key] += values[key][selected]
    require(sorted((p['lo'], p['hi']) for p in provenance) == list(WINDOWS),
            'missing/duplicate/overlapping HIGH bands')
    result = {key: value.copy() for key, value in base.items()}
    for key in CHANNELS:
        replacement = lc[key][low_rows] + total[key]
        require(np.isfinite(replacement).all(), 'non-finite assembled channel: ' + key)
        result[key][rows] = replacement
        unaffected = np.ones(len(bt), dtype=bool); unaffected[rows] = False
        require(result[key][unaffected].tobytes() == bc[key][unaffected].tobytes(),
                'unaffected values changed: ' + key)
    n_ell = np.full((len(bt), 8), 96, dtype=np.int32); n_ell[rows] = 384
    n_phi = np.full((len(bt), 8), 512, dtype=np.int32)
    require(not any(key in result for key in ('n_ell_by_row_band', 'n_phi_by_row_band',
                                             'refined_row_indices')), 'base already has overlay fields')
    result.update(n_ell_by_row_band=n_ell, n_phi_by_row_band=n_phi,
                  refined_row_indices=rows.astype(np.int64))
    overlay = dict(schema='row_band_quadrature_v1', base_sha256=initial[paths[0]],
                   low_sha256=initial[paths[1]], provenance_manifest_sha256=initial[paths[-1]],
                   row_indices=rows.tolist(), cosine_triples=wanted.tolist(), gamma_arcmin=gamma,
                   n_ell=384, n_phi=512, bands=provenance,
                   pk_table_high=manifest['pk_table_high'], pk_table_high_sha256=pk_hash,
                   source_sha256=source_hash, base_high_source_sha256=manifest['base_high_source_sha256'],
                   potential_convention=manifest['potential_convention'],
                   spin2_sign_convention=manifest['spin2_sign_convention'],
                   assembly_completed=True, numerical_acceptance='not_evaluated',
                   scope='Three rows, all six channels, all16 shells; other rows remain96/512',
                   low_rebuilt=False, uniform_quadrature=False)
    bm.update(producer='assemble_refined_overlay.py', band_quadrature=[[96, 512], [384, 512]],
              quadrature_schema='row_band_quadrature_v1', refined_overlay=overlay,
              uniform_quadrature=False)
    result['cosmo_meta'] = encode(bm)
    # Serialize before exclusively creating the output. Rehash all bound inputs
    # at the safe point; no input is rewritten and existing outputs are rejected.
    buffer = io.BytesIO(); np.savez(buffer, **result)
    for path, expected in initial.items():
        require(digest(path) == expected, 'input changed during assembly: ' + str(path))
    with out.open('xb') as stream:
        stream.write(buffer.getvalue())
    return overlay


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--base', type=Path, required=True)
    parser.add_argument('--low', type=Path, required=True)
    parser.add_argument('--high', type=Path, nargs='+', required=True)
    parser.add_argument('--provenance', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    overlay = assemble(args.base, args.low, args.high, args.provenance, args.out)
    print(json.dumps(dict(out=str(args.out), rows=overlay['row_indices'],
                          numerical_acceptance='not_evaluated')))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
