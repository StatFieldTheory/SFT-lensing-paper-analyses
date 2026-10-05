"""Bounded interpolation and public-contract checks, not FK acceptance."""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


@pytest.fixture
def candidate(monkeypatch):
    module = importlib.import_module('perm_aware_kappa3_callable_pchip')
    monkeypatch.setattr(module, '_CACHE', None)
    monkeypatch.setattr(module, 'ALLOW_INTERPOLATION', False)
    return module


def install(module, grid, rows):
    module._CACHE = {'lambda': np.asarray(grid, dtype=float),
                     'per_channel': {'zeta_TTT': np.asarray(rows, dtype=float)},
                     'pchippers': {}}
    return module._CACHE


def locate_arrays(grid, queries, indices, *, exact=True, angular_weights=None):
    queries = np.asarray(queries, dtype=float)
    indices = np.asarray(indices, dtype=int)
    high = np.minimum(np.searchsorted(grid, queries, side='right'), len(grid)-1)
    low = np.maximum(high-1, 0)
    weights = (queries-grid[low])/(grid[high]-grid[low])
    aw = np.ones(indices.shape) if angular_weights is None else angular_weights
    return {'idx': indices, 'i_lo': low[:, None], 'i_hi': high[:, None],
            'w_lam': weights[:, None], 'w': aw, 'w_sum': aw.sum(axis=1),
            'exact': np.full(len(queries), exact)}


def test_mixed_batch_sign_runs_and_linear_crossings(candidate):
    grid = np.array([1., 2., 4., 7., 11.])
    rows = [[1., 2., 16., 32., 128.], [-1., -2., -16., -32., -128.],
            [1., 4., -8., -32., 0.], [0., 4., 0., -8., -16.]]
    install(candidate, grid, rows)
    q = [3., 5., 3., 1.5]
    at = locate_arrays(grid, q, np.arange(4)[:, None])
    got = candidate._read('zeta_TTT', at)
    expected = [np.exp(PchipInterpolator(grid, np.log(rows[0]))(q[0])),
                -np.exp(PchipInterpolator(grid, np.log(-np.asarray(rows[1])))(q[1])),
                -2., 2.]
    np.testing.assert_allclose(got, expected, rtol=1e-14)
    # The curve uses more than the two bracketing knots.
    assert abs(got[0]-8.) > 0.1


def test_run_stops_at_zero_and_sign_change(candidate):
    grid = np.arange(1., 9.)
    values = np.array([1., 2., 8., 0., -4., -16., 32., 128.])
    cache = install(candidate, grid, [values])
    mapping = candidate._row_pchippers(cache, 'zeta_TTT', 0)
    assert set(mapping) == {0, 1, 4, 6}
    assert mapping[0] is mapping[1]
    assert mapping[1] is not mapping[4]
    q = [3.5, 4.5, 6.5, 5.5]
    got = candidate._read('zeta_TTT', locate_arrays(grid, q, np.zeros((4,1), int)))
    np.testing.assert_allclose(got, [4., -2., 8., -8.], rtol=1e-14)


@pytest.mark.parametrize('sign', [1., -1.])
def test_monotonic_run_no_interval_overshoot(candidate, sign):
    grid = np.array([1., 2., 5., 9., 10.])
    y = sign*np.array([1., 1.1, 20., 21., 100.])
    install(candidate, grid, [y])
    q = np.linspace(grid[0], grid[-1], 301)
    got = candidate._read('zeta_TTT', locate_arrays(grid, q, np.zeros((len(q),1), int)))
    assert np.all(sign*np.diff(got) >= 0)
    for i in range(len(grid)-1):
        subset = got[(q >= grid[i]) & (q <= grid[i+1])]
        assert np.all(subset >= min(y[i],y[i+1])-1e-12)
        assert np.all(subset <= max(y[i],y[i+1])+1e-12)


def test_signed_zero_and_all_native_endpoints_bit_exact(candidate):
    grid = np.arange(1., 8.)
    values = np.array([1.3, -0., 0., -2.1, -9.7, 0., 1e-320])
    install(candidate, grid, [values])
    got = candidate._read('zeta_TTT', locate_arrays(grid, grid, np.zeros((7,1), int)))
    np.testing.assert_array_equal(got.view(np.uint64), values.view(np.uint64))


def test_four_row_angular_blend_after_native_radial(candidate):
    grid = np.array([1., 2., 4., 7.])
    rows = np.array([[1.,2.,16.,64.], [-1.,-3.,-27.,-81.],
                     [1.,6.,-10.,-20.], [0.,4.,0.,3.]])
    install(candidate, grid, rows)
    aw = np.array([[1.,2.,3.,4.], [4.,3.,2.,1.]])
    q = np.array([3., 2.5])
    got = candidate._read('zeta_TTT', locate_arrays(
        grid, q, np.tile(np.arange(4),(2,1)), exact=False, angular_weights=aw))
    native = np.column_stack([
        np.exp(PchipInterpolator(grid,np.log(rows[0]))(q)),
        -np.exp(PchipInterpolator(grid,np.log(-rows[1]))(q)),
        rows[2,1]+(q-2)/2*(rows[2,2]-rows[2,1]),
        rows[3,1]+(q-2)/2*(rows[3,2]-rows[3,1])])
    np.testing.assert_allclose(got, (aw*native).sum(axis=1)/aw.sum(axis=1), rtol=1e-14)


def write_table(path, scale, *, signed_zero=False):
    grid = np.array([1.,2.,4.,7.])
    triples = np.array([[1.,1.,1.]])
    values = scale*np.array([[1.,2.,16.,64.]])
    if signed_zero:
        values[0,1] = -0.
    payload = {ch:values.copy() for ch in ('zeta_TTT','zeta_TTP','zeta_TPP',
                                          'zeta_PPP','zeta_Bmod','zeta_Dmod')}
    import json
    meta = np.frombuffer(json.dumps({'radial_measure':'lambda'}).encode(),dtype=np.uint8)
    np.savez(path, cosine_triples=triples, lambda_shells_Mpc=grid,
             cosmo_meta=meta, has_modulus_channels_flag=True, **payload)


def test_loaded_table_cache_isolation(candidate, tmp_path, monkeypatch):
    first, second = tmp_path/'first.npz', tmp_path/'second.npz'
    write_table(first, 1.)
    write_table(second, 10.)
    monkeypatch.setattr(candidate, 'TABLE_PATH', first)
    a = candidate._ensure_cache()
    at = candidate._locate(np.array([[1.,1.,1.]]), np.array([3.]))
    va = candidate._read('zeta_TTT', at)
    assert a['pchippers']
    candidate._CACHE = None
    monkeypatch.setattr(candidate, 'TABLE_PATH', second)
    b = candidate._ensure_cache()
    assert b['pchippers'] == {}
    vb = candidate._read('zeta_TTT', candidate._locate(np.array([[1.,1.,1.]]),np.array([3.])))
    np.testing.assert_allclose(vb, 10*va, rtol=1e-14)
    assert a['pchippers'] is not b['pchippers']


def test_loaded_native_signed_zero_is_preserved(candidate, tmp_path, monkeypatch):
    path = tmp_path/'zero.npz'
    write_table(path, 1., signed_zero=True)
    monkeypatch.setattr(candidate,'TABLE_PATH',path)
    cache = candidate._ensure_cache()
    assert np.signbit(cache['per_channel']['zeta_TTT'][0,1])
    got = candidate.coupling_fn([[0.,0.,1.]]*3,[2.]*3)
    assert np.signbit(got[0,0,0])


@pytest.fixture
def actual(candidate, monkeypatch):
    import hashlib
    import json

    package = Path(__file__).resolve().parents[5]
    manifest = json.loads((package / 'reproduce/active_products.json').read_text())
    table = (package / 'reproduce' / manifest['kappa3_input']).resolve()
    assert table.is_relative_to(package) and table.is_file()
    lineage = json.loads((package / 'reproduce/public_product_lineage.json').read_text())
    expected = lineage['members'][table.relative_to(package).as_posix()]['public_sha256']
    assert hashlib.sha256(table.read_bytes()).hexdigest() == expected
    monkeypatch.setattr(candidate,'TABLE_PATH',table)
    v2 = importlib.import_module('perm_aware_kappa3_callable_v2')
    monkeypatch.setattr(v2,'TABLE_PATH',table)
    monkeypatch.setattr(v2,'_CACHE',None)
    monkeypatch.setattr(v2,'ALLOW_INTERPOLATION',False)
    return candidate,v2


def actual_directions(cache):
    rows = cache['uniq']
    selected = rows[(rows[:,0]==1)&(rows[:,1]==rows[:,2])&(rows[:,1]>.9)&(rows[:,1]<1)]
    c = float(selected[0,1])
    x = np.array([0.,0.,1.]); y = np.array([np.sqrt(1-c*c),0.,c])
    return np.array([x,x,y])


def test_public_all_actual_shell_endpoints_bit_exact(actual):
    candidate,v2 = actual
    cache = candidate._ensure_cache()
    dirs = actual_directions(cache)
    for lam in cache['lambda']:
        got = candidate.coupling_fn(dirs,[lam]*3)
        expected = v2.coupling_fn(dirs,[lam]*3)
        np.testing.assert_array_equal(got.view(np.uint64),expected.view(np.uint64))


def test_public_mixed_batch_scalar_bounds_and_offgrid(actual):
    candidate,_ = actual
    cache = candidate._ensure_cache()
    dirs = actual_directions(cache)
    lam = np.array([cache['lam_min']-1,cache['lambda'][3],
                    (cache['lambda'][5]+cache['lambda'][6])/2,cache['lam_max']+1])
    got = candidate.coupling_fn_batch(np.repeat(dirs[:,None,:],4,axis=1),
                                      np.broadcast_to(lam,(3,4)))
    assert np.isfinite(got).all()
    np.testing.assert_array_equal(got[[0,3]],np.zeros((2,3,3,3)))
    for i,t in enumerate(lam):
        np.testing.assert_array_equal(got[i],candidate.coupling_fn(dirs,[t]*3))
    x = np.array([0.,0.,1.]); angle = np.pi/10800
    y = np.array([np.sin(angle),0.,np.cos(angle)])
    with pytest.raises(ValueError,match='not tabulated rows'):
        candidate.coupling_fn([x,x,y],[cache['lambda'][4]]*3)


def test_public_optin_four_neighbor_blend(candidate, tmp_path, monkeypatch):
    """Exercise KD lookup, native radial evaluation, then tensor assembly."""
    import json
    path = tmp_path/'angular.npz'
    grid = np.array([1.,2.,4.,7.])
    triples = np.repeat(np.array([.2,.4,.6,.8])[:,None],3,axis=1)
    rows = np.array([[1.,2.,16.,64.],[-1.,-3.,-27.,-81.],
                     [1.,6.,-10.,-20.],[0.,4.,0.,3.]])
    payload = {ch:rows.copy() for ch in candidate._CHANNELS}
    meta = np.frombuffer(json.dumps({'radial_measure':'lambda'}).encode(),dtype=np.uint8)
    np.savez(path,cosine_triples=triples,lambda_shells_Mpc=grid,
             cosmo_meta=meta,has_modulus_channels_flag=True,**payload)
    monkeypatch.setattr(candidate,'TABLE_PATH',path)
    # Three unit vectors with every pair cosine c.
    c = .5
    x = np.array([0.,0.,1.])
    y = np.array([np.sqrt(1-c*c),0.,c])
    zx = (c-c*c)/y[0]
    z = np.array([zx,np.sqrt(1-c*c-zx*zx),c])
    with pytest.raises(ValueError,match='not tabulated rows'):
        candidate.coupling_fn([x,y,z],[3.]*3)
    monkeypatch.setattr(candidate,'ALLOW_INTERPOLATION',True)
    got = candidate.coupling_fn([x,y,z],[3.]*3)
    native = np.array([np.exp(PchipInterpolator(grid,np.log(rows[0]))(3.)),
                       -np.exp(PchipInterpolator(grid,np.log(-rows[1]))(3.)),
                       -2.,2.])
    weights = 1/np.linalg.norm(triples-c,axis=1)
    expected = np.dot(weights,native)/weights.sum()
    assert got[0,0,0] == pytest.approx(expected,rel=1e-14)
    assert got[1,1,1] == pytest.approx(expected,rel=1e-14)
    assert got[0,0,2] == 0


def test_nonmonotone_run_has_no_interval_overshoot(candidate):
    grid = np.array([1.,2.,4.,7.,10.])
    values = np.array([1.,10.,2.,20.,3.])
    install(candidate,grid,[values])
    for i in range(len(grid)-1):
        q = np.linspace(grid[i],grid[i+1],101)
        got = candidate._read('zeta_TTT',locate_arrays(grid,q,np.zeros((len(q),1),int)))
        assert np.all(got >= min(values[i],values[i+1])-1e-12)
        assert np.all(got <= max(values[i],values[i+1])+1e-12)
