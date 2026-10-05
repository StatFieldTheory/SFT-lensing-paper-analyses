"""Small synthetic tests for the owned refinement-overlay assembler."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import assemble_refined_overlay as asm


class OverlayTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='stf_overlay_test_')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.base = self.root/'base.npz'
        self.low = self.root/'low.npz'
        self.out = self.root/'new.npz'
        self.manifest = self.root/'provenance.json'
        self.rows = (1, 3, 5)
        c = np.cos(np.deg2rad(.5/60))
        self.triples = np.array([[.1,.2,.3], [1.,c,c], [.4,.5,.6],
                                 [c,1.,c], [.7,.8,.9], [c,c,1.]])
        self.shells = dict(lambda_shells_Mpc=np.arange(1.,17.),
                           chi_shells_Mpc=np.arange(21.,37.), z_shells=np.arange(16.)/10)
        self.meta = dict(units='physical', lambda_convention='project', radial_measure='lambda',
                         Omega_m=.3160919980475834, h=.6711, pk_table_high='fid.txt',
                         pk_table_low='historical.txt', potential_convention='plus Psi0',
                         spin2_sign_convention='plus_exp2iphi', b_delta_model='tree',
                         band_windows=[list(w) for w in asm.WINDOWS], band_quadrature=[[96,512]])
        base = dict(cosine_triples=self.triples, **self.shells,
                    cosmo_meta=asm.encode(self.meta), ell_max=np.int64(15360),
                    has_modulus_channels_flag=np.bool_(True), unrelated=np.array([123],np.int32))
        for i,key in enumerate(asm.CHANNELS):
            base[key] = np.arange(96.,dtype=np.float64).reshape(6,16) + i
            base[key][0,0] = -0.
        np.savez(self.base, **base)
        self.original = asm.load(self.base)
        low = dict(cosine_triples=self.triples[::-1],
                   lambda_shells=self.shells['lambda_shells_Mpc'],
                   chi_shells=self.shells['chi_shells_Mpc'], z_shells=self.shells['z_shells'],
                   cosmo_meta=asm.encode(self.meta), ell_max=np.int64(60),
                   build=asm.encode(dict(ell_cut=60,pk_table='historical.txt')))
        for i,key in enumerate(asm.CHANNELS):
            alias={'zeta_Bmod':'bmod','zeta_Dmod':'dmod'}.get(key,key)
            low[alias]=np.full((6,16),5.+i)
        np.savez(self.low, **low)
        self.high=[]
        for lo,hi in asm.WINDOWS:
            p=self.root/f'band_{lo}_{hi}.npz'
            d=dict(cosine_triples=self.triples[list(self.rows)][::-1],
                   lambda_shells=self.shells['lambda_shells_Mpc'],
                   chi_shells=self.shells['chi_shells_Mpc'], z_shells=self.shells['z_shells'],
                   cosmo_meta=asm.encode(dict(self.meta,n_ell_quad=384,n_phi_quad=512)),
                   ell_max=np.int64(hi),build=asm.encode(dict(kind='band',lo=lo,hi=hi,
                       n_ell=384,n_phi=512,b_model='tree',leg_order='auto',radial_measure='lambda',
                       pk_table='fid.txt',pk_table_sha256='a'*64,canoes_commit='native-commit')))
            for i,key in enumerate(asm.CHANNELS):
                alias={'zeta_Bmod':'bmod','zeta_Dmod':'dmod'}.get(key,key)
                d[alias]=np.full((3,16),1.+i)
            np.savez(p,**d);self.high.append(p)
        self.identity=dict(schema='refined_overlay_inputs_v1',base_sha256=asm.digest(self.base),
                           low_sha256=asm.digest(self.low),source_sha256='b'*64,
                           base_high_source_sha256='b'*64,pk_table_high='fid.txt',
                           pk_table_high_sha256='a'*64,potential_convention='plus Psi0',
                           spin2_sign_convention='plus_exp2iphi')
        self.write_manifest()

    def write_manifest(self):
        self.identity['high_inputs']={asm.digest(p):dict(source_sha256='b'*64,
            canoes_commit='native-commit') for p in self.high}
        self.manifest.write_text(json.dumps(self.identity))

    def mutate(self,path,fn):
        d=asm.load(path);fn(d);np.savez(path,**d);self.write_manifest()

    def run_assembly(self):
        return asm.assemble(self.base,self.low,self.high,self.manifest,self.out,self.rows)

    def assert_reject(self,pattern):
        with self.assertRaisesRegex(ValueError,pattern):self.run_assembly()
        self.assertFalse(self.out.exists())

    def test_native_preserves_unaffected_bitwise(self):
        inputs={p:asm.digest(p) for p in [self.base,self.low,*self.high,self.manifest]}
        info=self.run_assembly();new=asm.load(self.out)
        for i,key in enumerate(asm.CHANNELS):
            self.assertEqual(new[key][[0,2,4]].tobytes(),self.original[key][[0,2,4]].tobytes())
            np.testing.assert_array_equal(new[key][list(self.rows)],np.full((3,16),13.+9*i))
        for key in self.original:
            if key not in (*asm.CHANNELS,'cosmo_meta'):
                self.assertEqual(new[key].dtype,self.original[key].dtype)
                self.assertEqual(new[key].tobytes(),self.original[key].tobytes())
        self.assertEqual(inputs,{p:asm.digest(p) for p in inputs})
        self.assertEqual(info['numerical_acceptance'],'not_evaluated')
        np.testing.assert_array_equal(new['n_ell_by_row_band'][list(self.rows)],384)
        np.testing.assert_array_equal(new['n_ell_by_row_band'][[0,2,4]],96)
        self.assertEqual(len(info['bands']),8)
        self.assertEqual(asm.metadata(new['cosmo_meta'])['band_quadrature'],[[96,512],[384,512]])

    def test_summed_canonical_format(self):
        d=asm.load(self.high[0]);build=asm.metadata(d['build'])
        build.update(kind='refined_high',band_windows=[list(w) for w in asm.WINDOWS],
            band_provenance=[dict(lo=lo,hi=hi,n_ell=384,n_phi=512,sha256=asm.digest(p),
                                 source_sha256='b'*64,pk_table_sha256='a'*64)
                             for (lo,hi),p in zip(asm.WINDOWS,self.high)])
        d['build']=asm.encode(build);d['ell_max']=np.int64(15360)
        d['lambda_shells_Mpc']=d.pop('lambda_shells');d['chi_shells_Mpc']=d.pop('chi_shells')
        for key in asm.CHANNELS:
            alias={'zeta_Bmod':'bmod','zeta_Dmod':'dmod'}.get(key,key)
            d[key]=d.pop(alias)*8
        summed=self.root/'summed.npz';np.savez(summed,**d);self.high=[summed];self.write_manifest()
        self.run_assembly()
        np.testing.assert_array_equal(asm.load(self.out)['zeta_Bmod'][list(self.rows)],49.)

    def test_missing_channel(self):
        self.mutate(self.high[0],lambda d:d.pop('bmod'));self.assert_reject('missing field')

    def test_duplicate_rows(self):
        def bad(d):d['cosine_triples'][1]=d['cosine_triples'][0]
        self.mutate(self.high[0],bad);self.assert_reject('duplicate cosine')

    def test_missing_exact_row(self):
        def bad(d):d['cosine_triples'][0,0]=np.nextafter(d['cosine_triples'][0,0],0.)
        self.mutate(self.high[0],bad);self.assert_reject('missing exact')

    def test_channel_shape(self):
        self.mutate(self.high[0],lambda d:d.update(zeta_TTT=np.zeros((3,15))))
        self.assert_reject('channel shape')

    def test_nonfinite_channel(self):
        def bad(d):d['dmod'][0,0]=np.nan
        self.mutate(self.high[0],bad);self.assert_reject('non-finite')

    def test_shell_order_and_exactness(self):
        def bad(d):d['z_shells']=d['z_shells'][::-1]
        self.mutate(self.high[0],bad);self.assert_reject('HIGH shell mismatch')

    def test_source_mismatch(self):
        data=json.loads(self.manifest.read_text());data['high_inputs'][asm.digest(self.high[0])]['source_sha256']='c'*64
        self.manifest.write_text(json.dumps(data));self.assert_reject('source identity')

    def test_base_source_mismatch(self):
        self.identity['base_high_source_sha256']='c'*64;self.write_manifest()
        self.assert_reject('base/refined HIGH source')

    def test_pk_mix(self):
        def bad(d):b=asm.metadata(d['build']);b['pk_table_sha256']='c'*64;d['build']=asm.encode(b)
        self.mutate(self.high[0],bad);self.assert_reject('pk_table_sha256')

    def test_sign_mismatch(self):
        def bad(d):m=asm.metadata(d['cosmo_meta']);m['potential_convention']='minus Psi0';d['cosmo_meta']=asm.encode(m)
        self.mutate(self.high[0],bad);self.assert_reject('HIGH sign/source')

    def test_duplicate_band(self):
        def bad(d):b=asm.metadata(d['build']);b.update(lo=60,hi=120);d['build']=asm.encode(b);d['ell_max']=np.int64(120)
        self.mutate(self.high[1],bad);self.assert_reject('missing/duplicate/overlapping')

    def test_resolution_mismatch(self):
        def bad(d):b=asm.metadata(d['build']);b['n_ell']=192;d['build']=asm.encode(b)
        self.mutate(self.high[0],bad);self.assert_reject('n_ell')

    def test_input_overwrite_and_existing_output(self):
        original=asm.digest(self.base)
        with self.assertRaisesRegex(ValueError,'overwrite'):
            asm.assemble(self.base,self.low,self.high,self.manifest,self.base,self.rows)
        self.assertEqual(original,asm.digest(self.base))
        self.out.write_bytes(b'keep existing artifact')
        with self.assertRaisesRegex(ValueError,'output exists'):self.run_assembly()
        self.assertEqual(self.out.read_bytes(),b'keep existing artifact')

    def test_commit_mismatch(self):
        data=json.loads(self.manifest.read_text())
        data['high_inputs'][asm.digest(self.high[0])]['canoes_commit']='different'
        self.manifest.write_text(json.dumps(data));self.assert_reject('HIGH commit mismatch')

    def test_summed_source_provenance_mismatch(self):
        d=asm.load(self.high[0]);build=asm.metadata(d['build'])
        build.update(kind='refined_high',band_windows=[list(w) for w in asm.WINDOWS],
            band_provenance=[dict(lo=lo,hi=hi,n_ell=384,n_phi=512,sha256=asm.digest(p),
                                 source_sha256='c'*64,pk_table_sha256='a'*64)
                             for (lo,hi),p in zip(asm.WINDOWS,self.high)])
        d['build']=asm.encode(build);d['ell_max']=np.int64(15360)
        summed=self.root/'summed_bad.npz';np.savez(summed,**d)
        self.high=[summed];self.write_manifest();self.assert_reject('summed band identity')

    def test_optional_source_file_hash(self):
        source=self.root/'source.py';source.write_text('changed source')
        self.identity['source_path']=str(source);self.write_manifest()
        self.assert_reject('provenance path hash')

    def test_wrong_row_count(self):
        with self.assertRaisesRegex(ValueError,'exactly three'):
            asm.assemble(self.base,self.low,self.high,self.manifest,self.out,(1,3))
        self.assertFalse(self.out.exists())

    def test_low_hash_mismatch(self):
        self.identity['low_sha256']='c'*64;self.write_manifest();self.assert_reject('LOW hash')


if __name__=='__main__':
    unittest.main()
