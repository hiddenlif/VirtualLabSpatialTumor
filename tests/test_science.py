import unittest
import numpy as np
from scipy.spatial.distance import cdist
from spatial_lab.data import load,region,MARKERS,digest
from spatial_lab.analysis import base_metrics,analyze,bh_adjust,terrain
from spatial_lab.reconstruction import fit,split_masks,mask_channel,spatial_predict,reconstruct
from spatial_lab.movement import simulate

class ScienceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a,cls.meta=load();cls.r=region('reg001_A');cls.model=fit(10)

    def test_source_roster_and_missing_annotations(self):
        self.assertEqual(len(self.a['cell_id']),258385)
        self.assertEqual(len(set(self.a['cell_id'].tolist())),258385)
        self.assertEqual(len(self.meta['regions']),140)
        self.assertEqual(len(set(x['patient'] for x in self.meta['regions'].values())),35)
        ix=self.meta['neighborhood_names'].index('Not assigned by authors')
        self.assertEqual(int((self.a['neighborhood']==ix).sum()),4566)
        self.assertEqual(int(sum(self.meta['exclusions'].values())),17831)

    def test_kdtree_metric_matches_direct_distances(self):
        r=self.r; xy=r['xy']; t=r['tumor']&r['eligible']; cd8=r['cd8']&r['eligible']
        direct=cdist(xy[cd8],xy[t]).min(axis=1)
        x=base_metrics(xy,t,cd8,r['eligible'])
        self.assertAlmostEqual(x['access_fraction'],float((direct<=80).mean()))
        self.assertAlmostEqual(x['median_distance_px'],float(np.median(direct)))

    def test_invalid_or_insufficient_regions_are_explicit(self):
        r=self.r.copy();r['tumor']=np.zeros(len(r['xy']),bool)
        x=analyze(r)
        self.assertFalse(x['eligible']);self.assertIsNone(x['access_fraction'])
        with self.assertRaises(ValueError):analyze(self.r,radius=81)

    def test_random_labelling_conserves_support_and_reproduces(self):
        a=analyze(self.r,permutations=199,seed=7);b=analyze(self.r,permutations=199,seed=7)
        self.assertEqual(a,b)
        expected=a['available_fraction']*a['cd8_cells']
        self.assertLess(abs(np.mean(a['global_null_counts'])-expected),2)
        for k in ('global_null_counts','block_null_counts'):
            self.assertTrue(all(0<=x<=a['cd8_cells'] for x in a[k]))
        self.assertGreaterEqual(a['block_null_p'],1/200)

    def test_bh_known_small_case(self):
        np.testing.assert_allclose(bh_adjust([.01,.04,.03]),[.03,.04,.04])

    def test_patient_splits_are_disjoint_and_complete(self):
        masks=split_masks(self.a)
        self.assertTrue(np.all(sum(m.astype(int) for m in masks)==1))
        sets=[set(self.a['patient'][m].tolist()) for m in masks]
        self.assertEqual([len(s) for s in sets],[21,7,7])
        self.assertFalse(sets[0]&sets[1] or sets[0]&sets[2] or sets[1]&sets[2])
        self.assertEqual(self.model['features'],MARKERS[1:13]);self.assertNotIn('CD8',self.model['features'])

    def test_missingness_count_and_contiguous_band(self):
        xy=np.c_[np.arange(100),np.zeros(100)]
        for kind in ('random','band'):
            m=mask_channel(100,xy,kind,101);self.assertEqual(int(m.sum()),30)
            np.testing.assert_array_equal(m,mask_channel(100,xy,kind,101))
            if kind=='band':self.assertEqual(np.ptp(np.flatnonzero(m)),29)

    def test_spatial_inverse_distance_independent_example(self):
        xy=np.array([[0,0],[2,0],[1,0.]])
        pred=spatial_predict(xy,np.array([2.,4.,999.]),np.array([False,False,True]),2)
        self.assertAlmostEqual(pred[0],3)

    def test_masked_target_cannot_influence_predictions(self):
        r=self.r;use=r['eligible'];mask=mask_channel(int(use.sum()),r['xy'][use],'random',101)
        changed=r.copy();changed['markers']=r['markers'].copy()
        ix=np.flatnonzero(use)[mask];changed['markers'][ix,0]=1e9
        for method in ('median','spatial','ridge'):
            a=reconstruct(r,method,model=self.model,include_cells=True)
            b=reconstruct(changed,method,model=self.model,include_cells=True)
            np.testing.assert_array_equal(a['reconstructed_log1p'],b['reconstructed_log1p'])
            self.assertTrue(all(x is None for x,m in zip(a['observed_log1p'],a['mask']) if m))

    def test_terrain_does_not_fill_unobserved_voids(self):
        _,support=terrain(self.r,np.array([[1e6,1e6]]));self.assertFalse(support[0])

    def test_motion_reproducible_and_preserves_source(self):
        before=self.r['xy'].copy()
        a=simulate(self.r,seed=3,agents=8,steps=10);b=simulate(self.r,seed=3,agents=8,steps=10)
        self.assertEqual(a,b);self.assertTrue(a['eligible'])
        np.testing.assert_array_equal(before,self.r['xy'])
        self.assertEqual(np.asarray(a['paths']).shape,(8,11,2))
        _,support=terrain(self.r,np.array(a['paths']).reshape(-1,2));self.assertTrue(support.all())
        with self.assertRaises(ValueError):simulate(self.r,permeability=-1)

if __name__=='__main__':unittest.main()
