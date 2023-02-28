import unittest
import sys
import os
import time

sys.path.append('..')
from unredundant import RedundantDimerStructures

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
test_data = f'{test_dir}/data/unredundant_dimer_structures'

config = {
    'paths': {
        'lib': f'{test_data}/lib',
        'usalign': f'{test_dir}/../../bin/USalign/USalign',
        'resolu_file': f'{test_data}/resolu.idx'
    }
}


class TestRedundantDimerStructures(unittest.TestCase):
    def test_init(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        self.assertEqual(rg.things[0], 'placeholder')
        self.assertEqual(rg.things[1], 'nothing')
        self.assertEqual(rg.threshold, 10)
        self.assertEqual(rg.config['paths']['lib'], os.path.realpath(f'{test_data}/lib'))


    def test_distance_unsimilar(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        dimer1name = '1xdl_a1_m1_cA-1xdl_a1_m1_cB'
        dimer2name = '1xdl_a1_m1_cA-1xdl_a1_m0_cD'
        expected = 1-0.50732 
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_distance_similar(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        dimer1name = '1xdl_a1_m1_cA-1xdl_a1_m1_cB'
        dimer2name = '1xdl_a1_m0_cC-1xdl_a1_m0_cD'
        expected = 1-0.97924
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_num_residues(self):
        pdbfile = f'{test_data}/lib/rcsb/xd/1xdl-a1-m1-cA.pdb'
        expected = 303
        result = RedundantDimerStructures._num_residues(pdbfile)
        self.assertEqual(result, expected)


    def test_dimer_coverage(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        dimername = '1xdl_a1_m1_cA-1xdl_a1_m1_cB'
        expected = 300.993355
        result = rg._dimer_coverage(dimername)
        self.assertTrue(abs(result-expected)<0.0001)


    def test_representative_first_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        cluster = ['4ij2_a1_m1_cB-4ij2_a1_m1_cC', 'zzzz_a1_m1_cA-zzzz_a1_m1_cB']
        rep = rg.representative(cluster)
        expected = '4ij2_a1_m1_cB-4ij2_a1_m1_cC'
        self.assertEqual(rep, expected)


    def test_representative_second_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        cluster = ['1xdl_a1_m1_cA-1xdl_a1_m0_cD', '1xdl_a1_m1_cA-1xdl_a1_m1_cB', '1xdl_a1_m0_cC-1xdl_a1_m0_cD']
        # A-B is most like the others and all should pass first criterion (coverage)
        rep = rg.representative(cluster)
        expected = '1xdl_a1_m1_cA-1xdl_a1_m1_cB'
        self.assertEqual(rep, expected)


    def test_representative_third_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        cluster = ['1xdl_a1_m1_cA-1xdl_a1_m0_cD', '1xdl_a1_m1_cA-1xdl_a1_m1_cB', '1xdl_a1_m0_cC-1xdl_a1_m0_cD', 'zzzz_a1_m1_cX-zzzz_a1_m1_cY']
        # zzzz_X-zzzz_Y is identical to 1xdl_a1_m1_cA-1xdl_a1_m1_cB, forcing a criterion 3 tiebreaker
        rep = rg.representative(cluster)
        expected = 'zzzz_a1_m1_cX-zzzz_a1_m1_cY'
        self.assertEqual(rep, expected)


    def test_tmp_assembly_file(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        with rg._tmp_assembly_file(f'{test_data}/lib/rcsb/xd/1xdl-a1-m1-cA.pdb', f'{test_data}/lib/rcsb/xd/1xdl-a1-m0-cD.pdb') as f:
            self.assertTrue(os.path.exists(f))
            self.assertEqual(RedundantDimerStructures._num_residues(f), 303+297)


    def test_prune_redundancy_homodimers_mono(self):
        rg = RedundantDimerStructures(['1xdl_a1_m1_cA-1xdl_a1_m0_cD', '1xdl_a1_m1_cA-1xdl_a1_m1_cB', '1xdl_a1_m0_cC-1xdl_a1_m0_cD', '1xdl_a1_m1_cA-1xdl_a1_m0_cC', '4fb8_a1_m1_cA-4fb8_a1_m1_cB'], 0.3, config)
        results = set(rg.prune_redundancy())
        expected = set(['1xdl_a1_m1_cA-1xdl_a1_m0_cD', '1xdl_a1_m1_cA-1xdl_a1_m1_cB', '1xdl_a1_m1_cA-1xdl_a1_m0_cC', '4fb8_a1_m1_cA-4fb8_a1_m1_cB'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_homodimers_multi(self):
        rg = RedundantDimerStructures(['1xdl_a1_m1_cA-1xdl_a1_m0_cD', '1xdl_a1_m1_cA-1xdl_a1_m1_cB', '1xdl_a1_m0_cC-1xdl_a1_m0_cD', '1xdl_a1_m1_cA-1xdl_a1_m0_cC', '4fb8_a1_m1_cA-4fb8_a1_m1_cB'], 0.3, config)
        results = set(rg.prune_redundancy(num_workers=3))
        expected = set(['1xdl_a1_m1_cA-1xdl_a1_m0_cD', '1xdl_a1_m1_cA-1xdl_a1_m1_cB', '1xdl_a1_m1_cA-1xdl_a1_m0_cC', '4fb8_a1_m1_cA-4fb8_a1_m1_cB'])
        self.assertEqual(results, expected)

    
    def test_get_chain_resolu_xray(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        result = rg._get_chain_resolu('1xdl_a1_m1_cA')
        expected = 3.0
        self.assertEqual(result, expected)


    def test_get_chain_resolu_notxray(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        result = rg._get_chain_resolu('1acz_a1_m1_cA')
        expected = -1.0
        self.assertEqual(result, expected)

    def test_get_chain_resolu_notxray2(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        result = rg._get_chain_resolu('9999_a1_m1_cA')
        expected = -1.0
        self.assertEqual(result, expected)

    def test_get_chain_resolu_fail(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        with self.assertRaises(KeyError):
            rg._get_chain_resolu('9zzz_a1_m1_cA')


    def test_get_dimer_avg_resolu_both_xray(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        result = rg._get_dimer_avg_resolu('1xdl_a1_m1_cA-4zza_a1_m1_cA')
        expected = 2.46171
        self.assertTrue(abs(result-expected)<0.0001)


    def test_get_dimer_avg_resolu_non_xray_right(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        result = rg._get_dimer_avg_resolu('1xdl_a1_m1_cA-1acz_a1_m1_cA')
        expected = -1.0
        self.assertEqual(result, expected)


    def test_get_dimer_avg_resolu_non_xray_left(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        result = rg._get_dimer_avg_resolu('1acz_a1_m1_cA-1xdl_a1_m1_cA')
        expected = -1.0
        self.assertEqual(result, expected)


    def test_get_dimer_avg_resolu_non_xray_both(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        result = rg._get_dimer_avg_resolu('4znf_a2_m4_cZ-1acz_a1_m1_cA')
        expected = -1.0
        self.assertEqual(result, expected)


    def test_representative_criterion_one_point_five(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        cluster = ['1xdl_a1_m1_cA-1xdl_a1_m0_cD', '9999_a1_m1_cA-9999_a1_m1_cB']
        # have manually made 9999 identical to 1xdl, but 9999 is non xray in resolu.idx
        rep = rg.representative(cluster, prefer_xray=True)
        expected = '1xdl_a1_m1_cA-1xdl_a1_m0_cD'
        self.assertEqual(rep, expected)


    def test_representative_criterion_two_point_five(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        cluster = ['1xdl_a1_m1_cA-1xdl_a1_m1_cB', 'zzzz_a1_m1_cX-zzzz_a1_m1_cY']
        # zzzz_X-zzzz_Y is identical to 1xdl_a1_m1_cA-1xdl_a1_m1_cB, forcing us to reach criterion 2.5 - zzzz has been given very bad resolution manually in the resolu.idx file
        rep = rg.representative(cluster, prefer_xray=True)
        expected = '1xdl_a1_m1_cA-1xdl_a1_m1_cB'
        self.assertEqual(rep, expected)


if __name__ == '__main__':
    unittest.main()

