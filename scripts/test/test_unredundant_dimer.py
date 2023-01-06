import unittest
import sys
import os

sys.path.append('..')
from unredundant import RedundantDimers

import yaml
with open('data/fake_config.yaml', 'r') as f:
    config = yaml.safe_load(f)


class TestRedundantDimers(unittest.TestCase):
    def test_init(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, config)
        self.assertEqual(rg.things[0], 'placeholder')
        self.assertEqual(rg.things[1], 'nothing')
        self.assertEqual(rg.threshold, 10)
        self.assertEqual(rg.config['paths']['lib'], 'data/fake_lib')


    def test_distance_unsimilar(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, config)
        dimer1name = '1xdl_A-1xdl_B'
        dimer2name = '1xdl_A-1xdl_D'
        expected = 1-0.50732 
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_distance_similar(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, config)
        dimer1name = '1xdl_A-1xdl_B'
        dimer2name = '1xdl_C-1xdl_D'
        expected = 1-0.97924
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_num_residues(self):
        pdbfile = 'data/fake_lib/rcsb_pdb/xd/1xdlA.pdb'
        expected = 303
        result = RedundantDimers._num_residues(pdbfile)
        self.assertEqual(result, expected)


    def test_dimer_coverage(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, config)
        dimername = '1xdl_A-1xdl_B'
        expected = 300.993355
        result = rg._dimer_coverage(dimername)
        self.assertTrue(abs(result-expected)<0.0001)


    def test_representative_first_criterion(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, config)
        cluster = ['4ij2_B-4ij2_C', 'zzzz_A-zzzz_B']
        rep = rg.representative(cluster)
        expected = '4ij2_B-4ij2_C'
        self.assertEqual(rep, expected)


    def test_representative_second_criterion(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, config)
        cluster = ['1xdl_A-1xdl_D', '1xdl_A-1xdl_B', '1xdl_C-1xdl_D']
        # A-B is most like the others and all should pass first criterion (coverage)
        rep = rg.representative(cluster)
        expected = '1xdl_A-1xdl_B'
        self.assertEqual(rep, expected)
    

    def test_representative_third_criterion(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, config)
        cluster = ['1xdl_A-1xdl_D', '1xdl_A-1xdl_B', '1xdl_C-1xdl_D', 'zzzz_X-zzzz_Y']
        # zzzz_X-zzzz_Y is identical to 1xdl_A-1xdl_B, forcing a criterion 3 tiebreaker
        rep = rg.representative(cluster)
        expected = 'zzzz_X-zzzz_Y'
        self.assertEqual(rep, expected)


    def test_tmp_assembly_file(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, config)
        with rg._tmp_assembly_file('data/fake_lib/rcsb_pdb/xd/1xdlA.pdb', 'data/fake_lib/rcsb_pdb/xd/1xdlD.pdb') as f:
            self.assertTrue(os.path.exists(f))
            self.assertEqual(RedundantDimers._num_residues(f), 303+297)


    def test_prune_redundancy_homodimers_mono(self):
        rg = RedundantDimers(['1xdl_A-1xdl_D', '1xdl_A-1xdl_B', '1xdl_C-1xdl_D', '1xdl_A-1xdl_C', '4fb8_A-4fb8_B'], 0.3, config)
        results = set(rg.prune_redundancy())
        expected = set(['1xdl_A-1xdl_D', '1xdl_C-1xdl_D', '1xdl_A-1xdl_C', '4fb8_A-4fb8_B'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_homodimers_multi(self):
        rg = RedundantDimers(['1xdl_A-1xdl_D', '1xdl_A-1xdl_B', '1xdl_C-1xdl_D', '1xdl_A-1xdl_C', '4fb8_A-4fb8_B'], 0.3, config)
        results = set(rg.prune_redundancy(num_workers=3))
        expected = set(['1xdl_A-1xdl_D', '1xdl_C-1xdl_D', '1xdl_A-1xdl_C', '4fb8_A-4fb8_B'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_mono(self):
        rg = RedundantDimers(['3lue_B-3lue_M', '3lue_C-3lue_K', '3lue_D-3lue_M', '3lue_A-3lue_K'], 0.2, config)
        results = set(rg.prune_redundancy())
        expected = set(['3lue_D-3lue_M', '3lue_A-3lue_K'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_multi(self):
        rg = RedundantDimers(['3lue_B-3lue_M', '3lue_C-3lue_K', '3lue_D-3lue_M', '3lue_A-3lue_K'], 0.2, config)
        results = set(rg.prune_redundancy(num_workers=4))
        expected = set(['3lue_D-3lue_M', '3lue_A-3lue_K'])
        self.assertEqual(results, expected)


if __name__ == '__main__':
    unittest.main()

