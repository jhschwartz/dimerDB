import unittest
import sys
import os
import time

sys.path.append('..')
from unredundant import RedundantDimerStructures

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()

config = {
    'paths': {
        'lib': f'{test_dir}/data/fake_lib',
        'usalign': f'{test_dir}/../../bin/USalign/USalign'
    }
}


class TestRedundantDimerStructures(unittest.TestCase):
    def test_init(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        self.assertEqual(rg.things[0], 'placeholder')
        self.assertEqual(rg.things[1], 'nothing')
        self.assertEqual(rg.threshold, 10)
        self.assertEqual(rg.config['paths']['lib'], os.path.realpath(f'{test_dir}/data/fake_lib'))


    def test_distance_unsimilar(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        dimer1name = '1xdl_A_1-1xdl_B_1'
        dimer2name = '1xdl_A_1-1xdl_D_0'
        expected = 1-0.50732 
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_distance_similar(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        dimer1name = '1xdl_A_1-1xdl_B_1'
        dimer2name = '1xdl_C_0-1xdl_D_0'
        expected = 1-0.97924
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_num_residues(self):
        pdbfile = f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdl-A-1.pdb'
        expected = 303
        result = RedundantDimerStructures._num_residues(pdbfile)
        self.assertEqual(result, expected)


    def test_dimer_coverage(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        dimername = '1xdl_A_1-1xdl_B_1'
        expected = 300.993355
        result = rg._dimer_coverage(dimername)
        self.assertTrue(abs(result-expected)<0.0001)


    def test_representative_first_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        cluster = ['4ij2_B_1-4ij2_C_1', 'zzzz_A_1-zzzz_B_1']
        rep = rg.representative(cluster)
        expected = '4ij2_B_1-4ij2_C_1'
        self.assertEqual(rep, expected)


    def test_representative_second_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        cluster = ['1xdl_A_1-1xdl_D_0', '1xdl_A_1-1xdl_B_1', '1xdl_C_0-1xdl_D_0']
        # A-B is most like the others and all should pass first criterion (coverage)
        rep = rg.representative(cluster)
        expected = '1xdl_A_1-1xdl_B_1'
        self.assertEqual(rep, expected)


    def test_representative_third_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        cluster = ['1xdl_A_1-1xdl_D_0', '1xdl_A_1-1xdl_B_1', '1xdl_C_0-1xdl_D_0', 'zzzz_X_1-zzzz_Y_1']
        # zzzz_X-zzzz_Y is identical to 1xdl_A_1-1xdl_B_1, forcing a criterion 3 tiebreaker
        rep = rg.representative(cluster)
        expected = 'zzzz_X_1-zzzz_Y_1'
        self.assertEqual(rep, expected)


    def test_tmp_assembly_file(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
        with rg._tmp_assembly_file(f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdl-A-1.pdb', f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdl-D-0.pdb') as f:
            self.assertTrue(os.path.exists(f))
            self.assertEqual(RedundantDimerStructures._num_residues(f), 303+297)


    def test_prune_redundancy_homodimers_mono(self):
        rg = RedundantDimerStructures(['1xdl_A_1-1xdl_D_0', '1xdl_A_1-1xdl_B_1', '1xdl_C_0-1xdl_D_0', '1xdl_A_1-1xdl_C_0', '4fb8_A_1-4fb8_B_1'], 0.3, config)
        results = set(rg.prune_redundancy())
        expected = set(['1xdl_A_1-1xdl_D_0', '1xdl_C_0-1xdl_D_0', '1xdl_A_1-1xdl_C_0', '4fb8_A_1-4fb8_B_1'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_homodimers_multi(self):
        rg = RedundantDimerStructures(['1xdl_A_1-1xdl_D_0', '1xdl_A_1-1xdl_B_1', '1xdl_C_0-1xdl_D_0', '1xdl_A_1-1xdl_C_0', '4fb8_A_1-4fb8_B_1'], 0.3, config)
        results = set(rg.prune_redundancy(num_workers=3))
        expected = set(['1xdl_A_1-1xdl_D_0', '1xdl_C_0-1xdl_D_0', '1xdl_A_1-1xdl_C_0', '4fb8_A_1-4fb8_B_1'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_mono(self):
        rg = RedundantDimerStructures(['3lue_B_1-3lue_M_1', '3lue_C_1-3lue_K_1', '3lue_D_1-3lue_M_1', '3lue_A_1-3lue_K_1'], 0.2, config)
        results = set(rg.prune_redundancy())
        expected = set(['3lue_D_1-3lue_M_1', '3lue_A_1-3lue_K_1'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_multi(self):
        rg = RedundantDimerStructures(['3lue_B_1-3lue_M_1', '3lue_C_1-3lue_K_1', '3lue_D_1-3lue_M_1', '3lue_A_1-3lue_K_1'], 0.2, config)
        results = set(rg.prune_redundancy(num_workers=4))
        expected = set(['3lue_D_1-3lue_M_1', '3lue_A_1-3lue_K_1'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_alternate(self):
        rg = RedundantDimerStructures(['3o8o_A_1-3o8o_B_1', '3o8o_C_1-3o8o_D_1', '3o8o_A_1-3o8o_D_1', '3o8o_B_1-3o8o_C_1'], 0.2, config)
        results = set(rg.prune_redundancy(num_workers=2))
        expected = set(['3o8o_C_1-3o8o_D_1', '3o8o_B_1-3o8o_C_1'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_parallel_spedup(self):
        rg = RedundantDimerStructures(['3o8o_A_1-3o8o_B_1', '3o8o_C_1-3o8o_D_1', '3o8o_A_1-3o8o_D_1', '3o8o_B_1-3o8o_C_1'], 0.2, config)

        t = time.time()
        results1 = set(rg.prune_redundancy(num_workers=1))
        elapsed1 = time.time() - t
        
        t = time.time()
        results2 = set(rg.prune_redundancy(num_workers=4))
        elapsed2 = time.time() - t
        
        expected = set(['3o8o_C_1-3o8o_D_1', '3o8o_B_1-3o8o_C_1'])
        self.assertEqual(results1, expected)
        self.assertEqual(results2, expected)

        # check speedup
        print(elapsed1, elapsed2)
        self.assertTrue(elapsed2 < elapsed1)


if __name__ == '__main__':
    unittest.main()

