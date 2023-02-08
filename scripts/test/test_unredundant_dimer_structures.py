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
<<<<<<< HEAD
        'mmalign_exe': f'{test_dir}/../../bin/MMalign'
=======
        'usalign': f'{test_dir}/../../bin/USalign/USalign'
>>>>>>> addmodelsplittopdb
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
<<<<<<< HEAD
        dimer1name = '1xdl_A-1xdl_B'
        dimer2name = '1xdl_A-1xdl_D'
=======
        dimer1name = '1xdl_A_1-1xdl_B_1'
        dimer2name = '1xdl_A_1-1xdl_D_0'
>>>>>>> addmodelsplittopdb
        expected = 1-0.50732 
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_distance_similar(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
<<<<<<< HEAD
        dimer1name = '1xdl_A-1xdl_B'
        dimer2name = '1xdl_C-1xdl_D'
=======
        dimer1name = '1xdl_A_1-1xdl_B_1'
        dimer2name = '1xdl_C_0-1xdl_D_0'
>>>>>>> addmodelsplittopdb
        expected = 1-0.97924
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_num_residues(self):
<<<<<<< HEAD
        pdbfile = f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdlA.pdb'
=======
        pdbfile = f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdl-A-1.pdb'
>>>>>>> addmodelsplittopdb
        expected = 303
        result = RedundantDimerStructures._num_residues(pdbfile)
        self.assertEqual(result, expected)


    def test_dimer_coverage(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
<<<<<<< HEAD
        dimername = '1xdl_A-1xdl_B'
=======
        dimername = '1xdl_A_1-1xdl_B_1'
>>>>>>> addmodelsplittopdb
        expected = 300.993355
        result = rg._dimer_coverage(dimername)
        self.assertTrue(abs(result-expected)<0.0001)


    def test_representative_first_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
<<<<<<< HEAD
        cluster = ['4ij2_B-4ij2_C', 'zzzz_A-zzzz_B']
        rep = rg.representative(cluster)
        expected = '4ij2_B-4ij2_C'
=======
        cluster = ['4ij2_B_1-4ij2_C_1', 'zzzz_A_1-zzzz_B_1']
        rep = rg.representative(cluster)
        expected = '4ij2_B_1-4ij2_C_1'
>>>>>>> addmodelsplittopdb
        self.assertEqual(rep, expected)


    def test_representative_second_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
<<<<<<< HEAD
        cluster = ['1xdl_A-1xdl_D', '1xdl_A-1xdl_B', '1xdl_C-1xdl_D']
        # A-B is most like the others and all should pass first criterion (coverage)
        rep = rg.representative(cluster)
        expected = '1xdl_A-1xdl_B'
=======
        cluster = ['1xdl_A_1-1xdl_D_0', '1xdl_A_1-1xdl_B_1', '1xdl_C_0-1xdl_D_0']
        # A-B is most like the others and all should pass first criterion (coverage)
        rep = rg.representative(cluster)
        expected = '1xdl_A_1-1xdl_B_1'
>>>>>>> addmodelsplittopdb
        self.assertEqual(rep, expected)


    def test_representative_third_criterion(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
<<<<<<< HEAD
        cluster = ['1xdl_A-1xdl_D', '1xdl_A-1xdl_B', '1xdl_C-1xdl_D', 'zzzz_X-zzzz_Y']
        # zzzz_X-zzzz_Y is identical to 1xdl_A-1xdl_B, forcing a criterion 3 tiebreaker
        rep = rg.representative(cluster)
        expected = 'zzzz_X-zzzz_Y'
=======
        cluster = ['1xdl_A_1-1xdl_D_0', '1xdl_A_1-1xdl_B_1', '1xdl_C_0-1xdl_D_0', 'zzzz_X_1-zzzz_Y_1']
        # zzzz_X-zzzz_Y is identical to 1xdl_A_1-1xdl_B_1, forcing a criterion 3 tiebreaker
        rep = rg.representative(cluster)
        expected = 'zzzz_X_1-zzzz_Y_1'
>>>>>>> addmodelsplittopdb
        self.assertEqual(rep, expected)


    def test_tmp_assembly_file(self):
        rg = RedundantDimerStructures(['placeholder', 'nothing'], 10, config)
<<<<<<< HEAD
        with rg._tmp_assembly_file(f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdlA.pdb', f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdlD.pdb') as f:
=======
        with rg._tmp_assembly_file(f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdl-A-1.pdb', f'{test_dir}/data/fake_lib/rcsb_pdb/xd/1xdl-D-0.pdb') as f:
>>>>>>> addmodelsplittopdb
            self.assertTrue(os.path.exists(f))
            self.assertEqual(RedundantDimerStructures._num_residues(f), 303+297)


    def test_prune_redundancy_homodimers_mono(self):
<<<<<<< HEAD
        rg = RedundantDimerStructures(['1xdl_A-1xdl_D', '1xdl_A-1xdl_B', '1xdl_C-1xdl_D', '1xdl_A-1xdl_C', '4fb8_A-4fb8_B'], 0.3, config)
        results = set(rg.prune_redundancy())
        expected = set(['1xdl_A-1xdl_D', '1xdl_C-1xdl_D', '1xdl_A-1xdl_C', '4fb8_A-4fb8_B'])
=======
        rg = RedundantDimerStructures(['1xdl_A_1-1xdl_D_0', '1xdl_A_1-1xdl_B_1', '1xdl_C_0-1xdl_D_0', '1xdl_A_1-1xdl_C_0', '4fb8_A_1-4fb8_B_1'], 0.3, config)
        results = set(rg.prune_redundancy())
        expected = set(['1xdl_A_1-1xdl_D_0', '1xdl_C_0-1xdl_D_0', '1xdl_A_1-1xdl_C_0', '4fb8_A_1-4fb8_B_1'])
>>>>>>> addmodelsplittopdb
        self.assertEqual(results, expected)


    def test_prune_redundancy_homodimers_multi(self):
<<<<<<< HEAD
        rg = RedundantDimerStructures(['1xdl_A-1xdl_D', '1xdl_A-1xdl_B', '1xdl_C-1xdl_D', '1xdl_A-1xdl_C', '4fb8_A-4fb8_B'], 0.3, config)
        results = set(rg.prune_redundancy(num_workers=3))
        expected = set(['1xdl_A-1xdl_D', '1xdl_C-1xdl_D', '1xdl_A-1xdl_C', '4fb8_A-4fb8_B'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_mono(self):
        rg = RedundantDimerStructures(['3lue_B-3lue_M', '3lue_C-3lue_K', '3lue_D-3lue_M', '3lue_A-3lue_K'], 0.2, config)
        results = set(rg.prune_redundancy())
        expected = set(['3lue_D-3lue_M', '3lue_A-3lue_K'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_multi(self):
        rg = RedundantDimerStructures(['3lue_B-3lue_M', '3lue_C-3lue_K', '3lue_D-3lue_M', '3lue_A-3lue_K'], 0.2, config)
        results = set(rg.prune_redundancy(num_workers=4))
        expected = set(['3lue_D-3lue_M', '3lue_A-3lue_K'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_alternate(self):
        rg = RedundantDimerStructures(['3o8o_A-3o8o_B', '3o8o_C-3o8o_D', '3o8o_A-3o8o_D', '3o8o_B-3o8o_C'], 0.2, config)
        results = set(rg.prune_redundancy(num_workers=2))
        expected = set(['3o8o_C-3o8o_D', '3o8o_B-3o8o_C'])
        self.assertEqual(results, expected)


    def test_prune_redundancy_heterodimers_parallel_spedup(self):
        rg = RedundantDimerStructures(['3o8o_A-3o8o_B', '3o8o_C-3o8o_D', '3o8o_A-3o8o_D', '3o8o_B-3o8o_C'], 0.2, config)

        t = time.time()
        results1 = set(rg.prune_redundancy(num_workers=1))
        elapsed1 = time.time() - t
        
        t = time.time()
        results2 = set(rg.prune_redundancy(num_workers=4))
        elapsed2 = time.time() - t
        
        expected = set(['3o8o_C-3o8o_D', '3o8o_B-3o8o_C'])
        self.assertEqual(results1, expected)
        self.assertEqual(results2, expected)

        # check speedup
        print(elapsed1, elapsed2)
        self.assertTrue(elapsed2 < elapsed1)
=======
        rg = RedundantDimerStructures(['1xdl_A_1-1xdl_D_0', '1xdl_A_1-1xdl_B_1', '1xdl_C_0-1xdl_D_0', '1xdl_A_1-1xdl_C_0', '4fb8_A_1-4fb8_B_1'], 0.3, config)
        results = set(rg.prune_redundancy(num_workers=3))
        expected = set(['1xdl_A_1-1xdl_D_0', '1xdl_C_0-1xdl_D_0', '1xdl_A_1-1xdl_C_0', '4fb8_A_1-4fb8_B_1'])
        self.assertEqual(results, expected)


#    def test_prune_redundancy_heterodimers_mono(self):
#        rg = RedundantDimerStructures(['3lue_B_1-3lue_M_1', '3lue_C_1-3lue_K_1', '3lue_D_1-3lue_M_1', '3lue_A_1-3lue_K_1'], 0.2, config)
#        results = set(rg.prune_redundancy())
#        expected = set(['3lue_D_1-3lue_M_1', '3lue_A_1-3lue_K_1'])
#        self.assertEqual(results, expected)
#
#
#    def test_prune_redundancy_heterodimers_multi(self):
#        rg = RedundantDimerStructures(['3lue_B_1-3lue_M_1', '3lue_C_1-3lue_K_1', '3lue_D_1-3lue_M_1', '3lue_A_1-3lue_K_1'], 0.2, config)
#        results = set(rg.prune_redundancy(num_workers=4))
#        expected = set(['3lue_D_1-3lue_M_1', '3lue_A_1-3lue_K_1'])
#        self.assertEqual(results, expected)
#
#
#    def test_prune_redundancy_heterodimers_alternate(self):
#        rg = RedundantDimerStructures(['3o8o_A_1-3o8o_B_1', '3o8o_C_1-3o8o_D_1', '3o8o_A_1-3o8o_D_1', '3o8o_B_1-3o8o_C_1'], 0.2, config)
#        results = set(rg.prune_redundancy(num_workers=2))
#        expected = set(['3o8o_C_1-3o8o_D_1', '3o8o_B_1-3o8o_C_1'])
#        self.assertEqual(results, expected)
#
#
#    def test_prune_redundancy_heterodimers_parallel_spedup(self):
#        rg = RedundantDimerStructures(['3o8o_A_1-3o8o_B_1', '3o8o_C_1-3o8o_D_1', '3o8o_A_1-3o8o_D_1', '3o8o_B_1-3o8o_C_1'], 0.2, config)
#
#        t = time.time()
#        results1 = set(rg.prune_redundancy(num_workers=1))
#        elapsed1 = time.time() - t
#        
#        t = time.time()
#        results2 = set(rg.prune_redundancy(num_workers=4))
#        elapsed2 = time.time() - t
#        
#        expected = set(['3o8o_C_1-3o8o_D_1', '3o8o_B_1-3o8o_C_1'])
#        self.assertEqual(results1, expected)
#        self.assertEqual(results2, expected)
#
#        # check speedup
#        print(elapsed1, elapsed2)
#        self.assertTrue(elapsed2 < elapsed1)
>>>>>>> addmodelsplittopdb


if __name__ == '__main__':
    unittest.main()

