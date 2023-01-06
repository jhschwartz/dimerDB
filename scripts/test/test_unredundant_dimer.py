import unittest
import sys
import os

sys.path.append('..')
from unredundant import RedundantDimers


class TestRedundantDimers(unittest.TestCase):
    def test_init(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, '/tmp')
        self.assertEqual(rg.things[0], 'placeholder')
        self.assertEqual(rg.things[1], 'nothing')
        self.assertEqual(rg.threshold, 10)
        self.assertEqual(rg.rcsb_base_path, '/tmp')


    def test_distance_unsimilar(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, 'data/fake_rcsb')
        dimer1name = '1xdl_A-1xdl_B'
        dimer2name = '1xdl_A-1xdl_D'
        expected = 0.50732 
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_distance_similar(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, 'data/fake_rcsb')
        dimer1name = '1xdl_A-1xdl_B'
        dimer2name = '1xdl_C-1xdl_D'
        expected = 0.97924
        result = rg.distance(dimer1name, dimer2name)
        self.assertEqual(result, expected)


    def test_num_residues(self):
        pdbfile = 'data/fake_rcsb/xd/1xdlA.pdb'
        expected = 303
        result = RedundantDimers._num_residues(pdbfile)
        self.assertEqual(result, expected)


    def test_dimer_coverage(self):
        rg = RedundantDimers(['placeholder', 'nothing'], 10, 'data/fake_rcsb')
        dimername = '1xdl_A-1xdl_B'
        expected = 300.993355
        result = rg._dimer_coverage(dimername)
        self.assertTrue(abs(result-expected)<0.0001)


    def test_representative_first_criterion(self):
        assert False
        dimers = ['1-2', '4-5', 'abc-123']
        representative = RedundantDimers.representative(dimers)
        expected = '1-2'
        self.assertEqual(result, expected)
    
    def test_representative_second_criterion(self):
        assert False


    def test_representative_third_criterion(self):
        assert False


    def test_tmp_assembly_file(self):
        assert False
        rg = RedundantDimers(['placeholder', 'nothing'], 10, 'data/fake_rcsb')
        with rg._tmp_assembly_file('1-2') as f:
            self.assertTrue(os.path.exists(f))
            self.assertEqual(RedundantDimers._num_residues(f), 1234)


    def test_filenames(self):
        assert False
        rg = RedundantDimers(['placeholder', 'nothing'], 10, 'data/fake_rcsb')
        dimername = '123'
        f0, f1 = rg._filenames(dimername)
        self.assertEqual(f0, 'data/fake_rcsb/??.pdb')
        self.assertEqual(f1, 'data/fake_rcsb/??.pdb')
        self.assertRaises(rg._filenames('notarealdimer.notpdb'), FileNotFoundError)



    def test_prune_redundancy_homodimers_1(self):
        assert False

    def test_prune_redundancy_homodimers_2(self):
        assert False

    def test_prune_redundancy_heterodimers_1(self):
        assert False

    def test_prune_redundancy_heterodimers_2(self):
        assert False
