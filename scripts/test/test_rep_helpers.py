import os
import sys
sys.path.append('..')

import unittest
import rep_helpers
import string
import numpy as np
import name_pdb

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = f'{test_dir}/data/rep_helpers'
lib_path = os.path.join(data_dir, 'lib')

resolu_path = os.path.join(lib_path, 'resolu.idx')
methods_file = os.path.join(lib_path, 'pdb_entry_type.txt')



class TestRepHelpers(unittest.TestCase):

    def test_dimer_most_like_group_fake_nums(self):
        all_letters = list(string.ascii_lowercase)
        group_letters = ['a', 'c', 'j', 'q']
        A = 0
        C = 2
        J = 9
        Q = 16

        dist_mat = np.random.rand(len(all_letters), len(all_letters))
        dist_mat[A, J] = 0.3
        dist_mat[A, Q] = 0.45
        dist_mat[A, C] = 0.16
        dist_mat[J, C] = 0.99
        dist_mat[Q, C] = 0
        dist_mat[J, Q] = 0.72
        dist_mat[J, A] = 0.3
        dist_mat[Q, A] = 0.45
        dist_mat[C, A] = 0.16
        dist_mat[C, J] = 0.99
        dist_mat[C, Q] = 0
        dist_mat[Q, J] = 0.72

        result = rep_helpers.dimers_of_lowest_distance_to_others(group_letters, dist_mat, all_letters)
        expected = ['c']
        self.assertEqual(result, expected)


    def test_dimer_most_like_group_fake_nums_two_results(self):
        all_letters = list(string.ascii_lowercase)
        group_letters = ['a', 'c', 'j', 'q']
        A = 0
        C = 2
        J = 9
        Q = 16

        dist_mat = np.random.rand(len(all_letters), len(all_letters))
        dist_mat[A, J] = 0.3
        dist_mat[A, Q] = 0.16
        dist_mat[A, C] = 0.16
        dist_mat[J, C] = 0.99
        dist_mat[Q, C] = 0
        dist_mat[J, Q] = 0.99
        dist_mat[J, A] = 0.3
        dist_mat[Q, A] = 0.16
        dist_mat[C, A] = 0.16
        dist_mat[C, J] = 0.99
        dist_mat[C, Q] = 0
        dist_mat[Q, J] = 0.99

        result = rep_helpers.dimers_of_lowest_distance_to_others(group_letters, dist_mat, all_letters)
        expected = ['c', 'q']
        self.assertEqual(result, expected)



    def test_num_residues(self):
        pdbfile = name_pdb.name_pdb_file('1xdl', '1', '1', 'A', lib_path)
        expected = 303
        result = rep_helpers._num_residues(pdbfile)
        self.assertEqual(result, expected)


    def test_dimer_coverage(self):
        dimername = '1xdl-a1-m1-cA_1xdl-a1-m1-cB'
        expected = 300.993355
        result = rep_helpers.dimer_coverage(dimername, lib_path)
        self.assertTrue(abs(result-expected)<0.001)



    def test_get_chain_resolu_xray(self):
        result = rep_helpers._get_chain_resolu('1xdl-a1-m1-cA', resolu_path)
        expected = 3.0
        self.assertEqual(result, expected)


    def test_get_chain_resolu_notxray(self):
        result = rep_helpers._get_chain_resolu('1acz-a1-m1-cA', resolu_path)
        expected = -1.0
        self.assertEqual(result, expected)

    def test_get_chain_resolu_notxray2(self):
        result = rep_helpers._get_chain_resolu('9999-a1-m1-cA', resolu_path)
        expected = -1.0
        self.assertEqual(result, expected)

    def test_get_chain_resolu_fail(self):
        with self.assertRaises(KeyError):
            rep_helpers._get_chain_resolu('9zzz-a1-m1-cA', resolu_path)


    def test_get_dimer_avg_resolu_both_xray(self):
        result = rep_helpers.get_dimer_avg_resolu('1xdl-a1-m1-cA_4zza-a1-m1-cA', resolu_path)
        expected = 2.46171
        self.assertTrue(abs(result-expected)<0.0001)


    def test_get_dimer_avg_resolu_non_xray_right(self):
        result = rep_helpers.get_dimer_avg_resolu('1xdl-a1-m1-cA_1acz-a1-m1-cA', resolu_path)
        expected = -1.0
        self.assertEqual(result, expected)


    def test_get_dimer_avg_resolu_non_xray_left(self):
        result = rep_helpers.get_dimer_avg_resolu('1acz-a1-m1-cA_1xdl-a1-m1-cA', resolu_path)
        expected = -1.0
        self.assertEqual(result, expected)


    def test_get_dimer_avg_resolu_non_xray_both(self):
        result = rep_helpers.get_dimer_avg_resolu('4znf-a2-m4-cZ_1acz-a1-m1-cA', resolu_path)
        expected = -1.0
        self.assertEqual(result, expected)


    def test_dimer_is_xray(self):
        dimer = '1b6q-a1-m1-cA_1b6q-a1-m1-cB'
        self.assertTrue(rep_helpers.dimer_is_xray(dimer, methods_file))


    def test_dimer_em_is_not_xray(self):
        dimer = '7jfo-a1-m1-cA_7jfo-a1-m1-cB'
        self.assertFalse(rep_helpers.dimer_is_xray(dimer, methods_file))


    def test_dimer_nmr_is_not_xray(self):
        dimer = '4znf-a2-m4-cZ_1acz-a1-m1-cA'
        self.assertFalse(rep_helpers.dimer_is_xray(dimer, methods_file))


    def test_dimer_other_is_not_xray(self):
        dimer = '4znn-a2-m4-cZ_4znn-a1-m1-cA'
        self.assertFalse(rep_helpers.dimer_is_xray(dimer, methods_file))


