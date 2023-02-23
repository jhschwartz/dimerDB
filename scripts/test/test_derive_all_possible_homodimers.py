import sys
sys.path.append('..')

import unittest
import pickle
import tempfile
from derive_all_possible_homodimers import \
        _group_chains, _derive_homodimers_from_groups, derive_homodimers


import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
test_data = f'{test_dir}/data/derive_all_possible_homodimers'
test_lib = f'{test_data}/lib'




class TestDeriveHomodimers(unittest.TestCase):
    def test_group_chains_1(self):
        chains = ['1abc_a1_m1_cA', '1abc_a1_m1_cB', '7cba_a1_m1_cB', '1abc_a1_m1_cC']
        expected = [['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cC'], ['7cba_a1_m1_cB']] 
        result = _group_chains(chains)
        self.assertEqual(expected, result)

    def test_group_chains_2(self):
        chains = ['1abd_a1_m1_cA', '1hbc_a1_m2_cB', '7cba_a1_m1_cB', '1abc_a1_m1_cC']
        expected = [['1abc_a1_m1_cC'], ['1abd_a1_m1_cA'], ['1hbc_a1_m2_cB'], ['7cba_a1_m1_cB']] #sorted
        result = _group_chains(chains)
        self.assertEqual(expected, result)

    def test_group_chains_3(self):
        chains = ['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cD', '1abc_a1_m2_cC', '1abc_a7_m2_cT','1abc_a7_m2_cR']
        expected = [['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cD', '1abc_a1_m2_cC'], ['1abc_a7_m2_cR', '1abc_a7_m2_cT']] 
        result = _group_chains(chains)
        self.assertEqual(expected, result)

    def test_derive_homodimers_from_groups_1(self):
        groups = [['1abc_a1_m2_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cC'], ['7cba_a1_m1_cB']]
        expected = [('1abc_a1_m2_cA', '1abc_a1_m1_cB'), ('1abc_a1_m2_cA', '1abc_a1_m1_cC'), ('1abc_a1_m1_cB', '1abc_a1_m1_cC')]
        result = _derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)
    
    def test_derive_homodimers_from_groups_2(self):
        groups = [['1abc_a1_m1_cC'], ['1abd_a1_m1_cA'], ['1hbc_a1_m3_cB'], ['7cba_a1_m1_cB6']]
        expected = []
        result = _derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)

    def test_derive_homodimers_from_groups_3(self):
        groups = [['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cC', '1abc_a1_m1_cD']]
        expected = [('1abc_a1_m1_cA', '1abc_a1_m1_cB'), ('1abc_a1_m1_cA', '1abc_a1_m1_cC'), ('1abc_a1_m1_cA', '1abc_a1_m1_cD'), ('1abc_a1_m1_cB', '1abc_a1_m1_cC'), ('1abc_a1_m1_cB', '1abc_a1_m1_cD'), ('1abc_a1_m1_cC', '1abc_a1_m1_cD')]
        result = _derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)

    def test_homodimers_pkl(self):
        infile = f'{test_data}/expanded_cleaned_uniparc2others.pkl'
        outfile = 'testout.tmp'
        with tempfile.TemporaryDirectory() as td:
            outfile = f'{td}/out.pkl'
            derive_homodimers(infile, outfile)

            expected = {
                'UPI00000C0390': ['1on3_a1_m1_cA-1on3_a1_m1_cB', '1on9_a1_m1_cA-1on9_a1_m1_cB', '1on9_a1_m1_cA-1on9_a1_m1_cC', '1on9_a1_m1_cB-1on9_a1_m1_cC'],
                'UPI00001BD6BE': ['6k3g_a1_m1_cB-6k3g_a1_m2_cB', '6k3g_a2_m14_cN-6k3g_a2_m2_cB']
            }
            
            with open(outfile, 'rb') as f:
                result = pickle.load(f)

            self.assertTrue(expected == result)






if __name__ == '__main__':
    unittest.main()
        


