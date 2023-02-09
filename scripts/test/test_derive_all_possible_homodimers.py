import sys
sys.path.append('..')

import os 
import unittest
import yaml
import tempfile
from derive_all_possible_homodimers import *

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
test_data = f'{test_dir}/data/derive_all_possible_homodimers'
test_lib = f'{test_data}/lib'


class TestDeriveHomodimers(unittest.TestCase):
    def test_group_chains_1(self):
        chains = ['1abc_a1_m1_cA', '1abc_a1_m1_cB', '7cba_a1_m1_cB', '1abc_a1_m1_cC']
        expected = [['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cC'], ['7cba_a1_m1_cB']] 
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_group_chains_2(self):
        chains = ['1abd_a1_m1_cA', '1hbc_a1_m2_cB', '7cba_a1_m1_cB', '1abc_a1_m1_cC']
        expected = [['1abc_a1_m1_cC'], ['1abd_a1_m1_cA'], ['1hbc_a1_m2_cB'], ['7cba_a1_m1_cB']] #sorted
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_group_chains_3(self):
        chains = ['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cD', '1abc_a1_m2_cC', ]
        expected = [['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cD', '1abc_a1_m2_cC', ]] 
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_derive_homodimers_from_groups_1(self):
        groups = [['1abc_a1_m2_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cC'], ['7cba_a1_m1_cB']]
        expected = [('1abc_a1_m2_cA', '1abc_a1_m1_cB'), ('1abc_a1_m2_cA', '1abc_a1_m1_cC'), ('1abc_a1_m1_cB', '1abc_a1_m1_cC')]
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)
    
    def test_derive_homodimers_from_groups_2(self):
        groups = [['1abc_a1_m1_cC'], ['1abd_a1_m1_cA'], ['1hbc_a1_m3_cB'], ['7cba_a1_m1_cB6']]
        expected = []
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)

    def test_derive_homodimers_from_groups_3(self):
        groups = [['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a1_m1_cC', '1abc_a1_m1_cD']]
        expected = [('1abc_a1_m1_cA', '1abc_a1_m1_cB'), ('1abc_a1_m1_cA', '1abc_a1_m1_cC'), ('1abc_a1_m1_cA', '1abc_a1_m1_cD'), ('1abc_a1_m1_cB', '1abc_a1_m1_cC'), ('1abc_a1_m1_cB', '1abc_a1_m1_cD'), ('1abc_a1_m1_cC', '1abc_a1_m1_cD')]
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)

    def test_homodimers_yaml(self):
        infile = f'{test_data}/fake_uniparc2others.yaml'
        outfile = 'testout.tmp'
        with tempfile.TemporaryDirectory() as td:
            outfile = f'{td}/out.yaml'
            homodimers(infile, outfile, test_lib)

            expected = {
                'UPI00000C0390': ['1on3_a1_m1_cA-1on3_a1_m1_cB', '1on9_a1_m1_cA-1on9_a1_m1_cB', '1on9_a1_m1_cA-1on9_a1_m1_cC', '1on9_a1_m1_cB-1on9_a1_m1_cC'],
                'UPI00001BD6BE': ['6k3g_a1_m1_cB-6k3g_a1_m2_cB', '6k3g_a1_m1_cB-6k3g_a2_m2_cB', '6k3g_a1_m2_cB-6k3g_a2_m2_cB']
            }
            
            with open(outfile, 'r') as f:
                result = yaml.safe_load(f)

            self.assertTrue(expected == result)



    def test_expand_chains_across_assemblies_models(self):
        chains = ['1on3_A', '1on3_B', '6k3g_B']
        result = expand_chains_across_assemblies_models(chains, test_lib)
        expected = ['1on3_a1_m1_cA', '1on3_a1_m1_cB', '6k3g_a1_m1_cB', '6k3g_a1_m2_cB', '6k3g_a2_m2_cB']
        self.assertEqual(set(result), set(expected))


if __name__ == '__main__':
    unittest.main()
        


