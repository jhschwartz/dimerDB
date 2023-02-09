import sys
sys.path.append('..')

import os 
import unittest
import yaml
from derive_all_possible_homodimers import *

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
testlib = f'{test_dir}/data/fake_lib'


class TestDeriveHomodimers(unittest.TestCase):
    def test_group_chains_1(self):
        chains = ['1abc_A_1', '1abc_B_1', '7cba_B_1', '1abc_C_1']
        expected = [['1abc_A_1', '1abc_B_1', '1abc_C_1'], ['7cba_B_1']] 
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_group_chains_2(self):
        chains = ['1abd_A_1', '1hbc_B_2', '7cba_B_1', '1abc_C_1']
        expected = [['1abc_C_1'], ['1abd_A_1'], ['1hbc_B_2'], ['7cba_B_1']] #sorted
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_group_chains_3(self):
        chains = ['1abc_A_1', '1abc_B_1', '1abc_C_2', '1abc_D_1']
        expected = [['1abc_A_1', '1abc_B_1', '1abc_C_2', '1abc_D_1']] 
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_derive_homodimers_from_groups_1(self):
        groups = [['1abc_A_2', '1abc_B_1', '1abc_C_1'], ['7cba_B_1']]
        expected = [('1abc_A_2', '1abc_B_1'), ('1abc_A_2', '1abc_C_1'), ('1abc_B_1', '1abc_C_1')]
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)
    
    def test_derive_homodimers_from_groups_2(self):
        groups = [['1abc_C_1'], ['1abd_A_1'], ['1hbc_B_3'], ['7cba_B_16']]
        expected = []
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)

    def test_derive_homodimers_from_groups_3(self):
        groups = [['1abc_A_1', '1abc_B_1', '1abc_C_1', '1abc_D_1']]
        expected = [('1abc_A_1', '1abc_B_1'), ('1abc_A_1', '1abc_C_1'), ('1abc_A_1', '1abc_D_1'), ('1abc_B_1', '1abc_C_1'), ('1abc_B_1', '1abc_D_1'), ('1abc_C_1', '1abc_D_1')]
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)

    def test_homodimers_yaml(self):
        infile = 'data/fake_uniparc2others.yaml'
        outfile = 'testout.tmp'
        homodimers(infile, outfile, testlib)

        expected = {
            'UPI00000C0390': ['1on3_A_1-1on3_B_1', '1on9_A_1-1on9_B_1', '1on9_A_1-1on9_C_1', '1on9_B_1-1on9_C_1'],
            'UPI00001BD6BE': ['6k3g_B_1-6k3g_B_2']
        }
        
        with open(outfile, 'r') as f:
            result = yaml.safe_load(f)
        os.remove(outfile)

        self.assertTrue(expected == result)


if __name__ == '__main__':
    unittest.main()
        


