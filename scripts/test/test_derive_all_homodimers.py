import sys
sys.path.append('..')

import os 
import unittest
import yaml
from derive_all_homodimers import *

class TestDeriveHomodimers(unittest.TestCase):
    def test_group_chains_1(self):
        chains = ['1abc_A', '1abc_B', '7cba_B', '1abc_C']
        expected = [['1abc_A', '1abc_B', '1abc_C'], ['7cba_B']] 
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_group_chains_2(self):
        chains = ['1abd_A', '1hbc_B', '7cba_B', '1abc_C']
        expected = [['1abc_C'], ['1abd_A'], ['1hbc_B'], ['7cba_B']] #sorted
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_group_chains_3(self):
        chains = ['1abc_A', '1abc_B', '1abc_C', '1abc_D']
        expected = [['1abc_A', '1abc_B', '1abc_C', '1abc_D']] 
        result = group_chains(chains)
        self.assertEqual(expected, result)

    def test_derive_homodimers_from_groups_1(self):
        groups = [['1abc_A', '1abc_B', '1abc_C'], ['7cba_B']]
        expected = [('1abc_A', '1abc_B'), ('1abc_A', '1abc_C'), ('1abc_B', '1abc_C')]
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)
    
    def test_derive_homodimers_from_groups_2(self):
        groups = [['1abc_C'], ['1abd_A'], ['1hbc_B'], ['7cba_B']]
        expected = []
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)

    def test_derive_homodimers_from_groups_3(self):
        groups = [['1abc_A', '1abc_B', '1abc_C', '1abc_D']]
        expected = [('1abc_A', '1abc_B'), ('1abc_A', '1abc_C'), ('1abc_A', '1abc_D'), ('1abc_B', '1abc_C'), ('1abc_B', '1abc_D'), ('1abc_C', '1abc_D')]
        result = derive_homodimers_from_groups(groups)
        self.assertEqual(expected, result)

    def test_homodimers_yaml(self):
        infile = 'data/fake_uniparc2others.yaml'
        outfile = 'testout.tmp'
        homodimers(infile, outfile)

        expected = {
            'UPI00000C0390': ['1ON3_A-1ON3_B', '1ON9_A-1ON9_B', '1ON9_A-1ON9_C', '1ON9_B-1ON9_C'],
            'UPI00001BD6BE': ['6K3G_B-6K3G_FAKE']
        }
        
        with open(outfile, 'r') as f:
            result = yaml.safe_load(f)
        os.remove(outfile)

        self.assertTrue(expected == result)


if __name__ == '__main__':
    unittest.main()
        


