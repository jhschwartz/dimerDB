import sys
sys.path.append('..')
import uniparc_to_uniprot_and_pdb
import os
import yaml
import unittest



class Uniparc2othersTestCase(unittest.TestCase):
    def test_uniparc2others(self):
        infile = 'data/idmapping_selected.tab.fake.gz'
        outfile = 'testout.yaml.tmp'
        uniparc_to_uniprot_and_pdb.make_uniparc2others(infile, outfile)

        expected = {
            'UPI00000C0390': {
                'uniprot': ['Q8GBW6'],
                'pdb': ['1ON3_A', '1ON3_B', '1ON3_C', '1ON3_D', '1ON3_E', '1ON3_F', '1ON9_A', '1ON9_B', '1ON9_C', '1ON9_D', '1ON9_E', '1ON9_F']
            },
            'UPI0000124DBC': {
                'pdb': ['2E9Q_A', '2EVX_A'],
                'uniprot': ['ABC123', 'P13744']
            },
            'UPI00001BD6BE': {
                'pdb': ['1234_Z', '5555_R', '6K3G_B', '6KJ5_A'],
                'uniprot': ['Q6V4H0', 'XYZ987']
            }
        }

        with open(outfile, 'r') as f:
            result = yaml.safe_load(f)
        os.remove(outfile)

        self.assertTrue(expected == result)


    def test_extract_chains_many(self):
        fakechains = 'abcd:A; efgh:A; zxvd:B; '
        result = uniparc_to_uniprot_and_pdb._extract_chains(fakechains)
        expected = {'abcd_A', 'efgh_A', 'zxvd_B'}
        self.assertEqual(result, expected)


    def test_extract_chains_one(self):
        fakechains = 'abcd:A; '
        result = uniparc_to_uniprot_and_pdb._extract_chains(fakechains)
        expected = {'abcd_A'}
        self.assertEqual(result, expected)


    def test_extract_chains_zero(self):
        fakechains = ''
        result = uniparc_to_uniprot_and_pdb._extract_chains(fakechains)
        expected = set()
        self.assertEqual(result, expected)   


    def test_extract_ids(self):
        fakeline = 'Q6V4H0	10HGO_CATRO			75324836	6K3G:B; 6KJ5:A	GO:0102311; GO:0008270; GO:0071704	UniRef100_Q6V4H0	UniRef90_Q6V4H0	UniRef50_A0A1B1FHP3	UPI00001BD6BE		4058				AY352047	AAQ55962.1				32181958'
        result = uniparc_to_uniprot_and_pdb._extract_ids(fakeline)
        expected = ('Q6V4H0', {'6K3G_B', '6KJ5_A'}, 'UPI00001BD6BE')
        self.assertEqual(result, expected)


    def test_sort_entries(self):
        fakeunsorted = {
            'a': {
                'uniprot': { 'z', 'a', 'g'},
                'pdb': {'5', '1', '19'}
            },
            'b': {
                'uniprot': { 'R', 'ZZZ', 'AAA'},
                'pdb': {'123', '456', '0000000'}
            }
        }
        uniparc_to_uniprot_and_pdb._sort_entries(fakeunsorted)
        expected = {
            'a': {
                'uniprot': ['a', 'g', 'z'],
                'pdb': ['1', '19', '5']
            },
            'b': {
                'uniprot': ['AAA', 'R', 'ZZZ'],
                'pdb': ['0000000', '123', '456']
            }
        }
        self.assertTrue(fakeunsorted == expected)


    def test_clear_entries_with_no_chains(self):
        fake = {
            'a': {
                'uniprot': {'a', 'g', 'z'},
                'pdb': {'1', '5', '19'}
            },
            'b': {
                'uniprot': {'AAA', 'R', 'ZZZ'},
                'pdb': {}
            }
        }
        result = uniparc_to_uniprot_and_pdb._clear_entries_with_no_chains(fake)
        expected = {
            'a': {
                'uniprot': {'a', 'g', 'z'},
                'pdb': {'1', '5', '19'}
            }
        }
        self.assertEqual(expected, result)
        

if __name__ == '__main__':
    unittest.main()
