import unittest
import sys
import yaml
import os

sys.path.append('..')
import derive_all_possible_heterodimers


class TestDeriveHeterodimers(unittest.TestCase):
    def test_heterodimers_yaml(self):
        infile = 'data/fake_uniparc2others.yaml'
        outfile = 'testout.tmp'
        derive_all_possible_heterodimers.heterodimers(infile, outfile)

        expected = {
            'UPI00000C0390-UPI00001BD6BE': ['1ZZZ_R-1ZZZ_D']
        }

        with open(outfile, 'r') as f:
            result = yaml.safe_load(f)
        os.remove(outfile)

        self.assertTrue(expected == result)


    def test_group_by_assembly(self):
        monomers = [
            ('UP123', 'ABCD', 'A'),
            ('UP123', 'ABCD', 'B'),
            ('UP456', 'ABCD', 'C'),
            ('UP789', 'DEFG', 'Z'),
            ('UP987', 'DEFG', 'R'),
            ('UP987', 'QWER', 'O')
        ]

        expected = [
            [ ('UP123', 'ABCD', 'A'), ('UP123', 'ABCD', 'B'), ('UP456', 'ABCD', 'C') ],
            [ ('UP789', 'DEFG', 'Z'), ('UP987', 'DEFG', 'R') ],
            [ ('UP987', 'QWER', 'O') ]
        ]

        result = derive_all_possible_heterodimers.group_by_assembly(monomers)

        self.assertTrue(expected == result)


    def test_derive_heterodimers_from_groups(self):
        monomers_grouped = [
            [ ('UP123', 'ABCD', 'A'), ('UP123', 'ABCD', 'B'), ('UP456', 'ABCD', 'C') ],
            [ ('UP789', 'DEFG', 'Z'), ('UP987', 'DEFG', 'R'),  ('UP000', 'DEFG', 'A') ],
            [ ('UP987', 'QWER', 'O') ]
        ]

        expected = {
            'UP123-UP456': ['ABCD_A-ABCD_C', 'ABCD_B-ABCD_C'],
            'UP789-UP987': ['DEFG_Z-DEFG_R'],
            'UP000-UP789': ['DEFG_A-DEFG_Z'],
            'UP000-UP987': ['DEFG_A-DEFG_R'] 
        }

        result = derive_all_possible_heterodimers.derive_heterodimers_from_assembly_groups(monomers_grouped)

        self.assertTrue(expected == result)


    def test_chains_as_tuples(self):
        data = {
            'UP123': {
                'uniprot': ['nothin', 'morenothin', 'whatver'],
                'pdb': ['ABCD_A', 'ABCD_B', 'DEFG_Z']
            },
            'UP456': {
                'uniprot': ['othernothing', 'notsomething'],
                'pdb': ['AAAA_D', 'AAAA_R', 'DEFG_V']
            },
            'UP789': {
                'uniprot': ['notaprot'],
                'pdb': ['FFFF_F', 'FFFF_Q']
            }
        }

        expected = [
            ('UP123', 'ABCD', 'A'),
            ('UP123', 'ABCD', 'B'),
            ('UP123', 'DEFG', 'Z'),
            ('UP456', 'AAAA', 'D'),
            ('UP456', 'AAAA', 'R'),
            ('UP456', 'DEFG', 'V'),
            ('UP789', 'FFFF', 'F'),
            ('UP789', 'FFFF', 'Q')
        ]

        result = derive_all_possible_heterodimers.chains_as_tuples(data)

        self.assertEqual(expected, result)






if __name__ == '__main__':
    unittest.main()
