import unittest
import sys

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()

test_data = f'{test_dir}/data/check_chains_contact'

labels_dir = f'{test_dir}/../../bin/labels'
sys.path.append(labels_dir)
sys.path.append('..')
import labels as label_funcs
import check_chains_contact


class TestCheckChainsContact(unittest.TestCase):
    def test_pos_pairwise(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-E-1.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-D-1.pdb'
        result = check_chains_contact._check_chains_contact_pairwise(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertTrue(result)

    def test_neg_pairwise(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-H-1.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-C-1.pdb'
        result = check_chains_contact._check_chains_contact_pairwise(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertFalse(result)

    def test_pos_all(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-E-1.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-D-1.pdb'
        result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertTrue(result)
        
    def test_neg_all(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-H-1.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-C-1.pdb'
        result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertFalse(result)

    def test_impossible_negcontact(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-A-1.pdb'  
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-Y-1.pdb'  
        is_impossible = check_chains_contact._check_chains_contact_impossible(chain1_pdb, chain2_pdb, label_funcs, 8)
        self.assertTrue(is_impossible)

    def test_impossible_poscontact(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-E-1.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1/chains/4jk1-D-1.pdb'
        is_impossible = check_chains_contact._check_chains_contact_impossible(chain1_pdb, chain2_pdb, label_funcs, 8)
        self.assertFalse(is_impossible)

    
    def test_check_many_neg(self):
        pairs = [('5u8t-5-1.pdb', '5u8t-6-1.pdb'), ('4xmn-F-1.pdb', '4xmn-H-1.pdb'), ('4ij2-A-1.pdb', '4ij2-H-1.pdb'), ('6qle-K-1.pdb', '6qle-Y-1.pdb'), ('6sw9-R-1.pdb', '6sw9-8-1.pdb')]
        for p in pairs:
            pdb1 = p[0].split('-')[0]
            pdb2 = p[1].split('-')[0]
            div1 = pdb1[1:3]
            div2 = pdb2[1:3]
            chain1_pdb = f'{test_data}/lib/rcsb/{div1}/{pdb1}/chains/{p[0]}' 
            chain2_pdb = f'{test_data}/lib/rcsb/{div2}/{pdb2}/chains/{p[1]}'
            result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
            self.assertFalse(result)


    def test_check_many_pos(self):
        pairs = [('4utt-A-1.pdb', '4utt-D-1.pdb'), ('4fb8-A-1.pdb', '4fb8-B-1.pdb'), ('5aa5-A-1.pdb', '5aa5-B-1.pdb'), ('3rht-A-1.pdb', '3rht-D-1.pdb')]
        for p in pairs:
            pdb1 = p[0].split('-')[0]
            pdb2 = p[1].split('-')[0]
            div1 = pdb1[1:3]
            div2 = pdb2[1:3]
            chain1_pdb = f'{test_data}/lib/rcsb/{div1}/{pdb1}/chains/{p[0]}' 
            chain2_pdb = f'{test_data}/lib/rcsb/{div2}/{pdb2}/chains/{p[1]}'
            result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
            self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()
