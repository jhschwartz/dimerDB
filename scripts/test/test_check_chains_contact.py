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
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cE.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cD.pdb'
        result = check_chains_contact._check_chains_contact_pairwise(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertTrue(result)

    def test_neg_pairwise(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cH.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cC.pdb'
        result = check_chains_contact._check_chains_contact_pairwise(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertFalse(result)

    def test_pos_all(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cE.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cD.pdb'
        result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertTrue(result)
        
    def test_neg_all(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cH.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cC.pdb'
        result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertFalse(result)

    def test_impossible_negcontact(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cA.pdb'  
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cY.pdb'  
        is_impossible = check_chains_contact._check_chains_contact_impossible(chain1_pdb, chain2_pdb, label_funcs, 8)
        self.assertTrue(is_impossible)

    def test_impossible_poscontact(self):
        chain1_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cE.pdb'
        chain2_pdb = f'{test_data}/lib/rcsb/jk/4jk1-a1-m1-cD.pdb'
        is_impossible = check_chains_contact._check_chains_contact_impossible(chain1_pdb, chain2_pdb, label_funcs, 8)
        self.assertFalse(is_impossible)

    
    def test_check_many_neg(self):
        pairs = [('5u8t-a1-m1-c5.pdb', '5u8t-a1-m1-c6.pdb'), ('4xmn-a1-m1-cF.pdb', '4xmn-a1-m1-cH.pdb'), ('4ij2-a1-m1-cA.pdb', '4ij2-a1-m1-cH.pdb'), ('6qle-a1-m1-cK.pdb', '6qle-a1-m1-cY.pdb'), ('6sw9-a1-m1-cR.pdb', '6sw9-a1-m1-c8.pdb')]
        for p in pairs:
            pdb1 = p[0].split('-')[0]
            pdb2 = p[1].split('-')[0]
            div1 = pdb1[1:3]
            div2 = pdb2[1:3]
            chain1_pdb = f'{test_data}/lib/rcsb/{div1}/{p[0]}' 
            chain2_pdb = f'{test_data}/lib/rcsb/{div2}/{p[1]}'
            result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
            self.assertFalse(result)


    def test_check_many_pos(self):
        pairs = [('4utt-a1-m1-cA.pdb', '4utt-a1-m1-cD.pdb'), ('4fb8-a1-m1-cA.pdb', '4fb8-a1-m1-cB.pdb'), ('5aa5-a1-m1-cA.pdb', '5aa5-a1-m1-cB.pdb'), ('3rht-a1-m1-cA.pdb', '3rht-a1-m1-cD.pdb')]
        for p in pairs:
            pdb1 = p[0].split('-')[0]
            pdb2 = p[1].split('-')[0]
            div1 = pdb1[1:3]
            div2 = pdb2[1:3]
            chain1_pdb = f'{test_data}/lib/rcsb/{div1}/{p[0]}' 
            chain2_pdb = f'{test_data}/lib/rcsb/{div2}/{p[1]}'
            result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
            self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()
