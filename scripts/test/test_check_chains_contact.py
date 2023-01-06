import unittest
import sys

labels_dir = '/scratch/sigbio_project_root/sigbio_project1/jaschwa/HDpred_working/preprocessing/bin/labels'
sys.path.append(labels_dir)
sys.path.append('..')
import labels as label_funcs
import check_chains_contact


class TestCheckChainsContact(unittest.TestCase):
    def test_pos_pairwise(self):
        chain1_pdb = 'data/fake_rcsb/jk/4JK1E.pdb'
        chain2_pdb = 'data/fake_rcsb/jk/4JK1D.pdb'
        result = check_chains_contact._check_chains_contact_pairwise(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertTrue(result)

    def test_neg_pairwise(self):
        chain1_pdb = 'data/fake_rcsb/jk/4JK1H.pdb'
        chain2_pdb = 'data/fake_rcsb/jk/4JK1C.pdb'
        result = check_chains_contact._check_chains_contact_pairwise(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertFalse(result)

    def test_pos_all(self):
        chain1_pdb = 'data/fake_rcsb/jk/4JK1E.pdb'
        chain2_pdb = 'data/fake_rcsb/jk/4JK1D.pdb'
        result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertTrue(result)
        
    def test_neg_all(self):
        chain1_pdb = 'data/fake_rcsb/jk/4JK1H.pdb'
        chain2_pdb = 'data/fake_rcsb/jk/4JK1C.pdb'
        result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
        self.assertFalse(result)

    def test_impossible_negcontact(self):
        chain1_pdb = 'data/fake_rcsb/jk/4JK1A.pdb'  
        chain2_pdb = 'data/fake_rcsb/jk/4JK1Y.pdb'  
        is_impossible = check_chains_contact._check_chains_contact_impossible(chain1_pdb, chain2_pdb, label_funcs, 8)
        self.assertTrue(is_impossible)

    def test_impossible_poscontact(self):
        chain1_pdb = 'data/fake_rcsb/jk/4JK1E.pdb'
        chain2_pdb = 'data/fake_rcsb/jk/4JK1D.pdb'
        is_impossible = check_chains_contact._check_chains_contact_impossible(chain1_pdb, chain2_pdb, label_funcs, 8)
        self.assertFalse(is_impossible)

    
    def test_check_many_neg(self):
        pairs = [('5u8t5.pdb', '5u8t6.pdb'), ('4xmnF.pdb', '4xmnH.pdb'), ('4ij2A.pdb', '4ij2H.pdb'), ('6qleK.pdb', '6qleY.pdb'), ('6sw9R.pdb', '6sw98.pdb')]
        for p in pairs:
            chain1_pdb = 'data/fake_rcsb/' + p[0][1:3] + '/' + p[0]
            chain2_pdb = 'data/fake_rcsb/' + p[1][1:3] + '/' + p[1]
            result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
            self.assertFalse(result)


    def test_check_many_pos(self):
        pairs = [('4uttA.pdb', '4uttD.pdb'), ('4fb8A.pdb', '4fb8B.pdb'), ('5aa5A.pdb', '5aa5B.pdb'), ('3rhtA.pdb', '3rhtD.pdb')]
        for p in pairs:
            chain1_pdb = 'data/fake_rcsb/' + p[0][1:3] + '/' + p[0]
            chain2_pdb = 'data/fake_rcsb/' + p[1][1:3] + '/' + p[1]
            result = check_chains_contact.check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, 8, 10)
            self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()
