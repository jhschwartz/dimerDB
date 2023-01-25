#!/nfs/turbo/umms-petefred/jaschwa/.conda/HDPRED/bin/python
import sys, os, subprocess
import unittest
import numpy as np
import pickle as pkl
from tools.alt_calcs import calc_dihedral, dist, calc_angle_b

sys.path.append('..')
import labels as label_funcs




def check_angles(testcase, nums_1, nums_2, res_dic_1, res_dic_2, labels):
    for numA in nums_1:
        for numB in nums_2:
            ca_A = res_dic_1[f'{numA}_CA'].flatten()
            ca_B = res_dic_2[f'{numB}_CA'].flatten()
            testcase.assertAlmostEqual(dist(ca_A, ca_B), labels['ca_dis'][int(numA), int(numB)], delta=1e-5)
            #try:
            #    testcase.assertAlmostEqual(dist(ca_A, ca_B), labels['ca_dis'][int(numA), int(numB)], delta=1e-5)
            #except:
            #    pass
        
            if f'{numA}_CB' in res_dic_1 and f'{numB}_CB' in res_dic_2:
                # phi
                a = ca_A
                b = res_dic_1[f'{numA}_CB'].flatten()
                c = res_dic_2[f'{numB}_CB'].flatten()
                phi = calc_angle_b(a, b, c)
                testcase.assertAlmostEqual(phi, labels['phi'][int(numA), int(numB)], delta=1e-5)
               
                # omega
                d = res_dic_2[f'{numB}_CA'].flatten()
                omega = calc_dihedral(a, b, c, d)
                testcase.assertAlmostEqual(omega, labels['omega'][int(numA), int(numB)], delta=1e-5)
                
                # theta
                a = res_dic_1[f'{numA}_N'].flatten()
                b = res_dic_1[f'{numA}_CA'].flatten()
                c = res_dic_1[f'{numA}_CB'].flatten()
                d = res_dic_2[f'{numB}_CB'].flatten()
                theta = calc_dihedral(a, b, c, d)
                testcase.assertAlmostEqual(theta, labels['theta'][int(numA), int(numB)], delta=1e-5)

            else:
                testcase.assertTrue(np.isnan(labels['phi'][int(numA), int(numB)]))
                testcase.assertTrue(np.isnan(labels['omega'][int(numA), int(numB)]))
                testcase.assertTrue(np.isnan(labels['theta'][int(numA), int(numB)]))



class TestLabelGen(unittest.TestCase):
    def test_extract_labels_fxn_monomer(self):
        res_dic, nums, _ = label_funcs.read_chain('data/exampleA.pdb')
        labels = label_funcs.extract_labels(res_dic, nums)
        check_angles(self, nums, nums, res_dic, res_dic, labels)

    def test_extract_labels_fxn_dimer(self):
        res_dic_1, nums_1, _ = label_funcs.read_chain('data/exampleA.pdb') 
        res_dic_2, nums_2, _ = label_funcs.read_chain('data/exampleB.pdb')

        
        intrachain_labels_1, intrachain_labels_2, interchain_labels, L1, L2 = label_funcs.extract_labels_dimer(res_dic_1, nums_1, res_dic_2, nums_2)
       
        self.assertEqual(L1, len(nums_1))
        self.assertEqual(L2, len(nums_2))

        # check intrachain labels
        check_angles(self, nums_1, nums_1, res_dic_1, res_dic_1, intrachain_labels_1)
        check_angles(self, nums_2, nums_2, res_dic_2, res_dic_2, intrachain_labels_2)

        # check interchain labels
        check_angles(self, nums_1, nums_2, res_dic_1, res_dic_2, interchain_labels)

    def test_dimer_label_gen_exe(self):
        subprocess.check_call('../gen_labels_dimer.py data/exampleA.pdb data/seqA.txt data/exampleB.pdb data/seqB.txt testout.tmp', shell=True)
        with open('testout.tmp' , 'rb') as f:
            data = pkl.load(f)
        os.remove('testout.tmp')
            
        intrachain_labels_1 = data['intrachain_labels_1']
        intrachain_labels_2 = data['intrachain_labels_2']
        interchain_labels = data['interchain_labels']
        L1 = data['L1']
        L2 = data['L2']
        seqA = data['seq1']
        seqB = data['seq2']
            
        res_dic_1_fxn, nums_1_fxn, seq1_fxn = label_funcs.read_chain('data/exampleA.pdb') 
        res_dic_2_fxn, nums_2_fxn, seq2_fxn = label_funcs.read_chain('data/exampleB.pdb')

        intrachain_labels_1_fxn, intrachain_labels_2_fxn, interchain_labels_fxn, L1_fxn, L2_fxn = label_funcs.extract_labels_dimer(res_dic_1_fxn, nums_1_fxn, res_dic_2_fxn, nums_2_fxn)
       
        self.assertEqual(L1, L1_fxn)
        self.assertEqual(L2, L2_fxn)
        self.assertEqual(seqA, seq1_fxn)
        self.assertEqual(seqB, seq2_fxn)

        check_angles(self, nums_1_fxn, nums_1_fxn, res_dic_1_fxn, res_dic_1_fxn, intrachain_labels_1)
        check_angles(self, nums_2_fxn, nums_2_fxn, res_dic_2_fxn, res_dic_2_fxn, intrachain_labels_2)
        check_angles(self, nums_1_fxn, nums_2_fxn, res_dic_1_fxn, res_dic_2_fxn, interchain_labels)
    

    def test_monomer_label_gen_exe(self):
        subprocess.check_call('../gen_labels_monomer.py data/exampleA.pdb data/seqA.txt testout.tmp', shell=True)
        with open('testout.tmp' , 'rb') as f:
            data = pkl.load(f)
        os.remove('testout.tmp')
            
        intrachain_labels = data['intrachain_labels']
        L = data['L']
        seqA = data['seq']
            
        res_dic_fxn, nums_fxn, seq1_fxn = label_funcs.read_chain('data/exampleA.pdb') 
        intrachain_labels_fxn = label_funcs.extract_labels(res_dic_fxn, nums_fxn)

        self.assertEqual(seqA, seq1_fxn)
        check_angles(self, nums_fxn, nums_fxn, res_dic_fxn, res_dic_fxn, intrachain_labels)


if __name__ == '__main__':
    unittest.main()
    #t = TestLabelGen()
    #t.test_extract_labels_monomer()
