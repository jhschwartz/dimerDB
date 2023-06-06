import os
import sys
sys.path.append('..')

import unittest
import wrap_tmscore


import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = f'{test_dir}/data/wrap_tmscore'
lib_path = os.path.join(data_dir, 'lib')

exe = '../../bin/USalign/USalign'

_4m1uE = os.path.join(lib_path, 'rcsb', 'm1', '4m1u-a1-m1-cE.pdb') 
_4m1uF = os.path.join(lib_path, 'rcsb', 'm1', '4m1u-a1-m1-cF.pdb')
_4m1uB = os.path.join(lib_path, 'rcsb', 'm1', '4m1u-a1-m1-cB.pdb')
_4m1uC = os.path.join(lib_path, 'rcsb', 'm1', '4m1u-a1-m1-cC.pdb')

_2gqqA1 = os.path.join(lib_path, 'rcsb', 'gq', '2gqq-a1-m1-cA.pdb') 
_2gqqA2 = os.path.join(lib_path, 'rcsb', 'gq', '2gqq-a1-m2-cA.pdb')
_2gqqB2 = os.path.join(lib_path, 'rcsb', 'gq', '2gqq-a1-m2-cB.pdb')
_2gqqC1 = os.path.join(lib_path, 'rcsb', 'gq', '2gqq-a1-m1-cC.pdb')
_2gqqC2 = os.path.join(lib_path, 'rcsb', 'gq', '2gqq-a1-m2-cC.pdb')
_2gqqD1 = os.path.join(lib_path, 'rcsb', 'gq', '2gqq-a1-m1-cD.pdb')
_2gqqD2 = os.path.join(lib_path, 'rcsb', 'gq', '2gqq-a1-m2-cD.pdb')



class TestWrapTMscore(unittest.TestCase):

    def test_compare_dimers(self):
        dimer1 = (_4m1uE, _4m1uF)
        dimer2 = (_4m1uB, _4m1uC)
        result = wrap_tmscore.calculate_dimers_TM_score(dimer1, dimer2, exe)
        expected = (0.98, 0.98)
        self.assertEqual(result, expected)


    def test_compare_many_dimers_linear(self):
        dimer_pairs = [ ('4m1u-a1-m1-cE_4m1u-a1-m1-cF', '4m1u-a1-m1-cB_4m1u-a1-m1-cC'),
                        ('2gqq-a1-m1-cA_2gqq-a1-m2-cB', '2gqq-a1-m2-cB_2gqq-a1-m1-cB'), 
                        ('5zfv-a1-m1-cC_5zfv-a1-m1-cD', '5zfv-a1-m1-cB_5zfv-a1-m1-cC') ]
        result = wrap_tmscore.calculate_many_dimers_TM_score(dimer_pairs, exe, lib_path, 1)
        expected = [ (0.98, 0.98), (0.52, 0.53), (0.95, 0.93) ]
        self.assertEqual(result, expected)


    def test_compare_many_dimers_parallel(self):
        dimer_pairs = [ ('4m1u-a1-m1-cE_4m1u-a1-m1-cF', '4m1u-a1-m1-cB_4m1u-a1-m1-cC'),
                        ('2gqq-a1-m1-cA_2gqq-a1-m2-cB', '2gqq-a1-m2-cB_2gqq-a1-m1-cB'), 
                        ('5zfv-a1-m1-cC_5zfv-a1-m1-cD', '5zfv-a1-m1-cB_5zfv-a1-m1-cC') ]
        result = wrap_tmscore.calculate_many_dimers_TM_score(dimer_pairs, exe, lib_path, 3)
        expected = [ (0.98, 0.98), (0.52, 0.53), (0.95, 0.93) ]
        
        self.assertEqual(result, expected)


    def test_not_fail_on_rechain_assembly2_case1(self):
        dimerA1B2 = (_2gqqA1, _2gqqB2)	
        dimerB2D2 = (_2gqqB2, _2gqqD2)	
        
        result = wrap_tmscore.calculate_dimers_TM_score(dimerA1B2, dimerB2D2, exe)
        expected = ( 0.96, 0.96)
        self.assertEqual(result, expected)


    def test_not_fail_on_rechain_assembly2_case2(self):
        dimerB2D2 = (_2gqqB2, _2gqqD2)	
        dimerA2C2 = (_2gqqA2, _2gqqC2)	
        
        result = wrap_tmscore.calculate_dimers_TM_score(dimerA2C2, dimerB2D2, exe)
        expected = ( 0.95, 0.97)
        self.assertEqual(result, expected)
    

    def test_not_fail_on_rechain_assembly2_case3(self):
        dimerC1D1 = (_2gqqC1, _2gqqD1)	
        dimerC2D2 = (_2gqqC2, _2gqqD2)	
        
        result = wrap_tmscore.calculate_dimers_TM_score(dimerC1D1, dimerC2D2, exe)
        expected = ( 1.0, 1.0)
        self.assertEqual(result, expected)

