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
