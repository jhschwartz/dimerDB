import os
import sys
sys.path.append('..')

import unittest
from align_tools import calc_nwalign_glocal

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = f'{test_dir}/data/align_tools'


pdbshorter = os.path.join(data_dir, '4fb8-a1-m1-cB.pdb')
pdblonger = os.path.join(data_dir, '4fb8-a1-m1-cA.pdb')

exe = '../../bin/USalign/NWalign'


class TestAlignTools(unittest.TestCase):

    def test_nwalign_glocal_refshorter(self):

        expected = 0.984
        result1 = calc_nwalign_glocal(USnw=exe, pdb1=pdbshorter, pdb2=pdblonger, refshorter=True, reflonger=False)
        result2 = calc_nwalign_glocal(USnw=exe, pdb2=pdbshorter, pdb1=pdblonger, refshorter=True, reflonger=False)

        self.assertEqual(result1, expected)
        self.assertEqual(result2, expected)
    
    
    def test_nwalign_glocal_reflonger(self):

        expected = 0.966
        result1 = calc_nwalign_glocal(USnw=exe, pdb1=pdbshorter, pdb2=pdblonger, refshorter=False, reflonger=True)
        result2 = calc_nwalign_glocal(USnw=exe, pdb2=pdbshorter, pdb1=pdblonger, refshorter=False, reflonger=True)

        self.assertEqual(result1, expected)
        self.assertEqual(result2, expected)


    def test_ref_both_raise_valueerror(self):

        with self.assertRaises(ValueError):
            calc_nwalign_glocal(USnw=exe, pdb1=pdbshorter, pdb2=pdblonger, refshorter=True, reflonger=True)
        
        with self.assertRaises(ValueError):
            calc_nwalign_glocal(USnw=exe, pdb2=pdbshorter, pdb1=pdblonger, refshorter=True, reflonger=True)

