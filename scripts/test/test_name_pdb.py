import sys
sys.path.append('..')

import unittest
import tempfile
from name_pdb import *

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()


class TestNamePDB(unittest.TestCase):
    def test_name_pdb_reg(self):
        lib = f'{test_dir}/data/fake_lib'
        pdb = '6qle'
        chain = 'N'
        result = name_pdb_file(pdb, chain, lib)
        self.assertEqual(result, f'{test_dir}/data/fake_lib/rcsb_pdb/ql/6qleN.pdb')

    
    def test_name_pdb_oversized(self):
        lib = f'{test_dir}/data/fake_lib'
        pdb = '7a6w'
        chain = 'AAA'
        result = name_pdb_file(pdb, chain, lib)
        self.assertEqual(result, f'{test_dir}/data/fake_lib/rcsb_oversized/a6/7a6w/split_renamed/7a6wAAA.pdb')


    def test_name_pdb_nonexist(self):
        lib = f'{test_dir}/data/fake_lib'
        pdb = '0xyz'
        chain = 'M'
        with self.assertRaises(FileNotFoundError):
            name_pdb_file(pdb, chain, lib)


    def test_dimer2pdbs(self):
        lib = f'{test_dir}/data/fake_lib'
        dimer_name = '6qle_N-7a6w_AAA' 
        result1, result2 = dimer2pdbs(dimer_name, lib)
        self.assertEqual(result1, f'{test_dir}/data/fake_lib/rcsb_pdb/ql/6qleN.pdb')
        self.assertEqual(result2, f'{test_dir}/data/fake_lib/rcsb_oversized/a6/7a6w/split_renamed/7a6wAAA.pdb')


