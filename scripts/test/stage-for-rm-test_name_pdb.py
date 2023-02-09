import sys
sys.path.append('..')

import unittest
import tempfile
from name_pdb import *

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()


class TestNamePDB(unittest.TestCase):
    def test_name_pdb(self):
        lib = f'{test_dir}/data/fake_lib'
        pdb = '6qle'
        chain = 'N'
        model = 1
        result = name_pdb_file(pdb, chain, model, lib)
        self.assertEqual(result, f'{test_dir}/data/fake_lib/rcsb_pdb/ql/6qle-N-1.pdb')

    
    def test_name_pdb_nonexist_allowFalse(self):
        lib = f'{test_dir}/data/fake_lib'
        pdb = '0xyz'
        chain = 'M'
        model = '4'
        with self.assertRaises(FileNotFoundError):
            name_pdb_file(pdb, chain, model, lib)


    def test_name_pdb_nonexist_allowTrue(self):
        lib = f'{test_dir}/data/fake_lib'
        pdb = '0xyz'
        chain = 'M'
        model = '4'
        result = name_pdb_file(pdb, chain, model, lib, True)
        self.assertEqual(result, f'{test_dir}/data/fake_lib/rcsb_pdb/xy/0xyz-M-4.pdb')


    def test_dimer2pdbs(self):
        lib = f'{test_dir}/data/fake_lib'
        dimer_name = '6qle_N_1-7a6w_AAA_42' 
        result1, result2 = dimer2pdbs(dimer_name, lib)
        self.assertEqual(result1, f'{test_dir}/data/fake_lib/rcsb_pdb/ql/6qle-N-1.pdb')
        self.assertEqual(result2, f'{test_dir}/data/fake_lib/rcsb_pdb/a6/7a6w-AAA-42.pdb')


    def test_get_models_of_chain(self):
        lib = f'{test_dir}/data/fake_lib'
        result = get_models_of_chain('7a6w', 'AAA', lib)
        self.assertEqual(result, ['1', '42'])



if __name__ == '__main__':
    unittest.main()

