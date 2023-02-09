import sys
sys.path.append('..')

import unittest
import tempfile
from name_pdb import *

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
lib = f'{test_dir}/data/name_pdb/lib'


class TestNamePDB(unittest.TestCase):
    def test_name_pdb_1(self):
        pdb = '6qle'
        chain = 'N'
        model = '1'
        assembly = '1'
        result = name_pdb_file(pdb, assembly, model, chain, lib)
        self.assertEqual(result, f'{lib}/rcsb/ql/6qle-a1-m1-cN.pdb')

    
    def test_name_pdb_2(self):
        pdb = '6qle'
        chain = 'N'
        model = '2'
        assembly = '4'
        result = name_pdb_file(pdb, assembly, model, chain, lib)
        self.assertEqual(result, f'{lib}/rcsb/ql/6qle-a4-m2-cN.pdb')

    
    def test_name_pdb_nonexist_allowFalse(self):
        pdb = '0xyz'
        chain = 'M'
        model = '4'
        assembly = '19'
        with self.assertRaises(FileNotFoundError):
            name_pdb_file(pdb, assembly, model, chain, lib)


    def test_name_pdb_nonexist_allowTrue(self):
        pdb = '0xyz'
        chain = 'M'
        model = '4'
        assembly = '19'
        result = name_pdb_file(pdb, assembly, model, chain, lib, True)
        self.assertEqual(result, f'{lib}/rcsb/xy/0xyz-a19-m4-cM.pdb')


    def test_name_pdb_noLib(self):
        pdb = '0xyz'
        chain = 'M'
        model = '4'
        assembly = '19'
        result = name_pdb_file(pdb, assembly, model, chain)
        self.assertEqual(result, '0xyz-a19-m4-cM.pdb')


    def test_name_chain(self):
        result = name_chain('abcd', '3', '22', 'G')
        expected = 'abcd_a3_m22_cG'
        self.assertEqual(result, expected)


    def test_dimer2pdbs(self):
        dimer_name = '6qle_a1_m1_cN-7a6w_a1_m42_cAAA' 
        result1, result2 = dimer2pdbs(dimer_name, lib)
        self.assertEqual(result1, f'{lib}/rcsb/ql/6qle-a1-m1-cN.pdb')
        self.assertEqual(result2, f'{lib}/rcsb/a6/7a6w-a1-m42-cAAA.pdb')


    def test_read_chain_name_filename_fullpath_1(self):
        name = '/path/that/is/full/7a6w-a155-m42-cA2fA.pdb'
        result = read_chain_names(name)
        expected = ('7a6w', '155', '42', 'A2fA')
        self.assertEqual(result, expected)


    def test_read_chain_name_filename_fullpath_2(self):
        name = '/path/that/is/full/4zbw-a1-m2-cZ.pdb'
        result = read_chain_names(name)
        expected = ('4zbw', '1', '2', 'Z')
        self.assertEqual(result, expected)


    def test_read_chain_name_filename_localpath_1(self):
        name = '7a6w-a155-m42-cA2fA.pdb'
        result = read_chain_names(name)
        expected = ('7a6w', '155', '42', 'A2fA')
        self.assertEqual(result, expected)


    def test_read_chain_name_filename_localpath_2(self):
        name = '4zbw-a1-m2-cZ.pdb'
        result = read_chain_names(name)
        expected = ('4zbw', '1', '2', 'Z')
        self.assertEqual(result, expected)


    def test_read_chain_name_nameonly_1(self):
        name = '7a6w_a155_m42_cA2fA'
        result = read_chain_names(name)
        expected = ('7a6w', '155', '42', 'A2fA')
        self.assertEqual(result, expected)


    def test_read_chain_name_nameonly_2(self):
        name = '4zbw_a1_m2_cZ'
        result = read_chain_names(name)
        expected = ('4zbw', '1', '2', 'Z')
        self.assertEqual(result, expected)


    def test_read_chain_name_mismatch_1(self):
        name = '7a6w_a155_m42_cA2fA.pdb'
        with self.assertRaises(ValueError):
            read_chain_names(name)


    def test_read_chain_name_mismatch_2(self):
        name = '7a6w-a155-m42-cA2fA'
        with self.assertRaises(ValueError):
            read_chain_names(name)


    def test_read_chain_name_mismatch_3(self):
        name = 'some/garbage__'
        with self.assertRaises(ValueError):
            read_chain_names(name)


    def test_instances_of_chain(self):
        result = get_instances_of_chain('7a6w', 'AAA', lib)

        expected_instances = [
                ('7a6w', '1', '1', 'AAA'),
                ('7a6w', '1', '42', 'AAA'),
                ('7a6w', '2', '42', 'AAA')
        ]

        for item in expected_instances:
            self.assertTrue(item in result)




if __name__ == '__main__':
    unittest.main()

