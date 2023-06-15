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

    
    def test_name_pdb_dashedChain(self):
        pdb = '6qle'
        chain = 'N-42'
        model = '1'
        assembly = '1'
        result = name_pdb_file(pdb, assembly, model, chain, lib)
        self.assertEqual(result, f'{lib}/rcsb/ql/6qle-a1-m1-cN-42.pdb')

    
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
        result = name_chain_from_filename('abcd-a3-m22-cG.pdb')
        expected = 'abcd-a3-m22-cG'
        self.assertEqual(result, expected)


    def test_name_chain_dashedChain(self):
        result = name_chain_from_filename('abcd-a3-m22-cG-42.pdb')
        expected = 'abcd-a3-m22-cG-42'
        self.assertEqual(result, expected)


    def test_dimer2pdbs(self):
        dimer_name = '6qle-a1-m1-cN_7a6w-a1-m42-cAAA' 
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


    def test_read_chain_name_filename_dashedChain(self):
        name = '4zbw-a1-m2-cZ-42.pdb'
        result = read_chain_names(name)
        expected = ('4zbw', '1', '2', 'Z-42')
        self.assertEqual(result, expected)


    def test_read_chain_name_nameonly_1(self):
        name = '7a6w-a155-m42-cA2fA'
        result = read_chain_names(name)
        expected = ('7a6w', '155', '42', 'A2fA')
        self.assertEqual(result, expected)


    def test_read_chain_name_nameonly_2(self):
        name = '4zbw-a1-m2-cZ'
        result = read_chain_names(name)
        expected = ('4zbw', '1', '2', 'Z')
        self.assertEqual(result, expected)


    def test_read_chain_name_nameonly_dashedChain(self):
        name = '4zbw-a1-m2-cZ-42'
        result = read_chain_names(name)
        expected = ('4zbw', '1', '2', 'Z-42')
        self.assertEqual(result, expected)


    def test_read_chain_name_mismatch_1(self):
        name = '7a6w_a155_m42_cA2fA.pdb'
        with self.assertRaises(ValueError):
            read_chain_names(name)


    def test_read_chain_name_mismatch_2(self):
        name = '7a6w_a155_m42_cA2fA'
        with self.assertRaises(ValueError):
            read_chain_names(name)


    def test_read_chain_name_mismatch_3(self):
        name = 'some/garbage__'
        with self.assertRaises(ValueError):
            read_chain_names(name)


    def test_get_div_no_ext(self):
        name = '7a6w-a155-m42-cA2fA'
        div = get_div(name)
        self.assertEqual(div, 'a6')


    def test_get_div_with_ext(self):
        name = '7a6w-a155-m42-cA2fA.pdb'
        div = get_div(name)
        self.assertEqual(div, 'a6')

    
    def test_get_div_abs_path_name(self):
        name = '/nfs/turbo/umms-helloworld/notreal/7a6w-a155-m42-cA2fA.pdb'
        div = get_div(name)
        self.assertEqual(div, 'a6')


    def test_name_dimer_sort_needed(self):
        c1 = '7a6w-a155-m42-cA2fA'
        c2 = '4zbw-a1-m2-cZ-42'
        result = name_dimer(c1, c2)
        self.assertEqual(result, '4zbw-a1-m2-cZ-42_7a6w-a155-m42-cA2fA')
    
                
    def test_name_dimer__nosort_needed(self):
        c2 = '7a6w-a155-m42-cA2fA'
        c1 = '4zbw-a1-m2-cZ-42'
        result = name_dimer(c1, c2)
        self.assertEqual(result, '4zbw-a1-m2-cZ-42_7a6w-a155-m42-cA2fA')


    def test_sort_key_dimer_in(self):
        result = sort_key('4zbw-a1-m2-cZ-42_7a6w-a155-m42-cA2fA')
        expected = ['4zbw-a1-m2-cZ-42', '7a6w-a155-m42-cA2fA']
        self.assertEqual(result, expected)


    def test_sort_key_chains_in(self):
        result = sort_key('4zbw-a1-m2-cZ-42', '7a6w-a155-m42-cA2fA')
        expected = ['4zbw-a1-m2-cZ-42', '7a6w-a155-m42-cA2fA']
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()

