import unittest
import sys
import tempfile
from distutils.dir_util import copy_tree
import glob

sys.path.append('..')
import rename_oversized_chains as oversized_pdb

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()


class OversizedPDB_Test(unittest.TestCase):
    def test_read_bundle_mappingline_1(self):
        result = oversized_pdb._bundle_name_from_mappingline('7w5z-pdb-bundle1.pdb:')
        self.assertEqual(result, 'bundle1')
  

    def test_read_bundle_mappingline_2(self):
        result = oversized_pdb._bundle_name_from_mappingline('7w5z-pdb-bundle12.pdb:')
        self.assertEqual(result, 'bundle12')


    def test_read_bundle_filename_1(self):
        result = oversized_pdb._bundle_name_from_filename('7w5z-pdb-bundle1A.pdb')
        self.assertEqual(result, 'bundle1')
  

    def test_read_bundle_filename_2(self):
        result = oversized_pdb._bundle_name_from_filename('7w5z-pdb-bundle14A.pdb')
        self.assertEqual(result, 'bundle14')


    def test_read_mappings(self):
        mappingfile = f'{test_dir}/data/fake_oversized_mapping.txt' 
        mappings = oversized_pdb._read_mappings(mappingfile)
        self.assertEqual(mappings['bundle1']['A'], 'U1')
        self.assertEqual(mappings['bundle1']['B'], 'U2')
        self.assertEqual(mappings['bundle1']['I'], 'C3')
        self.assertEqual(mappings['bundle1']['Q'], '7L')
        self.assertEqual(mappings['bundle2']['A'], 'u5')
        self.assertEqual(mappings['bundle2']['B'], 'u6')
        self.assertEqual(mappings['bundle2']['I'], '6l')
        self.assertEqual(mappings['bundle2']['P'], 'm3')
        self.assertEqual(mappings['bundle2']['T'], 't4')


    def test_chain_from_filename_1(self):
        result = oversized_pdb._chain_name_from_filename('7w5z-pdb-bundle1A.pdb')
        self.assertEqual(result, 'A')


    def test_chain_from_filename_2(self):
        result = oversized_pdb._chain_name_from_filename('7w5z-pdb-bundle11A.pdb')
        self.assertEqual(result, 'A')


    def test_derive_new_name(self):
        mappingfile = f'{test_dir}/data/fake_oversized_mapping.txt' 
        mappings = oversized_pdb._read_mappings(mappingfile)
        result = oversized_pdb._derive_new_name('7w5z-pdb-bundle1A.pdb', mappings)
        self.assertEqual(result, '7w5zU1.pdb')


    def test_file_is_single_chain_pos(self):
        mappingfile = f'{test_dir}/data/fake_oversized_mapping.txt'
        mappings = oversized_pdb._read_mappings(mappingfile)
        result = oversized_pdb._file_is_single_chain('7w5z-pdb-bundle2Z.pdb', mappings)
        self.assertTrue(result)


    def test_file_is_single_chain_neg(self):
        mappingfile = f'{test_dir}/data/fake_oversized_mapping.txt'
        mappings = oversized_pdb._read_mappings(mappingfile)
        result = oversized_pdb._file_is_single_chain('7w5z-pdb-bundle2.pdb', mappings)
        self.assertFalse(result)


    def test_rename(self):
        data_folder = f'{test_dir}/data/stored-7a01-pdb-bundle'
        with tempfile.TemporaryDirectory() as td:
            copy_tree(data_folder, td)
            oversized_pdb._rename_chainfiles_in_dir(td, '7a01')
            rechained_files = [fn for fn in glob.glob(f'{td}/split_renamed/*.pdb') if not 'bundle' in fn]
            self.assertEqual(len(rechained_files), 85)




if __name__ == '__main__':
    unittest.main()
