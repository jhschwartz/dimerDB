import unittest
import tempfile
import shutil
import glob
import sys
import subprocess

sys.path.append('..')
from parallel_convert_split_cif import *

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()

cif2pdb_exe = f'{test_dir}/../../bin/USalign/cif2pdb'


class TestConvertSplitCif(unittest.TestCase):
    def test_run_cif2pdb(self):
        cif = f'{test_dir}/data/splitchains/1aa0-assembly1.cif'
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(cif, td)
            tmpcif = f'{td}/{os.path.basename(cif)}'
            run_cif2pdb(target=tmpcif, outprefix='tmp', cif2pdb_exe=cif2pdb_exe)
            files = glob.glob(f'{td}/*')
            self.assertEqual(len(files), 4)
            results = [os.path.basename(r) for r in glob.glob(f'{td}/*.pdb')]
            expected = ['tmpA.pdb', 'tmpA-2.pdb', 'tmpA-3.pdb']
            self.assertEqual(sorted(results), sorted(expected))
    

    def test_rename_resulting_pdbs(self):
        cif = f'{test_dir}/data/splitchains/1aa1-assembly1.cif'
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(cif, td)
            tmpcif = f'{td}/{os.path.basename(cif)}'
            run_cif2pdb(target=tmpcif, outprefix='tmp', cif2pdb_exe=cif2pdb_exe)
            rename_resulting_pdbs(td, 'tmp')
            all_pdb = glob.glob(f'{td}/*.pdb')
            model_named_pdb = glob.glob(f'{td}/tmp-*-*.pdb')
            self.assertEqual(all_pdb, model_named_pdb)
            self.assertTrue(os.path.exists(f'{td}/tmp-S-1.pdb'))
            self.assertFalse(os.path.exists(f'{td}/tmpS.pdb'))
            self.assertTrue(os.path.exists(f'{td}/tmp-S-2.pdb'))


    def test_parallel_process_helper(self):
        cif = f'{test_dir}/data/splitchains/1aa7-assembly1.cif'
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(cif, td)
            tmpcif = f'{td}/{os.path.basename(cif)}'
            process_helper(tmpcif, cif2pdb_exe)
            results = set(glob.glob(f'{td}/*'))
            expected = set([tmpcif, f'{td}/1aa7-A-1.pdb', f'{td}/1aa7-B-1.pdb'])
            self.assertEqual(results, expected)


    def test_parallel_convert_split(self):
        cifs = glob.glob(f'{test_dir}/data/splitchains/*.cif')
        self.assertEqual(len(cifs), 24)
        with tempfile.TemporaryDirectory() as td:
            for cif in cifs:
                shutil.copy(cif, td)
            cifs_in_td = glob.glob(f'{td}/*.cif')
            self.assertEqual(len(cifs_in_td), 24)
            parallel_convert_split_rename_cifs(td, cif2pdb_exe, 4)
            all_files = glob.glob(f'{td}/*')
            self.assertEqual(len(all_files), 73)
            pdbs = glob.glob(f'{td}/*-*-*.pdb')
            self.assertEqual(len(pdbs), 49)


    def test_parallel_convert_split_with_recursive(self):
        cifs = glob.glob(f'{test_dir}/data/splitchains/*.cif')
        self.assertEqual(len(cifs), 24)
        with tempfile.TemporaryDirectory() as td:
            for cif in cifs:
                shutil.copy(cif, td)
            os.makedirs(f'{td}/sub')
            rname = os.path.basename(cifs[-1])
            rfile = f'{td}/{rname}'
            os.rename(rfile, f'{td}/sub/{rname}')
            cifs_in_td = glob.glob(f'{td}/*.cif')
            self.assertEqual(len(cifs_in_td), 23)
            parallel_convert_split_rename_cifs(td, cif2pdb_exe, 4)
            for f in glob.glob(f'{td}/sub/*'):
                name = os.path.basename(f)
                os.rename(f, f'{td}/{name}')
            all_files = glob.glob(f'{td}/*')
            self.assertEqual(len(all_files), 74)
            pdbs = glob.glob(f'{td}/*-*-*.pdb')
            self.assertEqual(len(pdbs), 49)


    def test_parallel_convert_split_exe(self):
        cifs = glob.glob(f'{test_dir}/data/splitchains/*.cif')
        self.assertEqual(len(cifs), 24)
        with tempfile.TemporaryDirectory() as td:
            for cif in cifs:
                shutil.copy(cif, td)
            cifs_in_td = glob.glob(f'{td}/*.cif')
            self.assertEqual(len(cifs_in_td), 24)
            subprocess.run(f'python ../parallel_convert_split_cif.py -p {td} -e {cif2pdb_exe} -t 7', shell=True, check=True)
            all_files = glob.glob(f'{td}/*')
            self.assertEqual(len(all_files), 73)
            pdbs = glob.glob(f'{td}/*-*-*.pdb')
            self.assertEqual(len(pdbs), 49)
        


if __name__ == '__main__':
    unittest.main()

