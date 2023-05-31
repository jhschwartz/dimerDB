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
test_lib = f'{test_dir}/data/parallel_convert_split_cif/lib'

cif2pdb_exe = f'{test_dir}/../../bin/USalign/cif2pdb'
python = '/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python'


class TestConvertSplitCif(unittest.TestCase):
    def test_run_cif2pdb(self):
        cif = f'{test_lib}/rcsb/aa/1aa0-assembly1.cif'
        with tempfile.TemporaryDirectory() as td1, tempfile.TemporaryDirectory() as td2:
            shutil.copy(cif, td1)
            tmpcif = f'{td1}/{os.path.basename(cif)}'
            run_cif2pdb(cif=tmpcif, outdir=td2, prefix='tmp', cif2pdb_exe=cif2pdb_exe)
            pdbs = glob.glob(f'{td2}/*')
            cifs = glob.glob(f'{td1}/*')
            self.assertEqual(len(pdbs), 3)
            self.assertEqual(len(cifs), 1)
            results = [os.path.basename(r) for r in pdbs]
            expected = ['tmpA.pdb', 'tmpA-2.pdb', 'tmpA-3.pdb']
            self.assertEqual(sorted(results), sorted(expected))
    

    def test_rename_resulting_pdbs_1(self):
        cif = f'{test_lib}/rcsb/aa/1aa1-assembly1.cif'
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(cif, td)
            tmpcif = f'{td}/{os.path.basename(cif)}'
            run_cif2pdb(cif=tmpcif, outdir=td, prefix='tmp', cif2pdb_exe=cif2pdb_exe)
            rename_resulting_pdbs(dir_=td, assembly_ID='1', prefix='tmp')
            all_pdb = glob.glob(f'{td}/*.pdb')
            model_named_pdb = glob.glob(f'{td}/tmp-*-*.pdb')
            self.assertEqual(all_pdb, model_named_pdb)
            self.assertTrue(os.path.exists(f'{td}/tmp-a1-m1-cS.pdb'))
            self.assertFalse(os.path.exists(f'{td}/tmpS.pdb'))
            self.assertTrue(os.path.exists(f'{td}/tmp-a1-m1-cS.pdb'))


    def test_rename_resulting_pdbs_2(self):
        cif = f'{test_lib}/rcsb/aa/1aa1-assembly1.cif'
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(cif, td)
            tmpcif = f'{td}/{os.path.basename(cif)}'
            run_cif2pdb(cif=tmpcif, outdir=td, prefix='tmp', cif2pdb_exe=cif2pdb_exe)
            rename_resulting_pdbs(dir_=td, assembly_ID='2', prefix='tmp')
            all_pdb = glob.glob(f'{td}/*.pdb')
            model_named_pdb = glob.glob(f'{td}/tmp-a*-m*-c*.pdb')
            self.assertEqual(all_pdb, model_named_pdb)
            self.assertTrue(os.path.exists(f'{td}/tmp-a2-m1-cS.pdb'))
            self.assertFalse(os.path.exists(f'{td}/tmpS.pdb'))
            self.assertTrue(os.path.exists(f'{td}/tmp-a2-m1-cS.pdb'))


    def test_parallel_process_helper(self):
        cif = f'{test_lib}/rcsb/aa/1aa7-assembly96.cif'
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(cif, td)
            tmpcif = f'{td}/{os.path.basename(cif)}'
            process_helper(tmpcif, cif2pdb_exe)
            results = set(glob.glob(f'{td}/*'))
            expected = set([tmpcif, f'{td}/1aa7-a96-m1-cA.pdb', f'{td}/1aa7-a96-m1-cB.pdb'])
            self.assertEqual(results, expected)


    def test_parallel_convert_split(self):
        cifs = glob.glob(f'{test_lib}/rcsb/**/*.cif', recursive=True)
        self.assertEqual(len(cifs), 24)
        with tempfile.TemporaryDirectory() as td:
            for cif in cifs:
                shutil.copy(cif, td)
            cifs_in_td = glob.glob(f'{td}/*.cif')
            self.assertEqual(len(cifs_in_td), 24)
            parallel_convert_split_rename_cifs(td, cif2pdb_exe, 4)
            all_files = glob.glob(f'{td}/*')
            self.assertEqual(len(all_files), 73)
            pdbs = glob.glob(f'{td}/*-*-*-*.pdb')
            self.assertEqual(len(pdbs), 49)


    def test_parallel_convert_split_with_recursive(self):
        cifs = glob.glob(f'{test_lib}/rcsb/**/*.cif', recursive=True)
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
            pdbs = glob.glob(f'{td}/*-*-*-*.pdb')
            self.assertEqual(len(pdbs), 49)


    def test_parallel_convert_split_with_recursive_whole_lib(self):
        parent = f'{test_lib}/rcsb'
        with tempfile.TemporaryDirectory() as td:
            shutil.copytree(parent, f'{td}/rcsb')
            before_cifs = glob.glob(f'{td}/**/*.cif', recursive=True)
            before_pdbs = glob.glob(f'{td}/**/*.pdbs', recursive=True)
            parallel_convert_split_rename_cifs(td, cif2pdb_exe, 4)
            after_cifs = glob.glob(f'{td}/**/*.cif', recursive=True)
            after_pdbs = glob.glob(f'{td}/**/*.pdb', recursive=True)
            after_pdbs_formatted = glob.glob(f'{td}/**/*-a*-m*-c*.pdb', recursive=True)
            self.assertEqual(len(before_cifs), 24)
            self.assertEqual(len(before_pdbs), 0)
            self.assertEqual(len(after_cifs), 24)
            self.assertEqual(len(after_pdbs), 49)
            self.assertEqual(len(after_pdbs_formatted), 49)


    def test_parallel_convert_split_exe(self):
        parent = f'{test_lib}/rcsb'
        with tempfile.TemporaryDirectory() as td:
            shutil.copytree(parent, f'{td}/rcsb')
            before_cifs = glob.glob(f'{td}/**/*.cif', recursive=True)
            before_pdbs = glob.glob(f'{td}/**/*.pdbs', recursive=True)
            subprocess.run(f'{python} ../parallel_convert_split_cif.py -p {td} -e {cif2pdb_exe} -t 7', shell=True, check=True)
            after_cifs = glob.glob(f'{td}/**/*.cif', recursive=True)
            after_pdbs = glob.glob(f'{td}/**/*.pdb', recursive=True)
            after_pdbs_formatted = glob.glob(f'{td}/**/*-a*-m*-c*.pdb', recursive=True)
            self.assertEqual(len(before_cifs), 24)
            self.assertEqual(len(before_pdbs), 0)
            self.assertEqual(len(after_cifs), 24)
            self.assertEqual(len(after_pdbs), 49)
            self.assertEqual(len(after_pdbs_formatted), 49)


    def test_fill_empty_pdb(self):
        with tempfile.TemporaryDirectory() as td:
            shutil.copytree(f'{test_lib}/rcsb/aa/', f'{td}/aa')
            parallel_convert_split_rename_cifs(td, cif2pdb_exe, 1)
            pdb_path = f'{td}/aa/1aa0-a1-m3-cA.pdb'
            with open(pdb_path, 'r') as f:
                old_pdb = f.read()
            os.remove(pdb_path) 
            fill_empty_pdb(filepath=pdb_path, cif2pdb_exe=cif2pdb_exe)
            with open(pdb_path, 'r') as f:
                new_pdb = f.read()
            self.assertEqual(old_pdb, new_pdb)


if __name__ == '__main__':
    unittest.main()

