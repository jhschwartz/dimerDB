import sys
sys.path.append('..')

import unittest
import tempfile
import tarfile
import os 
import subprocess

from tar_pdb_chains import tar_pdb_chains 

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
indexfile = f'{test_dir}/data/tar_pdb_chains/index.txt'
chains_lib = f'{test_dir}/data/tar_pdb_chains/lib'


class TestTarPDBChains(unittest.TestCase):
    def test_tar_one_div(self):
        with tempfile.TemporaryDirectory() as td:
            contents = tar_pdb_chains(indexfile=indexfile, div='aa', chains_lib=chains_lib, target_dir=td)
            expected = ['5aa5-a1-m1-cA.pdb', \
                        '5aa5-a1-m1-cD.pdb', \
                        '5aa5-a1-m1-cG.pdb', \
                        '5aa5-a1-m1-cK.pdb', \
                        '5aa5-a1-m1-cB.pdb', \
                        '5aa5-a1-m1-cE.pdb', \
                        '5aa5-a1-m1-cH.pdb', \
                        '5aa5-a1-m1-cL.pdb', \
                        '5aa5-a1-m1-cC.pdb', \
                        '5aa5-a1-m1-cF.pdb', \
                        '5aa5-a1-m1-cI.pdb', \
                        '5aa5-a1-m1-cM.pdb' ]
            self.assertEqual(sorted(expected), sorted(contents))

            # try to extract one
            subprocess.run(f'cd {td} && tar -xf {td}/aa.tar.gz 5aa5-a1-m1-cK.pdb', shell=True, check=True)
            self.assertTrue(os.path.exists(f'{td}/5aa5-a1-m1-cK.pdb'))


