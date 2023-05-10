import sys
sys.path.append('..')

import unittest
import tempfile
import tarfile
import os 
from tar_pdb_chains import tar_pdb_chains 

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
indexfile = f'{test_dir}/data/tar_pdb_chains/index.txt'
chains_lib = f'{test_dir}/data/tar_pdb_chains/chains'


class TestTarPDBChains(unittest.TestCase):
    def test_tar_simple(self):
        with tempfile.TemporaryDirectory() as td:
            results = tar_pdb_chains(indexfile=indexfile, chains_lib=chains_lib, target_dir=td)
            self.assertEqual(10, len(results))

            self.assertTrue('aa.tar.gz' in [os.path.basename(f) for f in results])

