import unittest
import subprocess
import glob
import tempfile
import shutil

class TestParallelSplitChain(unittest.TestCase):
    def test_split(self):
        pdbs = glob.glob('data/splitchains/*.pdb*')
        self.assertEqual(len(pdbs), 26)
        with tempfile.TemporaryDirectory() as td:
            for pdb in pdbs:
                shutil.copy(pdb, td)
            subprocess.run(f'../parallel_split_pdb.pl {td} 4 ../../bin/PDBParser/split_chain', shell=True)
            pdbs = glob.glob(f'{td}/*.pdb*')
            self.assertEqual(len(pdbs), 67)
            subprocess.run(f'../parallel_split_pdb.pl {td} 4 ../../bin/PDBParser/split_chain', shell=True)
            pdbs = glob.glob(f'{td}/*.pdb*')
            self.assertEqual(len(pdbs), 67)


