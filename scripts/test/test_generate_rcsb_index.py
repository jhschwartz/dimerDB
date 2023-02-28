import sys
sys.path.append('..')

import unittest
import tempfile
from generate_rcsb_index import *
import pickle
import shutil
import os


import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
lib_path = f'{test_dir}/data/generate_rcsb_index/lib'
rcsb_path = f'{test_dir}/data/generate_rcsb_index/lib/rcsb'


class TestGenerateRcsbIndex(unittest.TestCase):
    def test_gen_index(self):
        with tempfile.TemporaryDirectory() as td:
            index_file = f'{td}/index.txt'
            generate_rcsb_index(rcsb_path, index_file)
            with open(index_file) as f:
                pdbs = [line.strip() for line in f]
            expected = ['6qle-a1-m1-cH.pdb', 
                        '6qle-a4-m2-cN.pdb',
                        '7a6w-a1-m1-cAAA.pdb',
                        '7a6w-a1-m42-cAAA.pdb',
                        '7a6w-a2-m42-cAAA.pdb']
            self.assertEqual(expected, pdbs)


    def test_convert_index_to_pkl(self):
        with tempfile.TemporaryDirectory() as td:
            index_pkl = f'{td}/index.pkl'
            index_txt = f'{lib_path}/rcsb_index.txt' 
            rcsb_index_to_pkl(index_txt, index_pkl)
            with open(index_pkl, 'rb') as f:
                result = pickle.load(f)
            expected = {
                '6qle': {
                    'H': ['6qle-a1-m1-cH.pdb'],
                    'N': ['6qle-a4-m2-cN.pdb']
                },
                '7a6w': {
                    'AAA': ['7a6w-a1-m1-cAAA.pdb', '7a6w-a1-m42-cAAA.pdb', '7a6w-a2-m42-cAAA.pdb']
                }
            }
            self.assertEqual(result, expected)


    def test_detect_empty_pdbs(self):
        with tempfile.TemporaryDirectory() as td:
            index_file = f'{td}/index.txt'
            shutil.copytree(rcsb_path, f'{td}/rcsb')
            os.remove(f'{td}/rcsb/ql/6qle-a1-m1-cH.pdb')
            pathlib.Path(f'{td}/rcsb/ql/6qle-a1-m1-cH.pdb').touch()
            missing_files = generate_rcsb_index(f'{td}/rcsb', index_file)
            expected_missing = [f'{td}/rcsb/ql/6qle-a1-m1-cH.pdb']

            with open(index_file) as f:
                pdbs = [line.strip() for line in f]
            expected_pdbs = ['6qle-a1-m1-cH.pdb',
                        '6qle-a4-m2-cN.pdb',
                        '7a6w-a1-m1-cAAA.pdb',
                        '7a6w-a1-m42-cAAA.pdb',
                        '7a6w-a2-m42-cAAA.pdb']

            self.assertEqual(expected_missing, missing_files)
            self.assertEqual(expected_pdbs, pdbs)

