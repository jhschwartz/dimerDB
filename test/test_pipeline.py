import unittest
import subprocess
from generate_test_lib import generate_test_lib
import filecmp
import shutil
import random

import os
lib_path = os.environ['lib_path']

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
expected_dir = os.path.join(test_dir, 'expected')

preptest = os.path.join(test_dir, 'prep_test.sh')
runsmk = os.path.join(test_dir, 'run_test_smk.sh')

class TestPipeline(unittest.TestCase):

    
    def setUp(self):
        self.workdir = 'mytmp'+str(random.randint(1,100000))
        self.workdir = os.path.realpath(self.workdir)
        os.mkdir(self.workdir)
        self.lib = os.path.join(self.workdir, 'lib')

        generate_test_lib(  source_lib=lib_path,
                            test_entries_file='test_entries.txt', 
                            dest_lib=self.lib   )

        subprocess.run(f'{preptest} {self.workdir}', shell=True, check=True)

        env = {
            'basepath': self.workdir,
            'conda_env': os.path.join(test_dir, '..', 'env'),
            'lib_path': self.lib, 
            'smk_profile': os.environ['smk_profile']
        }

        env = {**os.environ.copy(), **env}

        subprocess.run(runsmk, shell=True, cwd=self.workdir, env=env)


    def tearDown(self):
        shutil.rmtree(self.workdir)


    def test_end_to_end(self):

        expected_dimers = os.path.join(expected_dir, 'nonred_dimers.txt')
        expected_info = os.path.join(expected_dir, 'nonred_dimers_info.tsv')

        outdir = os.path.join(self.workdir, 'homodimerDB', 'nonredundant')
        result_dimers = os.path.join(outdir, 'dimers.txt')
        result_info = os.path.join(outdir, 'dimers_info.tsv')

        self.assertTrue(filecmp.cmp(expected_dimers, result_dimers, shallow=False))
        self.assertTrue(filecmp.cmp(expected_info, result_info, shallow=False))


