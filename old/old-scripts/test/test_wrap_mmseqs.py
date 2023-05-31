import sys
sys.path.append('..')

import unittest
import tempfile
import os
import filecmp
from wrap_mmseqs import *

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
datadir = f'{test_dir}/data/wrap_mmseqs'


class TestWrapMMseqs(unittest.TestCase):
    def test_cluster_30(self):
        infasta = os.path.join(datadir, 'in.fasta')
        expected_tsv = os.path.join(datadir, 'expected30.tsv')
        expected_fasta = os.path.join(datadir, 'expected30.fasta')
        with tempfile.TemporaryDirectory() as td:
            pref = os.path.join(td, 'out')
            result_tsv, result_fasta = \
                        run_mmseqs_cluster(  infasta=infasta,
                                             outprefix=pref,
                                             seq_id=30,
                                             mmseqs_exe='mmseqs',
                                             cores=8              )
            self.assertTrue(filecmp.cmp(expected_tsv, result_tsv))
            self.assertTrue(filecmp.cmp(expected_fasta, result_fasta))


    def test_cluster_80(self):
        infasta = os.path.join(datadir, 'in.fasta')
        expected_tsv = os.path.join(datadir, 'expected80.tsv')
        expected_fasta = os.path.join(datadir, 'expected80.fasta')
        with tempfile.TemporaryDirectory() as td:
            pref = os.path.join(td, 'out')
            result_tsv, result_fasta = \
                        run_mmseqs_cluster(  infasta=infasta,
                                             outprefix=pref,
                                             seq_id=80,
                                             mmseqs_exe='mmseqs',
                                             cores=8              )
            self.assertTrue(filecmp.cmp(expected_tsv, result_tsv))
            self.assertTrue(filecmp.cmp(expected_fasta, result_fasta))



    def test_read_reps_30(self):
        expected = ['A0A086N7H4', 'A0A418GFQ3', 'A0A484Y3A2', 'K7B3U4', 'O95905', 'P04439', 'P08246', 'Q5RCA9']
        tsv = os.path.join(datadir, 'expected30.tsv')
        result = derive_mmseqs_reps_from_tsv(tsv)
        print('expected', expected)
        print('resulted', result)
        self.assertEqual(expected, result)


    def test_read_reps_80(self):
        expected = ['A0A086N7H4', 'A0A0Q5HG13', 'A0A2K5L455', 'A0A418GFQ3', 'A0A484Y3A2', 'A0A4Q1R4Z3', 'A0A5N3WSI5', 'P04439', 'P08246', 'Q5RCA9']
        tsv = os.path.join(datadir, 'expected80.tsv')
        result = derive_mmseqs_reps_from_tsv(tsv)
        self.assertEqual(expected, result)

