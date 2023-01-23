import sys
sys.path.append('..')

from name_fasta import *
import unittest
import tempfile
from read_fasta import read_prot_from_fasta

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()

class TestNameFasta(unittest.TestCase):
    def test_name_homodimer_pos(self):
        filtering_dir = f'{test_dir}/data/fake_homodimer_filtering'
        dimer_name = 'UPI0000000F08'
        fasta = get_homodimer_fasta(dimer_name, filtering_dir)
        self.assertEqual(fasta, f'{test_dir}/data/fake_homodimer_filtering/08/UPI0000000F08/seq.fasta')


    def test_name_homodimer_neg(self):
        filtering_dir = f'{test_dir}/data/fake_homodimer_filtering'
        dimer_name = 'UPI0000000XYZ'
        with self.assertRaises(FileNotFoundError):
            get_homodimer_fasta(dimer_name, filtering_dir)      

    def test_name_heterodimer_pos(self):
        filtering_dir = f'{test_dir}/data/fake_heterodimer_filtering'
        dimer_name = 'UPI0000000005-UPI00000F8202'
        fasta1, fasta2 = get_heterodimer_fasta(dimer_name, filtering_dir)
        expected1 = f'{test_dir}/data/fake_heterodimer_filtering/5-2/UPI0000000005-UPI00000F8202/seq1.fasta'
        expected2 = f'{test_dir}/data/fake_heterodimer_filtering/5-2/UPI0000000005-UPI00000F8202/seq2.fasta'
        self.assertEqual(fasta1, expected1)
        self.assertEqual(fasta2, expected2)
    
    
    def test_name_heterodimer_neg(self):
        filtering_dir = f'{test_dir}/data/fake_heterodimer_filtering'
        dimer_name = 'UPI0000000XYZ-UPI0000000ABC'
        with self.assertRaises(FileNotFoundError):
            get_homodimer_fasta(dimer_name, filtering_dir)      


