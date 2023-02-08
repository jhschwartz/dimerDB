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

