import sys
sys.path.append('..')

from name_fasta import uniparc_fasta
import unittest
import tempfile
from read_fasta import read_prot_from_fasta

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
test_lib = f'{test_dir}/data/name_fasta/lib'

class TestNameFasta(unittest.TestCase):
    def test_name_uniparc_pos(self):
        uniparc_id = 'UPI0000000F08'
        result = uniparc_fasta(uniparc_id=uniparc_id, lib_path=test_lib)
        expected = f'{test_dir}/data/name_fasta/lib/uniparc/08/UPI0000000F08.fasta'
        self.assertEqual(expected, result)

    def test_name_uniparc_neg(self):
        bad_uniparc_id = 'UPI0000000XYZ'
        with self.assertRaises(FileNotFoundError):
            uniparc_fasta(uniparc_id=bad_uniparc_id, lib_path=test_lib)
