import unittest
import sys

sys.path.append('..')
from align_tools import *
from read_fasta import read_prot_from_fasta

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
test_data = f'{test_dir}/data/align_tools'

nw = f'{test_dir}/../../bin/NWalign/align'
pdb2fasta = f'{test_dir}/../../bin/USalign/pdb2fasta'


class TestAlignTools(unittest.TestCase):
    def test_calc_nwalign(self):
        fasta1 = f'{test_data}/O15205.fasta'
        fasta2 = f'{test_data}/P63072.fasta'
        expected1 = 0.71
        expected2 = 0.697
        result1 = calc_nwalign(nw, fasta1, fasta2)
        result2 = calc_nwalign(nw, fasta2, fasta1)
        self.assertEqual(result1, expected1)
        self.assertEqual(result2, expected2)


    def test_fasta_of_pdb(self):
        pdb = f'{test_data}/mychain.pdb'
        with fasta_of_pdb(pdb, pdb2fasta) as fasta:
            header, seq = next(read_prot_from_fasta(fasta))
        self.assertEqual(seq, 'GIVEQCCASVCSLY')


    def test_nw_fasta_to_pdb(self):
        fasta = f'{test_data}/myseq.fasta'
        pdb = f'{test_data}/mychain.pdb'
        result = nw_fasta_to_pdb(fasta, pdb, nw, pdb2fasta)
        expected = 0.786
        self.assertEqual(result, expected)



if __name__ == '__main__':
    unittest.main()
