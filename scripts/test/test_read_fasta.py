import os
import sys
sys.path.append('..')

import unittest
from read_fasta import read_prot_from_fasta


import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = f'{test_dir}/data/read_fasta'
lib_path = os.path.join(data_dir, 'lib')


class TestReadFasta(unittest.TestCase):
    
    def test_read_fasta_one_seq(self):

        infile = os.path.join(data_dir, 'oneseq.fasta')

        header, seq = next(read_prot_from_fasta(infile))

        self.assertEqual(header, '>THIS is A HeAdEr')
        self.assertEqual(seq, 'AABBCCDDEEFFG123HJKLOPMMMMMMMMMMMMMMMMMMMMMMM')
    
    
    
    def test_read_fasta_multi_seq(self):

        infile = os.path.join(data_dir, 'multiseq.fasta')

        headers = [ '>i am a protein', '>i am also a protein', '>yet another protein' ]
        seqs = [ 'A',
                 'ABCHJKKSAKLJLAKDJFLKAJDFLKDSJFLDKJF',
                 ''.join(['A' for _ in range(60*4)])    ]

        for i, (header, seq) in enumerate(read_prot_from_fasta(infile)):
            self.assertEqual(header, headers[i])
            self.assertEqual(seq, seqs[i])


