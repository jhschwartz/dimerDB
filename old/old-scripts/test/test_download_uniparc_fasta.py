import sys
sys.path.append('..')

from download_uniparc_fasta import *
import unittest
import tempfile
from read_fasta import read_prot_from_fasta


class TestDownloadUniparcFasta(unittest.TestCase):
    def test_download(self):
        uniparc = 'UPI000000002B'
        with tempfile.NamedTemporaryFile(mode='w') as f:
            download_fasta(uniparc, f.name)
            f.seek(0)
            header, seq = next(read_prot_from_fasta(f.name))
        self.assertEqual(header, '>UPI000000002B status=active')
        self.assertEqual(seq, 'MSRGALIVFEGLDKSGKTTQCMNIMESIPANTIKYLNFPQRSTVTGKMIDDYLTRKKTYN' + \
                              'DHIVNLLFCANRWEFASFIQEQLEQGITLIVDRYAFSGVAYAAAKGASMTLSKSYESGLP' + \
                              'KPDLVIFLESGSKEINRNVGEEIYEDVTFQQKVLQEYKKMIEEGDIHWQIISSEFEEDVK' + \
                              'KELIKNIVIEAIHTVTGPVGQLWM')




    def test_download_many_singlethread(self):
        uniparcs = ['UPI000000002B', 'UPI000000000B', 'UPI0000000123']
        with tempfile.TemporaryDirectory() as td:
            lib = f'{td}/lib'
            resulting_fastas = download_many_fasta(uniparcs=uniparcs, lib_path=lib, num_workers=1)
            resulting_seqs = []
            for fasta in resulting_fastas:
                _, seq = next(read_prot_from_fasta(fasta))
                resulting_seqs.append(seq)

            self.assertTrue(resulting_seqs[0].startswith('MSRGALIVFEGLDKSGKTTQCMNIMESIPA'))
            self.assertTrue(resulting_seqs[1].startswith('MGTAATIQTPTKLMNKENAEMILEKIVDHI'))
            self.assertTrue(resulting_seqs[2].startswith('MDNSSVCSPNATFCEGDSCLVTESNFNAIL'))
            self.assertTrue(resulting_seqs[2].endswith('GSNILAYWIDGDMDL'))



    def test_download_many_multithread(self):
        uniparcs = ['UPI000000002B', 'UPI000000000B', 'UPI0000000123']
        with tempfile.TemporaryDirectory() as td:
            lib = f'{td}/lib'
            resulting_fastas = download_many_fasta(uniparcs=uniparcs, lib_path=lib, num_workers=3)
            resulting_seqs = []
            for fasta in resulting_fastas:
                _, seq = next(read_prot_from_fasta(fasta))
                resulting_seqs.append(seq)

            self.assertTrue(resulting_seqs[0].startswith('MSRGALIVFEGLDKSGKTTQCMNIMESIPA'))
            self.assertTrue(resulting_seqs[1].startswith('MGTAATIQTPTKLMNKENAEMILEKIVDHI'))
            self.assertTrue(resulting_seqs[2].startswith('MDNSSVCSPNATFCEGDSCLVTESNFNAIL'))
            self.assertTrue(resulting_seqs[2].endswith('GSNILAYWIDGDMDL'))

