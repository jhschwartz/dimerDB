import sys
sys.path.append('..')

from download_uniparc_fasta import download_fasta
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
