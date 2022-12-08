import unittest
import sys
sys.path.append('..')
from read_fasta import read_prot_from_fasta

fastafile1 = "data/example_db.fasta"
fastafile2 = "data/example_multiline_seqs.fasta"

class TestReadFasta(unittest.TestCase):
    def test_read_each_prot(self):
        prots = read_prot_from_fasta(fastafile1)

        expected_header = '>AFDB:AF-H5USR7-F1 NAD(P)-bd_dom domain-containing protein UA=H5USR7 UI=H5USR7_9MICO OS=Mobilicoccus pelagius NBRC 104925 OX=1089455 GN=MOPEL_080_00540'
        expected_seq = 'MNDSSRTALVTGVTGYIGGSLVRPLLDAGWDVRALTRSAAKLDRRPWRDEIDVVEGDATKSDDLDRALQGVDVAYYLVHSMDAQGDFERRDAELAEDFAAAAERAGVSRLVYLGGLHPDVPDDELSPHLRSRVEVGEIFLESDVPAACLQAGIVLGAGSASFQMLASLTERLPVMIAPKWLRNRIQPIGIDDAVYYLVHAADLPSDVNRTFDIGGADVLTYTEMIQRYAHIAGLRERVIAPVPVFTPDLASHWVGLVTPIDAGLARPLVGSLVHEVICKEHDIAEYVPDPPDGLIGFDEAVEEALADPTS'
        header, seq = next(prots)
        self.assertEqual(expected_header, header)
        self.assertEqual(expected_seq, seq)


        expected_header = '>AFDB:AF-D7WF56-F1 DHHA1 domain protein UA=D7WF56 UI=D7WF56_9CORY OS=Corynebacterium genitalium ATCC 33030 OX=585529 GN=HMPREF0291_11392'
        expected_seq = 'MAAMPGLFAPGEGQFAAVAQRLSEAPVVHVVTHIRPDADAIGSASALALGLQSIGVEATVRIGQHEEIPENLATIPGVTDVILGEPLPADGLVVTTDCASIDRTGSYRSVIQQERDRVVVIDHHESNPGFGSMNLVLPSESTTVIIRELFGHLGVEITPEMAYCLYAGLVTDTGSFRWGTPRMHLLAAELMGLGVEPRSAAMALMDAVSAEDLRLMGAVMAGMETRVAGHHTMSVLAIDDEHLHGMNQTAVESIIEYSRALTGSDIGVVFKELHPGYWSVSLRSSVVNVAEIAAVHGGGGHIPAAGYSARGPVADAIDALHATVAALDATAL'
        header, seq = next(prots)
        self.assertEqual(expected_header, header)
        self.assertEqual(expected_seq, seq)



    def test_read_each_multiline_prot(self):
        prots = read_prot_from_fasta(fastafile2)

        expected_header = '>this is a fake protein'
        expected_seq = 'ABCDEFGHIJKLZZZZZZZZZZZZ'
        header, seq = next(prots)
        self.assertEqual(expected_header, header)
        self.assertEqual(expected_seq, seq)

        expected_header = '>yet another fake protein'
        expected_seq = 'KJHDSKLJFHKLJBEKJLFNSKJNKEBFKENDCSKJLHLKWDJ'
        header, seq = next(prots)
        self.assertEqual(expected_header, header)
        self.assertEqual(expected_seq, seq)

        expected_header = '>one more fakey'
        expected_seq = 'ABCDEELKFNKLENFKLSKLWKLWLKLKMWSLKMWSMKLMWDSLKMWLKMWLKM'
        header, seq = next(prots)
        self.assertEqual(expected_header, header)
        self.assertEqual(expected_seq, seq)

if __name__=='__main__':
    t = TestReadFasta()
    t.test_read_each_multiline_prot()
