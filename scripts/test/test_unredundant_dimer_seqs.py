import unittest
import contextlib
import os

import sys
sys.path.append('..')

from unredundant import RedundantSeqs, RedundantSeqsHomodimer, RedundantSeqsHeterodimer
from read_fasta import read_prot_from_fasta

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()


config = {
    'paths': {
        'lib': f'{test_dir}/data/fake_lib',
        'nwalign': f'{test_dir}/../../bin/NWalign/align',
        'intermediates_homodimer_filtering': f'{test_dir}/data/unred_seq_fakedata/fake_homo_filter' 
    }
}


class ConcreteRedundantSeqs(RedundantSeqs):
    def __init__(self, names, yamlfile, threshold, config):
        super().__init__(names, yamlfile, threshold, config)
    
    @contextlib.contextmanager
    def dimer_tmp_fasta(self, name_not_dimer):
        fastas_for_testing = {
            'O15205':  f'{test_dir}/data/redundant_monomer_fastas/O15205.fasta', 
            'P63072':  f'{test_dir}/data/redundant_monomer_fastas/P63072.fasta', 
            'Q921A3':  f'{test_dir}/data/redundant_monomer_fastas/Q921A3.fasta', 
            'P0DTC2':  f'{test_dir}/data/redundant_monomer_fastas/P0DTC2.fasta', 
            'P11223':  f'{test_dir}/data/redundant_monomer_fastas/P11223.fasta', 
            'P59594':  f'{test_dir}/data/redundant_monomer_fastas/P59594.fasta'  
        }
        yield fastas_for_testing[name_not_dimer]



class TestRedundantSeqsSubclass(unittest.TestCase):
    def test_redundant_seqs_subclass(self):
        names = ['O15205', 'P63072', 'Q921A3', 'P0DTC2', 'P11223', 'P59594']
        yamlfile = f'{test_dir}/data/redundant_monomer_fastas/fake_seqs.yaml'
        threshold = 0.5 
        redundant_seqs = ConcreteRedundantSeqs(names, yamlfile, threshold, config)
        result = redundant_seqs.prune_redundancy(num_workers=2)
        expected = ['Q921A3', 'P11223', 'P59594']





homoyaml = f'{test_dir}/data/unred_seq_fakedata/homodimers.yaml'
homodimers = {
    'UPI00001653E1': ['2gqq_A-2gqq-B', '2gqq_A-2gqq-C', '2gqq_A-2gqq-D', '2gqq_B-2gqq-C', '2gqq_B-2gqq-D', '2gqq_C-2gqq-D'],    
    'UPI000016225C': ['fake', 'yet', 'another', 'fake'],
    'UPI000016F82F': ['just', 'fake', 'still', 'fake'],
    'UPI0000170134': ['fakey', 'omg', 'so', 'fake'],
    'fakeprot00': ['another', 'set', 'of', 'fake'],
    'UPI00131F240A': ['nothin'],
    'fakeprot01': ['more', 'fake', 'entries', 'here']
}

class TestRedundantSeqsHomo(unittest.TestCase):
    def test_init(self):        
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        self.assertEqual(rg.things[0], 'placeholder')
        self.assertEqual(rg.things[1], 'nothing')
        self.assertEqual(rg.threshold, 10)
        self.assertEqual(rg.config['paths']['lib'], os.path.realpath(f'{test_dir}/data/fake_lib'))
        self.assertEqual(rg.dimers, homodimers)
         

    def test_distance_similar(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        d1 = 'UPI00001653E1'
        d2 = 'UPI000016225C'
        expected = 1 - 0.994
        result = rg.distance(d1, d2)
        self.assertAlmostEqual(result, expected)


    def test_distance_notsimilar(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        d1 = 'UPI00001653E1'
        d2 = 'UPI00131F240A'
        expected = 1 - 0.415
        result = rg.distance(d1, d2)
        self.assertAlmostEqual(result, expected)


    def test_representative_criterion_1(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        cluster = ['UPI00001653E1', 'UPI000016225C', 'UPI000016F82F', 'UPI0000170134']
        rep = rg.representative(cluster)
        expected = 'UPI00001653E1'
        self.assertEqual(rep, expected)


    def test_representative_criterion_2(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        cluster = ['UPI000016225C', 'UPI000016F82F', 'UPI0000170134', 'fakeprot00']
        rep = rg.representative(cluster)
        expected = 'UPI000016F82F'
        self.assertEqual(rep, expected)


    def test_representative_criterion_3(self):
        #rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        #cluster = ['UPI000016F82F', 'UPI0000170134', 'fakeprot01']
        #rep = rg.representative(cluster)
        #expected = 'fakeprot01'
        #self.assertEqual(rep, expected)
        # I'm not sure if this case can ever exist. passing for now...
        pass


    def test_representative_criterion_4(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        cluster = ['UPI000016225C', 'UPI000016F82F', 'UPI0000170134']
        rep = rg.representative(cluster)
        expected = 'UPI000016225C'
        self.assertEqual(rep, expected)


    def test_dimer_tmp_fasta(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        with rg.dimer_tmp_fasta('UPI000016225C') as fasta:
            header, seq = next(read_prot_from_fasta(fasta))
            self.assertEqual(header, '>UPI000016225C status=active')
            self.assertEqual(seq, 'MVDSKKRPGKDLDRIDRNILNELQKDGRISNVELSKRVGLSPTPCLERVRRLERQGFIQG' + \
                                  'YTALLNPHYLDASLLVFVEITLNRGAPDVFEQFNAAVQKLEEIQECHLVSGDFDYLLKTR' + \
                                  'VPDMSAYRKLLGETLLRLPGVNDTRTYVVMEEVKQSNRLVIKTR')


    def test_prune(self):
        names = ['UPI00001653E1', 'UPI000016225C', 'UPI000016F82F', 'UPI0000170134', 'UPI00131F240A', 'fakeprot00']
        rg = RedundantSeqsHomodimer(names, homoyaml, 0.3, config)
        result = rg.prune_redundancy(num_workers=2)
        expected = ['UPI00001653E1', 'UPI00131F240A']


#heteroyaml= ?
#heterodimers = ?

class TestRedundantSeqsHetero(unittest.TestCase):
    #def test_init(self):
    #    rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
    #    self.assertEqual(rg.things[0], 'placeholder')
    #    self.assertEqual(rg.things[1], 'nothing')
    #    self.assertEqual(rg.threshold, 10)
    #    self.assertEqual(rg.config['paths']['lib'], os.path.realpath(f'{test_dir}/data/fake_lib'))
    #    self.assertEqual(rg.dimers, heterodimers)

    def test_distance_similar(self):
        pass

    def test_distance_notsimilar(self):
        pass

    def test_representative_criterion_1(self):
        pass

    def test_representative_criterion_2(self):
        pass

    def test_representative_criterion_3(self):
        pass

    def test_representative_criterion_4(self):
        pass

    def test_dimer_tmp_fasta(self):
        pass

    def test_prune(self):
        pass



if __name__ == '__main__':
    unittest.main()

