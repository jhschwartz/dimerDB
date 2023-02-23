import unittest
import contextlib
import os

import sys
sys.path.append('..')

from unredundant import RedundantSeqs, RedundantSeqsHomodimer
from read_fasta import read_prot_from_fasta

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
test_data = f'{test_dir}/data/unredundant_dimer_seqs'


nw = f'{test_dir}/../../bin/NWalign/align'
config = {
    'paths': {
        'nwalign': nw,
        #'intermediates_homodimer_filtering': f'{test_data}/homodimer_filtering'
        'lib': f'{test_data}/lib'
    }
}


class ConcreteRedundantSeqs(RedundantSeqs):
    def __init__(self, names, pklfile, threshold, config):
        super().__init__(names, pklfile, threshold, config)
    
    @staticmethod
    def __get_fasta_for_test(name_not_dimer):
        fastas_for_testing = {
            'O15205':  f'{test_data}/redundant_monomer_fastas/O15205.fasta', 
            'P63072':  f'{test_data}/redundant_monomer_fastas/P63072.fasta', 
            'Q921A3':  f'{test_data}/redundant_monomer_fastas/Q921A3.fasta', 
            'P0DTC2':  f'{test_data}/redundant_monomer_fastas/P0DTC2.fasta', 
            'P11223':  f'{test_data}/redundant_monomer_fastas/P11223.fasta', 
            'P59594':  f'{test_data}/redundant_monomer_fastas/P59594.fasta'  
        }
        return fastas_for_testing[name_not_dimer]
        
    def distance(self, d1name, d2name):
        fasta1 = self.__get_fasta_for_test(d1name)
        fasta2 = self.__get_fasta_for_test(d2name)
        return 1 - super().max_both_ways_nw(self.config['paths']['nwalign'], fasta1, fasta2)


class TestRedundantSeqsSubclass(unittest.TestCase):
    def test_redundant_seqs_subclass(self):
        names = ['O15205', 'P63072', 'Q921A3', 'P0DTC2', 'P11223', 'P59594']
        pklfile = f'{test_data}/redundant_monomer_fastas/fake_seqs.pkl'
        threshold = 0.5 
        redundant_seqs = ConcreteRedundantSeqs(names, pklfile, threshold, config)
        result = redundant_seqs.prune_redundancy(num_workers=2)
        expected = ['Q921A3', 'P11223', 'P59594']



    
    def test_max_both_ways_nw(self):
        fasta1 = f'{test_data}/redundant_monomer_fastas/O15205.fasta'
        fasta2 = f'{test_data}/redundant_monomer_fastas/P63072.fasta'
        expected = 0.71
        result1 = ConcreteRedundantSeqs.max_both_ways_nw(nw, fasta1, fasta2)
        result2 = ConcreteRedundantSeqs.max_both_ways_nw(nw, fasta2, fasta1)
        self.assertEqual(result1, expected)
        self.assertEqual(result1, result2)



homopkl = f'{test_data}/homodimers.pkl'
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
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homopkl, 10, config)
        self.assertEqual(rg.things[0], 'placeholder')
        self.assertEqual(rg.things[1], 'nothing')
        self.assertEqual(rg.threshold, 10)
        self.assertEqual(rg.dimers, homodimers)
         

    def test_distance_similar(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homopkl, 10, config)
        d1 = 'UPI00001653E1'
        d2 = 'UPI000016225C'
        expected = 1 - 0.994
        result = rg.distance(d1, d2)
        self.assertAlmostEqual(result, expected)


    def test_distance_notsimilar(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homopkl, 10, config)
        d1 = 'UPI00001653E1'
        d2 = 'UPI00131F240A'
        expected = 1 - 0.415
        result = rg.distance(d1, d2)
        self.assertAlmostEqual(result, expected)


    def test_representative_criterion_1(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homopkl, 10, config)
        cluster = ['UPI00001653E1', 'UPI000016225C', 'UPI000016F82F', 'UPI0000170134']
        rep = rg.representative(cluster)
        expected = 'UPI00001653E1'
        self.assertEqual(rep, expected)


    def test_representative_criterion_2(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homopkl, 10, config)
        cluster = ['UPI000016225C', 'UPI000016F82F', 'UPI0000170134', 'fakeprot00']
        rep = rg.representative(cluster)
        expected = 'UPI000016F82F'
        self.assertEqual(rep, expected)


    def test_representative_criterion_3(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homopkl, 10, config)
        cluster = ['fakeprot00', 'fakeprot01']
        rep = rg.representative(cluster)
        expected = 'fakeprot01'
        self.assertEqual(rep, expected)


    def test_representative_criterion_4(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homopkl, 10, config)
        cluster = ['UPI000016225C', 'UPI000016F82F', 'UPI0000170134']
        rep = rg.representative(cluster)
        expected = 'UPI000016225C'
        self.assertEqual(rep, expected)


    def test_prune(self):
        names = ['UPI00001653E1', 'UPI000016225C', 'UPI000016F82F', 'UPI0000170134', 'UPI00131F240A', 'fakeprot00']
        rg = RedundantSeqsHomodimer(names, homopkl, 0.3, config)
        result = rg.prune_redundancy(num_workers=2)
        expected = ['UPI00001653E1', 'UPI00131F240A']





if __name__ == '__main__':
    unittest.main()

