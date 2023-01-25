import unittest
import contextlib
import os

import sys
sys.path.append('..')

from unredundant import RedundantSeqs, RedundantSeqsHomodimer, RedundantSeqsHeterodimer
from read_fasta import read_prot_from_fasta

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()


nw = f'{test_dir}/../../bin/NWalign/align'
config = {
    'paths': {
        'lib': f'{test_dir}/data/fake_lib',
        'nwalign': nw,
        'intermediates_homodimer_filtering': f'{test_dir}/data/unred_seq_fakedata/fake_homo_filter',
        'intermediates_heterodimer_filtering': f'{test_dir}/data/unred_seq_fakedata/fake_hetero_filter' 
    }
}


class ConcreteRedundantSeqs(RedundantSeqs):
    def __init__(self, names, yamlfile, threshold, config):
        super().__init__(names, yamlfile, threshold, config)
    
    @staticmethod
    def __get_fasta_for_test(name_not_dimer):
        fastas_for_testing = {
            'O15205':  f'{test_dir}/data/redundant_monomer_fastas/O15205.fasta', 
            'P63072':  f'{test_dir}/data/redundant_monomer_fastas/P63072.fasta', 
            'Q921A3':  f'{test_dir}/data/redundant_monomer_fastas/Q921A3.fasta', 
            'P0DTC2':  f'{test_dir}/data/redundant_monomer_fastas/P0DTC2.fasta', 
            'P11223':  f'{test_dir}/data/redundant_monomer_fastas/P11223.fasta', 
            'P59594':  f'{test_dir}/data/redundant_monomer_fastas/P59594.fasta'  
        }
        return fastas_for_testing[name_not_dimer]
        
    def distance(self, d1name, d2name):
        fasta1 = self.__get_fasta_for_test(d1name)
        fasta2 = self.__get_fasta_for_test(d2name)
        return 1 - super().max_both_ways_nw(self.config['paths']['nwalign'], fasta1, fasta2)


class TestRedundantSeqsSubclass(unittest.TestCase):
    def test_redundant_seqs_subclass(self):
        names = ['O15205', 'P63072', 'Q921A3', 'P0DTC2', 'P11223', 'P59594']
        yamlfile = f'{test_dir}/data/redundant_monomer_fastas/fake_seqs.yaml'
        threshold = 0.5 
        redundant_seqs = ConcreteRedundantSeqs(names, yamlfile, threshold, config)
        result = redundant_seqs.prune_redundancy(num_workers=2)
        expected = ['Q921A3', 'P11223', 'P59594']


    def test_calc_nw(self):
        fasta1 = f'{test_dir}/data/redundant_monomer_fastas/O15205.fasta'
        fasta2 = f'{test_dir}/data/redundant_monomer_fastas/P63072.fasta'
        expected1 = 0.71
        expected2 = 0.697
        result1 = ConcreteRedundantSeqs._calc_nw(nw, fasta1, fasta2)
        result2 = ConcreteRedundantSeqs._calc_nw(nw, fasta2, fasta1)
        self.assertEqual(result1, expected1)
        self.assertEqual(result2, expected2)

    
    def test_max_both_ways_nw(self):
        fasta1 = f'{test_dir}/data/redundant_monomer_fastas/O15205.fasta'
        fasta2 = f'{test_dir}/data/redundant_monomer_fastas/P63072.fasta'
        expected = 0.71
        result1 = ConcreteRedundantSeqs.max_both_ways_nw(nw, fasta1, fasta2)
        result2 = ConcreteRedundantSeqs.max_both_ways_nw(nw, fasta2, fasta1)
        self.assertEqual(result1, expected)
        self.assertEqual(result1, result2)



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
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        cluster = ['fakeprot00', 'fakeprot01']
        rep = rg.representative(cluster)
        expected = 'fakeprot01'
        self.assertEqual(rep, expected)


    def test_representative_criterion_4(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], homoyaml, 10, config)
        cluster = ['UPI000016225C', 'UPI000016F82F', 'UPI0000170134']
        rep = rg.representative(cluster)
        expected = 'UPI000016225C'
        self.assertEqual(rep, expected)


    def test_prune(self):
        names = ['UPI00001653E1', 'UPI000016225C', 'UPI000016F82F', 'UPI0000170134', 'UPI00131F240A', 'fakeprot00']
        rg = RedundantSeqsHomodimer(names, homoyaml, 0.3, config)
        result = rg.prune_redundancy(num_workers=2)
        expected = ['UPI00001653E1', 'UPI00131F240A']


heteroyaml = f'{test_dir}/data/unred_seq_fakedata/heterodimers.yaml'
heterodimers = {
    'UPI0000052DE7-UPI0000052FA5': ['3o8o_C-3o8o_D', '3o8o_B-3o8o_C'],
    'UPI0006C47D63-UPI0006BF9280': ['not', 'real'],
    'UPI000012DB37-UPI00003BB37D': ['nope'],
    'UPI000013D4D2-UPI0000135942': [1],
    'UPI000012E857-UPI0000135948': [1],
    'UPI000012E856-UPI0001E427C0': [1],
    'UPI000209B964-UPI0003CD15F4': [1],
    'UPI000013D4D2-UPI000018FE19': [0],
    'UPI000012E857-UPI00131F240A': [0],
    'UPI000012E856-UPI0003C96339': [0],
    'UPI0003C96339-UPI000012E856': [0],
    'fakex-dimerx': [0]
}

class TestRedundantSeqsHetero(unittest.TestCase):
    def test_init(self):
        rg = RedundantSeqsHomodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        self.assertEqual(rg.things[0], 'placeholder')
        self.assertEqual(rg.things[1], 'nothing')
        self.assertEqual(rg.threshold, 10)
        self.assertEqual(rg.config['paths']['lib'], os.path.realpath(f'{test_dir}/data/fake_lib'))
        self.assertEqual(rg.dimers, heterodimers)


    def test_distance_similar_forward(self):
        rg = RedundantSeqsHeterodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        d1 = 'UPI0000052DE7-UPI0000052FA5' # PFK1-PFK2, S. cerivisae
        d2 = 'UPI0006C47D63-UPI0006BF9280' # PKF1-PFK2, S. eubayanus
        expected = 1 - (0.976 + 0.970)/2
        result = rg.distance(d1, d2)
        self.assertAlmostEqual(result, expected)


    def test_distance_similar_backward(self):
        rg = RedundantSeqsHeterodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        d1 = 'UPI0000052DE7-UPI0000052FA5' # PFK1-PFK2, S. cerivisae
        d2 = 'UPI0006BF9280-UPI0006C47D63' # PKF2-PFK1, S. eubayanus
        expected = 1 - (0.976 + 0.970)/2
        result = rg.distance(d1, d2)
        self.assertAlmostEqual(result, expected)


    def test_distance_notsimilar(self):
        rg = RedundantSeqsHeterodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        d1 = 'UPI000012DB37-UPI00003BB37D' # PFK1-PFK2, K. lactis
        d2 = 'UPI000012E856-UPI0003C96339' # Lactase(rabbit)-Spike(Bat SARS)
        expected = 1 - (0.284 + 0.320)/2
        result = rg.distance(d1, d2)
        self.assertAlmostEqual(result, expected)
    
    
    def test_distance_halfsimilar(self):
        rg = RedundantSeqsHeterodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        d1 = 'UPI000209B964-UPI0003CD15F4' # SonicHedgehog(sheep)-Lactase(sheep)
        d2 = 'UPI0003C96339-UPI000012E856' # Spike(Bat SARS)-Lactase(rabbit)
        expected = 1 - (0.191 + 0.8)/2
        result = rg.distance(d1, d2)
        self.assertAlmostEqual(result, expected)


    def test_representative_criterion_1(self):
        rg = RedundantSeqsHeterodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        cluster = ['UPI0000052DE7-UPI0000052FA5', 'UPI0006C47D63-UPI0006BF9280', 'UPI000012DB37-UPI00003BB37D']
        rep = rg.representative(cluster)
        expected = 'UPI0000052DE7-UPI0000052FA5'
        self.assertEqual(rep, expected)


    def test_representative_criterion_2(self):
        rg = RedundantSeqsHeterodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        cluster = ['UPI000012E857-UPI0000135948', 'UPI000013D4D2-UPI0000135942', 'UPI000012E856-UPI0001E427C0', 'UPI000209B964-UPI0003CD15F4'] # lactase/sonic fake dimers
        rep = rg.representative(cluster)
        expected = 'UPI000013D4D2-UPI0000135942' # human
        self.assertEqual(rep, expected)


    def test_representative_criterion_3(self):
        rg = RedundantSeqsHeterodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        cluster = ['UPI000012E857-UPI0000135948', 'UPI000013D4D2-UPI0000135942'] # lactase/sonic fake dimers
        rep = rg.representative(cluster)
        expected = 'UPI000013D4D2-UPI0000135942' # human
        self.assertEqual(rep, expected)


    def test_representative_criterion_4(self):
        rg = RedundantSeqsHeterodimer(['placeholder', 'nothing'], heteroyaml, 10, config)
        cluster = ['UPI000013D4D2-UPI0000135942', 'fakex-dimerx']
        rep = rg.representative(cluster)
        expected = 'UPI000013D4D2-UPI0000135942'
        self.assertEqual(rep, expected)

    
    def test_cluster(self):
        dimers = ['UPI0000052DE7-UPI0000052FA5', 'UPI0006C47D63-UPI0006BF9280', 'UPI000012DB37-UPI00003BB37D', 'UPI000013D4D2-UPI0000135942', 'UPI000012E857-UPI0000135948', 'UPI000012E856-UPI0001E427C0', 'UPI000209B964-UPI0003CD15F4', 'UPI000013D4D2-UPI000018FE19', 'UPI000012E857-UPI00131F240A', 'UPI000012E856-UPI0003C96339', 'UPI0003C96339-UPI000012E856', 'fakex-dimerx']
        rg = RedundantSeqsHeterodimer(dimers, heteroyaml, 0.3, config)
        rg.initiate_distance_matrix(num_workers=8)
        num_clusters = rg.initiate_clusters()
        self.assertEqual(num_clusters, 3)
        cluster_pkf1_pfk2 = rg.cluster_index_of_thing('UPI0000052DE7-UPI0000052FA5')
        cluster_lastase_sonic = rg.cluster_index_of_thing('UPI000013D4D2-UPI0000135942')
        cluster_lactase_spike = rg.cluster_index_of_thing('UPI0003C96339-UPI000012E856')
        self.assertEqual(rg.retrieve_cluster(cluster_pkf1_pfk2), dimers[:3])
        self.assertEqual(rg.retrieve_cluster(cluster_lastase_sonic), dimers[3:7]+[dimers[-1]])
        self.assertEqual(rg.retrieve_cluster(cluster_lactase_spike), dimers[7:-1])


    def test_prune(self):
        dimers = ['UPI0000052DE7-UPI0000052FA5', 'UPI0006C47D63-UPI0006BF9280', 'UPI000012DB37-UPI00003BB37D', 'UPI000013D4D2-UPI0000135942', 'UPI000012E857-UPI0000135948', 'UPI000012E856-UPI0001E427C0', 'UPI000209B964-UPI0003CD15F4', 'UPI000013D4D2-UPI000018FE19', 'UPI000012E857-UPI00131F240A', 'UPI000012E856-UPI0003C96339', 'UPI0003C96339-UPI000012E856', 'fakex-dimerx']
        rg = RedundantSeqsHeterodimer(dimers, heteroyaml, 0.3, config)
        result = rg.prune_redundancy(num_workers=6)
        expected = ['UPI0000052DE7-UPI0000052FA5', 'UPI000013D4D2-UPI0000135942', 'UPI000012E857-UPI00131F240A']
        self.assertEqual(set(result), set(expected))



if __name__ == '__main__':
    unittest.main()

