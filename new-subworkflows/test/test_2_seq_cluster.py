import unittest
import os
import filecmp

from run_tmp_snakemake import run_tmp_snakemake

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = os.path.join(test_dir, 'data', '2_seq_cluster')


snakefile = os.path.join('..', '2_seq_cluster.smk')
snake_exe = '/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/snakemake'
config = os.path.abspath('config_for_test_2.yaml')
outint = 'intermediates'
outlib = 'lib'
inlib = os.path.abspath(os.path.join(data_dir, 'lib'))
scripts = os.path.abspath('../../scripts')
bin = os.path.abspath('../../bin')



class TestSubworkflow2Rules(unittest.TestCase):
    def test_process_fastas(self):
        rule = 'process_fastas'
        
        expected_all_seqs = os.path.join(data_dir, 'expected_lib', 'all_seqs.fasta')
        expected_fasta_index = os.path.join(data_dir, 'expected_lib', 'fasta_index.tsv')
        expected_individual_fastas = []
        with open(expected_fasta_index, 'r') as f:
            for line in f:
                fasta = os.path.join(data_dir, 'expected_lib', 'fasta', line.split()[0])
                expected_individual_fastas.append(fasta)
        
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, inlib) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outlib)

            result_all_seqs = os.path.join(outdir, 'all_seqs.fasta')
            result_fasta_index = os.path.join(outdir, 'fasta_index.tsv')

            self.assertTrue(filecmp.cmp(result_all_seqs, expected_all_seqs, shallow=False))
            self.assertTrue(filecmp.cmp(result_fasta_index, expected_fasta_index, shallow=False))

            result_individual_fastas = []
            with open(result_fasta_index, 'r') as f:
                for line in f:
                    fasta = os.path.join(outdir, 'fasta', line.split()[0])
                    expected_individual_fastas.append(fasta)
            
            for expected, result in zip(expected_individual_fastas, result_individual_fastas):
                self.assertTrue(filecmp.cmp(expected, result, shallow=False))



    def test_write_homodimers_chains_fastas(self):
        rule = 'write_homodimers_chains_fastas'
        expected_fasta = os.path.join(data_dir, 'expected_intermediates', 'homodimer_seqs.fasta')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, inlib) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outint)
            result_fasta = os.path.join(outdir, 'homodimer_seqs.fasta')
            
            self.assertTrue(filecmp.cmp(result_fasta, expected_fasta, shallow=False))
    



    def test_run_seq_cluster(self):
        rule = 'run_seq_cluster'
        clust1 = '6qle-a1-m1-cN-42'
        clust2 = '7a6w-a1-m1-cAAA'
        expected_cluster_index = os.path.join(data_dir, 'expected_intermediates', 
                                    'cluster', 'cluster_index.txt')
        expected_clust1_chains = os.path.join(data_dir, 'expected_intermediates',
                                    'cluster', 'seq_clusters', clust1, 'members_chains.txt')
        expected_clust1_dimers = os.path.join(data_dir, 'expected_intermediates',
                                    'cluster', 'seq_clusters', clust1, 'members_dimers.txt')
        expected_clust2_chains = os.path.join(data_dir, 'expected_intermediates',
                                    'cluster', 'seq_clusters', clust2, 'members_chains.txt')
        expected_clust2_dimers = os.path.join(data_dir, 'expected_intermediates',
                                    'cluster', 'seq_clusters', clust2, 'members_dimers.txt') 

        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, inlib) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outint)
            
            result_cluster_index = os.path.join(outdir, 'cluster', 'cluster_index.txt')
            result_clust1_chains = os.path.join(outdir, 'cluster', 
                                    'seq_clusters', clust1, 'members_chains.txt')
            result_clust1_dimers = os.path.join(outdir, 'cluster',
                                    'seq_clusters', clust1, 'members_dimers.txt')
            result_clust2_chains = os.path.join(outdir, 'cluster',
                                    'seq_clusters', clust2, 'members_chains.txt')
            result_clust2_dimers = os.path.join(outdir, 'cluster',
                                    'seq_clusters', clust2, 'members_dimers.txt')


            self.assertTrue(filecmp.cmp(expected_cluster_index, result_cluster_index, shallow=False))
            self.assertTrue(filecmp.cmp(expected_clust1_chains, result_clust1_chains, shallow=False))
            self.assertTrue(filecmp.cmp(expected_clust1_dimers, result_clust1_dimers, shallow=False))
            self.assertTrue(filecmp.cmp(expected_clust2_chains, result_clust2_chains, shallow=False))
            self.assertTrue(filecmp.cmp(expected_clust2_dimers, result_clust2_dimers, shallow=False))


