import unittest
import os
import subprocess
import shutil
import filecmp
from datetime import datetime
import tempfile

from run_tmp_snakemake import run_tmp_snakemake

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = os.path.join(test_dir, 'data', '4_prep_final_data')

snakefile = os.path.join('..', '4_prep_final_data.smk')
snake_exe = '/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/snakemake'
config = os.path.abspath('config_for_test_4.yaml')
outname = 'tmpout'   
scripts = os.path.abspath('../../scripts')
intermediates = os.path.join(data_dir, 'shared', 'intermediates')
lib = os.path.join(data_dir, 'shared', 'lib')
bin = os.path.abspath('../../bin')

class TestSubworkflow4Rules(unittest.TestCase):
    def test_rule_out_pdb_div_tar(self):
        rule = 'all' # cannot run rule alone because of wilcards
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            a6 = os.path.join(outdir, 'pdb', 'div', 'a6.tar.gz')
            ql = os.path.join(outdir, 'pdb', 'div', 'ql.tar.gz')
            self.assertTrue(os.path.exists(a6))
            self.assertTrue(os.path.exists(ql))


    def test_rule_out_extra_clusters(self):
        rule = 'all' # cannot run rule alone because of wilcards
        expected_fasta_80 = os.path.join(data_dir, 'out_extra_clusters', 'expected80_seqs.fasta')
        expected_chains_80 = os.path.join(data_dir, 'out_extra_clusters', 'expected80_chains.txt')
        expected_dimers_80 = os.path.join(data_dir, 'out_extra_clusters', 'expected80_dimers.txt')
        expected_fasta_30 = os.path.join(data_dir, 'out_extra_clusters', 'expected30_seqs.fasta')
        expected_chains_30 = os.path.join(data_dir, 'out_extra_clusters', 'expected30_chains.txt')
        expected_dimers_30 = os.path.join(data_dir, 'out_extra_clusters', 'expected30_dimers.txt')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_fasta_80 = os.path.join(outdir, 'extra', 'cluster80_seqs.fasta')
            result_chains_80 = os.path.join(outdir, 'extra', 'cluster80_chains.txt') 
            result_dimers_80 = os.path.join(outdir, 'extra', 'cluster80_dimers.txt')  
            result_fasta_30 = os.path.join(outdir, 'extra', 'cluster30_seqs.fasta')
            result_chains_30 = os.path.join(outdir, 'extra', 'cluster30_chains.txt') 
            result_dimers_30 = os.path.join(outdir, 'extra', 'cluster30_dimers.txt')  
            
            self.assertTrue(filecmp.cmp(expected_fasta_80, result_fasta_80, shallow=False))
            self.assertTrue(filecmp.cmp(expected_chains_80, result_chains_80, shallow=False))
            self.assertTrue(filecmp.cmp(expected_dimers_80, result_dimers_80, shallow=False))
            self.assertTrue(filecmp.cmp(expected_fasta_30, result_fasta_30, shallow=False))
            self.assertTrue(filecmp.cmp(expected_chains_30, result_chains_30, shallow=False))
            self.assertTrue(filecmp.cmp(expected_dimers_30, result_dimers_30, shallow=False))
            



    def test_rule_out_all_seqs(self):
        rule = 'out_all_seqs'
        expected_file = os.path.join(data_dir, 'out_all_seqs', 'expected.fasta')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_fasta = os.path.join(outdir, 'all', 'seqs.fasta')
            self.assertTrue(filecmp.cmp(expected_file, result_fasta, shallow=False))



    def test_rule_out_all_dimers_chains(self):
        rule = 'out_all_dimers_chains'
        expected_dimers_file = os.path.join(data_dir, 
                                'out_all_dimers_chains', 'expected_dimers.txt')
        expected_chains_file = os.path.join(data_dir, 
                                'out_all_dimers_chains', 'expected_chains.txt')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_dimers_file = os.path.join(outdir, 'all', 'dimers.txt') 
            result_chains_file = os.path.join(outdir, 'all', 'chains.txt')
            
            self.assertTrue(filecmp.cmp(expected_dimers_file, 
                                result_dimers_file, shallow=False))
            self.assertTrue(filecmp.cmp(expected_chains_file, 
                                result_chains_file, shallow=False))



    def test_rule_out_cluster_info(self):
        rule = 'out_cluster_info'
        expected_clusters_file = os.path.join(data_dir, 'out_cluster_info', 'expected_clusters.txt')
        expected_membership_file = os.path.join(data_dir, 'out_cluster_info', 'expected_membership.tsv')

        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_clusters_file = os.path.join(outdir, 'cluster', 'clusters.txt') 
            result_membership_file = os.path.join(outdir, 'cluster', 'membership.tsv')

            self.assertTrue(filecmp.cmp(result_clusters_file, expected_clusters_file, shallow=False))
            self.assertTrue(filecmp.cmp(result_membership_file, expected_membership_file, shallow=False))

   


    def test_rule_out_nonredundant_dimers_chains(self):
        rule = 'out_nonredundant_dimers_chains'
        expected_dimers_file = os.path.join(data_dir, 
                            'out_nonredundant_dimers_chains', 'expected_dimers.txt')
        expected_chains_file = os.path.join(data_dir, 
                            'out_nonredundant_dimers_chains', 'expected_chains.txt')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_dimers_file = os.path.join(outdir, 'nonredundant', 'dimers.txt') 
            result_chains_file = os.path.join(outdir, 'nonredundant', 'chains.txt')
            self.assertTrue(filecmp.cmp(expected_dimers_file, 
                                result_dimers_file, shallow=False))
            self.assertTrue(filecmp.cmp(expected_chains_file, 
                                result_chains_file, shallow=False))
            


    
    def test_rule_out_nonredundant_seqs(self):
        rule = 'out_nonredundant_seqs'
        expected_seqs_file = os.path.join(data_dir, 
                            'out_nonredundant_seqs', 'expected.fasta')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_seqs_file = os.path.join(outdir, 'nonredundant', 'seqs.fasta')
            self.assertTrue(filecmp.cmp(expected_seqs_file,
                                result_seqs_file, shallow=False))



    def test_rule_out_summary(self):
        rule = 'out_summary'
        time_before_run = datetime.utcnow()
        expected_summary_file = os.path.join(data_dir, 'out_summary', 'expected.txt')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_summary_file = os.path.join(outdir, 'summary.txt')
            with open(expected_summary_file, 'r') as ef, \
              open(result_summary_file, 'r') as rf:

                # check the first line, keyword and timestamp
                expected_time_line = ef.readline()
                result_time_line = rf.readline()
                self.assertEqual(expected_time_line.split()[0], 
                                 result_time_line.split()[0]
                                 )
                spl = result_time_line.split()
                time_of_run_iso = f'{spl[1]} {spl[2]}'
                time_of_run = datetime.fromisoformat(time_of_run_iso)
                time_after_run = datetime.utcnow()
                self.assertTrue(time_of_run > time_before_run)
                self.assertTrue(time_of_run < time_after_run)
                
                # check counts in remaining lines, using tempfiles
                with tempfile.NamedTemporaryFile('w+t') as e_temp, \
                  tempfile.NamedTemporaryFile('w+t') as r_temp:
                
                    for eline in ef:
                        e_temp.write(eline)
                    for rline in rf:
                        r_temp.write(rline)

                    e_temp.seek(0)
                    r_temp.seek(0)
                    self.assertTrue(filecmp.cmp(e_temp.name, r_temp.name, shallow=False))



    def test_bulk_count_contacts(self):
        rule = 'bulk_count_contacts'
        expected_counts_file = os.path.join(data_dir, 'bulk_count_contacts', 'expected.tsv')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_counts_file = os.path.join(outdir, 'info', 'contact_counts.tsv')
            self.assertTrue(filecmp.cmp(expected_counts_file, result_counts_file, shallow=False))



    def test_bulk_calc_seqid(self):
        rule = 'bulk_calc_seqid'
        expected_seqid_file = os.path.join(data_dir, 'bulk_calc_seqid', 'seq_ids.tsv')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_seqid_file = os.path.join(outdir, 'info', 'seq_ids.tsv')
            self.assertTrue(filecmp.cmp(expected_seqid_file, result_seqid_file, shallow=False))



    def test_out_dimers_info(self):
        rule = 'out_dimers_info'
        expected_seqid_file = os.path.join(data_dir, 'out_dimers_info', 'expected.tsv')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_seqid_file = os.path.join(outdir, 'nonredundant', 'dimers_info.tsv')
            self.assertTrue(filecmp.cmp(expected_seqid_file, result_seqid_file, shallow=False))





