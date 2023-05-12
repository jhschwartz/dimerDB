import unittest
import os
import subprocess
import shutil

from run_tmp_snakemake import run_tmp_snakemake


import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = os.path.join(test_dir, 'data', '4_prep_final_data')

snakefile = os.path.join('..', '4_prep_final_data.smk')
snake_exe = '/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/snakemake'
config = os.path.abspath('config_for_test_4.yaml')
outname = 'tmpout'   
scripts = os.path.abspath('../../scripts')

class TestSubworkflow4Rules(unittest.TestCase):
    def test_rule_out_pdb_div_tar(self):
        rule = 'all' # cannot run rule alone because of wilcards
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            a6 = os.path.join(outdir, 'pdb', 'div', 'a6.tar.gz')
            ql = os.path.join(outdir, 'pdb', 'div', 'ql.tar.gz')
            self.assertTrue(os.path.exists(a6))
            self.assertTrue(os.path.exists(ql))


    def test_rule_out_extra_clusters(self):
        assert False
        rule = 'all' # cannot run rule alone because of wilcards
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            

    def test_rule_out_pdb_one_fasta(self):
        rule = 'out_pdb_one_fasta'
        expected_file = os.path.join(test_data, 'out_extra_clusters', 'expected.fasta')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            out_fasta = os.path.join(outdir, out_fasta)
            with open(out_fasta, 'r') as f:
                result_text = f.read()
            with open(expected_file, 'r') as f:
                expected_text = f.read()
            self.assertEqual(result_text, expected_text)


    def test_rule_out_all_dimers_chains(self):
        assert False
        rule = 'out_all_dimers_chains'
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)

    
    def test_rule_out_all_seqs(self):
        assert False
        run = 'out_all_seqs'
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
    
    
    def test_rule_out_cluster_info(self):
        assert False
        rule = 'out_cluster_info'
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
    
    
    def test_rule_out_nonredundant_dimers_chains(self):
        assert False
        rule = 'out_nonredundant_dimers_chains'
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
    
    
    def test_rule_out_nonredundant_seqs(self):
        assert False
        rule = 'out_nonredundant_seqs'
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)


    def test_rule_out_summarize(self):
        assert False
        rule = 'out_summarize'
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, rule) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)


