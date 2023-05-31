import unittest
import os
import sys
import filecmp


import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = os.path.join(test_dir, 'data', '1_derive_dimers')

sys.path.append(test_dir)
from run_tmp_snakemake import run_tmp_snakemake


snakefile = os.path.join('..', '1_derive_dimers.smk')
snake_exe = '/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/snakemake'
config = os.path.abspath('config_for_test_1.yaml')
outname = 'intermediates/check_pairs'
scripts = os.path.abspath('../../scripts')
lib = os.path.join(data_dir, 'shared', 'lib')
bin = os.path.abspath('../../bin')



class TestSubworkflow1Rules(unittest.TestCase):
    
    # rules pair_possible_chains, check_contacts, categorize_dimers use samples
    # so cannot be run separately; thus we run subworkflow 1 end-to-end and check
    # the results of these rules.
    def test_rule_all(self):
        rule = 'all'

        expected_6qle_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_6qle.txt')
        expected_7a6w_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_7a6w.txt')
        expected_6qle_contacting = os.path.join(data_dir, 'all', 'expected_6qle_contacts.txt')
        expected_7a6w_contacting = os.path.join(data_dir, 'all', 'expected_7a6w_contacts.txt') 
        expected_6qle_homodimers = os.path.join(data_dir, 'all', 'expected_6qle_homodimers.txt')
        expected_7a6w_homodimers = os.path.join(data_dir, 'all', 'expected_7a6w_homodimers.txt')
        expected_all_homodimers = os.path.join(data_dir, 'all', 'expected_all_homodimers.txt') 

        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)

            result_6qle_poss = os.path.join(outdir, 'div', 'ql', '6qle',
                                                'all_intra_assembly_chain_pairs.txt')
            result_7a6w_poss = os.path.join(outdir, 'div', 'a6', '7a6w', 
                                                'all_intra_assembly_chain_pairs.txt')
            self.assertTrue(filecmp.cmp(result_6qle_poss, expected_6qle_poss_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_7a6w_poss, expected_7a6w_poss_chains, shallow=False))

            result_6qle_cont = os.path.join(outdir, 'div', 'ql', '6qle',
                                                'contacting_chain_pairs.txt')
            result_7a6w_cont = os.path.join(outdir, 'div', 'a6', '7a6w', 
                                                'contacting_chain_pairs.txt')
            self.assertTrue(filecmp.cmp(result_6qle_cont, expected_6qle_contacting, shallow=False))
            self.assertTrue(filecmp.cmp(result_7a6w_cont, expected_7a6w_contacting, shallow=False))

            result_6qle_homodimer = os.path.join(outdir, 'div', 'ql', '6qle',
                                                'homodimers.txt')
            result_7a6w_homodimer = os.path.join(outdir, 'div', 'a6', '7a6w', 
                                                'homodimers.txt')
            self.assertTrue(filecmp.cmp(result_6qle_homodimer, expected_6qle_homodimers, shallow=False))
            self.assertTrue(filecmp.cmp(result_7a6w_homodimer, expected_7a6w_homodimers, shallow=False))

            result_all_homodimers = os.path.join(outdir, 'all_homodimers.txt')
            self.assertTrue(filecmp.cmp(result_all_homodimers, expected_all_homodimers, shallow=False))


    def test_rule_write_chains(self):
        rule = 'write_chains'
        expected_ql = os.path.join(data_dir, 'write_chains', 'expected_6qle.txt')
        expected_a6 = os.path.join(data_dir, 'write_chains', 'expected_7a6w.txt')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib) as tmpsnake:
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            result_ql = os.path.join(outdir, 'div', 'ql', '6qle', 'chains.txt')
            result_a6 = os.path.join(outdir, 'div', 'a6', '7a6w', 'chains.txt')
            self.assertTrue(filecmp.cmp(result_ql, expected_ql, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6, expected_a6, shallow=False))

