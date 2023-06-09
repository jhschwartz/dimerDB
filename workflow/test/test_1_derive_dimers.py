import unittest
import os
import sys
import filecmp


import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = os.path.join(test_dir, 'data', '1_derive_dimers')

sys.path.append(test_dir)
from run_tmp_snakemake import run_tmp_snakemake

conda_env = os.environ['conda_env']


snakefile = os.path.join('..', '1_derive_dimers.smk')
snake_exe = os.path.join(conda_env, 'bin/snakemake')
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

        expected_ql_chains = os.path.join(data_dir, 'all', 'expected_ql_chains.txt')
        expected_a6_chains = os.path.join(data_dir, 'all', 'expected_a6_chains.txt')
        expected_ql_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_ql.txt')
        expected_a6_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_a6.txt')
        expected_ql_contacting = os.path.join(data_dir, 'all', 'expected_ql_contacts.txt')
        expected_a6_contacting = os.path.join(data_dir, 'all', 'expected_a6_contacts.txt') 
        expected_ql_ids = os.path.join(data_dir, 'all', 'expected_ql_ids.tsv')
        expected_a6_ids = os.path.join(data_dir, 'all', 'expected_a6_ids.tsv')

        expected_all_homodimers = os.path.join(data_dir, 'all', 'expected_all_homodimers.txt') 

        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib) as tmpsnake:
            outinter = os.path.join(os.path.dirname(tmpsnake), 'intermediates')
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)

            result_ql_chains = os.path.join(outdir, 'ql', 'chains.txt')
            result_a6_chains = os.path.join(outdir, 'a6', 'chains.txt')
            self.assertTrue(filecmp.cmp(result_ql_chains, expected_ql_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_chains, expected_a6_chains, shallow=False))
           
            result_ql_poss = os.path.join(outdir, 'ql', 'intra_assembly_chain_pairs.txt')
            result_a6_poss = os.path.join(outdir, 'a6', 'intra_assembly_chain_pairs.txt')
            self.assertTrue(filecmp.cmp(result_ql_poss, expected_ql_poss_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_poss, expected_a6_poss_chains, shallow=False))

            result_ql_cont = os.path.join(outdir, 'ql', 'dimers.txt')
            result_a6_cont = os.path.join(outdir, 'a6', 'dimers.txt')
            self.assertTrue(filecmp.cmp(result_ql_cont, expected_ql_contacting, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_cont, expected_a6_contacting, shallow=False))


            result_ql_ids = os.path.join(outdir, 'ql', 'dimer_seq_ids.tsv')
            result_a6_ids = os.path.join(outdir, 'a6', 'dimer_seq_ids.tsv')
            self.assertTrue(filecmp.cmp(result_ql_ids, expected_ql_ids, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_ids, expected_a6_ids, shallow=False))

            result_all_homodimers = os.path.join(outinter, 'all_homodimers.txt')
            self.assertTrue(filecmp.cmp(result_all_homodimers, expected_all_homodimers, shallow=False))


