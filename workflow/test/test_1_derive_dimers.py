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
outname = 'intermediates/check_pairs'
scripts = os.path.abspath('../../scripts')
lib = os.path.join(data_dir, 'shared', 'lib')
bin = os.path.abspath('../../bin')

sys.path.append(scripts)
from generate_rcsb_index import generate_rcsb_index


class TestSubworkflow1Rules(unittest.TestCase):
    def setUp(self):
        # recreate pdb_index for each test.
        self.indexfile = os.path.join(lib, 'pdb_index.txt')
        self.rcsbdir = os.path.join(lib, 'rcsb')
        generate_rcsb_index(self.rcsbdir, self.indexfile)


    def tearDown(self):
        os.remove(self.indexfile)

    
    # rules pair_possible_chains, check_contacts, categorize_dimers use samples
    # so cannot be run separately; thus we run subworkflow 1 end-to-end and check
    # the results of these rules.
    def test_rule_all_nodbs(self):
        rule = 'all'
        config = os.path.abspath('config_for_test_1_nodbs.yaml')

        expected_ql_chains = os.path.join(data_dir, 'all', 'expected_ql_chains.txt')
        expected_a6_chains = os.path.join(data_dir, 'all', 'expected_a6_chains.txt')
        expected_ql_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_ql.txt')
        expected_a6_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_a6.txt')
        expected_ql_contacting = os.path.join(data_dir, 'all', 'expected_ql_contacts.txt')
        expected_a6_contacting = os.path.join(data_dir, 'all', 'expected_a6_contacts.txt') 
        expected_ql_ids = os.path.join(data_dir, 'all', 'expected_ql_ids.tsv')
        expected_a6_ids = os.path.join(data_dir, 'all', 'expected_a6_ids.tsv')
        expected_db_ql = os.path.join(data_dir, 'shared', 'lib', 'contactsdb_full', 'ql.tsv')
        expected_db_a6 = os.path.join(data_dir, 'shared', 'lib', 'contactsdb_full', 'a6.tsv')
        

        expected_all_homodimers = os.path.join(data_dir, 'all', 'expected_all_homodimers.txt') 

        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib) as tmpsnake:
            outinter = os.path.join(os.path.dirname(tmpsnake), 'intermediates')
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            outlib = os.path.join(os.path.dirname(tmpsnake), 'lib')

            result_ql_chains = os.path.join(outdir, 'ql', 'chains.txt')
            result_a6_chains = os.path.join(outdir, 'a6', 'chains.txt')
            self.assertTrue(filecmp.cmp(result_ql_chains, expected_ql_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_chains, expected_a6_chains, shallow=False))
           
            result_ql_poss = os.path.join(outdir, 'ql', 'intra_assembly_chain_pairs.txt')
            result_a6_poss = os.path.join(outdir, 'a6', 'intra_assembly_chain_pairs.txt')
            self.assertTrue(filecmp.cmp(result_ql_poss, expected_ql_poss_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_poss, expected_a6_poss_chains, shallow=False))

            # empty lookup expected
            result_ql_lookup = os.path.join(outdir, 'ql', 'contacts_lookup.txt')
            result_a6_lookup = os.path.join(outdir, 'a6', 'contacts_lookup.txt')
            self.assertEqual(os.stat(result_ql_lookup).st_size, 0)
            self.assertEqual(os.stat(result_a6_lookup).st_size, 0)

            # we expect contacts_needcalc to be the same as the possible pairs file, but ordered incorrectly.
            result_ql_needcalc = os.path.join(outdir, 'ql', 'contacts_needcalc.txt')
            result_a6_needcalc = os.path.join(outdir, 'a6', 'contacts_needcalc.txt')
            with open(result_ql_needcalc) as rf, open(expected_ql_poss_chains) as ef:
                data_r = rf.readlines()
                data_e = ef.readlines()
                self.assertEqual(sorted(data_r), sorted(data_e))
            with open(result_a6_needcalc) as rf, open(expected_a6_poss_chains) as ef:
                data_r = rf.readlines()
                data_e = ef.readlines()
                self.assertEqual(sorted(data_r), sorted(data_e))

            # because no contactslib, dimers and contacts_calc are same, although contacts_calc is out of order so can't use filecmp.
            result_ql_cont = os.path.join(outdir, 'ql', 'contacts_calc.txt')
            result_a6_cont = os.path.join(outdir, 'a6', 'contacts_calc.txt')
            with open(result_ql_cont) as rf, open(expected_ql_contacting) as ef:
                data_r = rf.readlines()
                data_e = ef.readlines()
                self.assertEqual(sorted(data_r), sorted(data_e))
            with open(result_a6_cont) as rf, open(expected_a6_contacting) as ef:
                data_r = rf.readlines()
                data_e = ef.readlines()
                self.assertEqual(sorted(data_r), sorted(data_e))
            result_ql_dimers = os.path.join(outdir, 'ql', 'dimers.txt')
            result_a6_dimers = os.path.join(outdir, 'a6', 'dimers.txt')
            self.assertTrue(filecmp.cmp(result_ql_dimers, expected_ql_contacting, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_dimers, expected_a6_contacting, shallow=False))

            result_ql_ids = os.path.join(outdir, 'ql', 'dimer_seq_ids.tsv')
            result_a6_ids = os.path.join(outdir, 'a6', 'dimer_seq_ids.tsv')
            self.assertTrue(filecmp.cmp(result_ql_ids, expected_ql_ids, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_ids, expected_a6_ids, shallow=False))

            result_all_homodimers = os.path.join(outinter, 'all_homodimers.txt')
            self.assertTrue(filecmp.cmp(result_all_homodimers, expected_all_homodimers, shallow=False))
            
            result_db_ql = os.path.join(outlib, 'notarealdb', 'ql.tsv')
            result_db_a6 = os.path.join(outlib, 'notarealdb', 'a6.tsv')
            self.assertTrue(filecmp.cmp(result_db_ql, expected_db_ql, shallow=False)) 
            self.assertTrue(filecmp.cmp(result_db_a6, expected_db_a6, shallow=False)) 


    def test_rule_all_fulldbs(self):
        rule = 'all'
        config = os.path.abspath('config_for_test_1_fulldbs.yaml')

        expected_ql_chains = os.path.join(data_dir, 'all', 'expected_ql_chains.txt')
        expected_a6_chains = os.path.join(data_dir, 'all', 'expected_a6_chains.txt')
        expected_ql_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_ql.txt')
        expected_a6_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_a6.txt')
        expected_ql_dimers = os.path.join(data_dir, 'all', 'expected_ql_contacts.txt')
        expected_a6_dimers = os.path.join(data_dir, 'all', 'expected_a6_contacts.txt') 
        expected_ql_ids = os.path.join(data_dir, 'all', 'expected_ql_ids.tsv')
        expected_a6_ids = os.path.join(data_dir, 'all', 'expected_a6_ids.tsv')
        expected_db_ql = os.path.join(data_dir, 'shared', 'lib', 'contactsdb_full', 'ql.tsv')
        expected_db_a6 = os.path.join(data_dir, 'shared', 'lib', 'contactsdb_full', 'a6.tsv')

        expected_all_homodimers = os.path.join(data_dir, 'all', 'expected_all_homodimers.txt') 

        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib) as tmpsnake:
            outinter = os.path.join(os.path.dirname(tmpsnake), 'intermediates')
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            outlib = os.path.join(os.path.dirname(tmpsnake), 'lib')

            result_ql_chains = os.path.join(outdir, 'ql', 'chains.txt')
            result_a6_chains = os.path.join(outdir, 'a6', 'chains.txt')
            self.assertTrue(filecmp.cmp(result_ql_chains, expected_ql_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_chains, expected_a6_chains, shallow=False))
           
            result_ql_poss = os.path.join(outdir, 'ql', 'intra_assembly_chain_pairs.txt')
            result_a6_poss = os.path.join(outdir, 'a6', 'intra_assembly_chain_pairs.txt')
            self.assertTrue(filecmp.cmp(result_ql_poss, expected_ql_poss_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_poss, expected_a6_poss_chains, shallow=False))

            # non empty lookup expected
            result_ql_lookup = os.path.join(outdir, 'ql', 'contacts_lookup.txt')
            result_a6_lookup = os.path.join(outdir, 'a6', 'contacts_lookup.txt')
            self.assertNotEqual(os.stat(result_ql_lookup).st_size, 0)
            self.assertNotEqual(os.stat(result_a6_lookup).st_size, 0)

            # empty calc expected
            result_ql_calc = os.path.join(outdir, 'ql', 'contacts_calc.txt')
            result_a6_calc = os.path.join(outdir, 'a6', 'contacts_calc.txt')
            self.assertEqual(os.stat(result_ql_calc).st_size, 0)
            self.assertEqual(os.stat(result_a6_calc).st_size, 0)

            # because empty calc, contacts need calc also empty
            result_ql_need = os.path.join(outdir, 'ql', 'contacts_needcalc.txt')
            result_a6_need = os.path.join(outdir, 'a6', 'contacts_needcalc.txt')
            self.assertEqual(os.stat(result_ql_need).st_size, 0)
            self.assertEqual(os.stat(result_a6_need).st_size, 0)

            result_ql_dimers = os.path.join(outdir, 'ql', 'dimers.txt')
            result_a6_dimers = os.path.join(outdir, 'a6', 'dimers.txt')
            self.assertTrue(filecmp.cmp(result_ql_dimers, expected_ql_dimers, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_dimers, expected_a6_dimers, shallow=False))

            result_ql_ids = os.path.join(outdir, 'ql', 'dimer_seq_ids.tsv')
            result_a6_ids = os.path.join(outdir, 'a6', 'dimer_seq_ids.tsv')
            self.assertTrue(filecmp.cmp(result_ql_ids, expected_ql_ids, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_ids, expected_a6_ids, shallow=False))

            result_all_homodimers = os.path.join(outinter, 'all_homodimers.txt')
            self.assertTrue(filecmp.cmp(result_all_homodimers, expected_all_homodimers, shallow=False))

            result_db_ql = os.path.join(outlib, 'contactsdb_full', 'ql.tsv')
            result_db_a6 = os.path.join(outlib, 'contactsdb_full', 'a6.tsv')
            self.assertTrue(filecmp.cmp(result_db_ql, expected_db_ql, shallow=False)) 
            self.assertTrue(filecmp.cmp(result_db_a6, expected_db_a6, shallow=False)) 


    def test_rule_all_partialdbs(self):
        rule = 'all'
        config = os.path.abspath('config_for_test_1_partialdbs.yaml')

        expected_ql_chains = os.path.join(data_dir, 'all', 'expected_ql_chains.txt')
        expected_a6_chains = os.path.join(data_dir, 'all', 'expected_a6_chains.txt')
        expected_ql_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_ql.txt')
        expected_a6_poss_chains = os.path.join(data_dir, 'all', 'expected_poss_chains_a6.txt')
        expected_ql_dimers = os.path.join(data_dir, 'all', 'expected_ql_contacts.txt')
        expected_a6_dimers = os.path.join(data_dir, 'all', 'expected_a6_contacts.txt') 
        expected_ql_ids = os.path.join(data_dir, 'all', 'expected_ql_ids.tsv')
        expected_a6_ids = os.path.join(data_dir, 'all', 'expected_a6_ids.tsv')
        expected_db_ql = os.path.join(data_dir, 'shared', 'lib', 'contactsdb_full', 'ql.tsv')
        expected_db_a6 = os.path.join(data_dir, 'shared', 'lib', 'contactsdb_full', 'a6.tsv')

        expected_all_homodimers = os.path.join(data_dir, 'all', 'expected_all_homodimers.txt') 

        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib) as tmpsnake:
            outinter = os.path.join(os.path.dirname(tmpsnake), 'intermediates')
            outdir = os.path.join(os.path.dirname(tmpsnake), outname)
            outlib = os.path.join(os.path.dirname(tmpsnake), 'lib')

            result_ql_chains = os.path.join(outdir, 'ql', 'chains.txt')
            result_a6_chains = os.path.join(outdir, 'a6', 'chains.txt')
            self.assertTrue(filecmp.cmp(result_ql_chains, expected_ql_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_chains, expected_a6_chains, shallow=False))
           
            result_ql_poss = os.path.join(outdir, 'ql', 'intra_assembly_chain_pairs.txt')
            result_a6_poss = os.path.join(outdir, 'a6', 'intra_assembly_chain_pairs.txt')
            self.assertTrue(filecmp.cmp(result_ql_poss, expected_ql_poss_chains, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_poss, expected_a6_poss_chains, shallow=False))

            # non empty lookup expected
            result_ql_lookup = os.path.join(outdir, 'ql', 'contacts_lookup.txt')
            result_a6_lookup = os.path.join(outdir, 'a6', 'contacts_lookup.txt')
            self.assertNotEqual(os.stat(result_ql_lookup).st_size, 0)
            self.assertEqual(os.stat(result_a6_lookup).st_size, 0) # but this was empty

            # non empty calc expected
            result_ql_calc = os.path.join(outdir, 'ql', 'contacts_calc.txt')
            result_a6_calc = os.path.join(outdir, 'a6', 'contacts_calc.txt')
            self.assertNotEqual(os.stat(result_ql_calc).st_size, 0)
            self.assertNotEqual(os.stat(result_a6_calc).st_size, 0)

            # because non empty calc, contacts need calc also non empty
            result_ql_need = os.path.join(outdir, 'ql', 'contacts_needcalc.txt')
            result_a6_need = os.path.join(outdir, 'a6', 'contacts_needcalc.txt')
            self.assertNotEqual(os.stat(result_ql_need).st_size, 0)
            self.assertNotEqual(os.stat(result_a6_need).st_size, 0) 

            result_ql_dimers = os.path.join(outdir, 'ql', 'dimers.txt')
            result_a6_dimers = os.path.join(outdir, 'a6', 'dimers.txt')
            self.assertTrue(filecmp.cmp(result_ql_dimers, expected_ql_dimers, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_dimers, expected_a6_dimers, shallow=False))

            result_ql_ids = os.path.join(outdir, 'ql', 'dimer_seq_ids.tsv')
            result_a6_ids = os.path.join(outdir, 'a6', 'dimer_seq_ids.tsv')
            self.assertTrue(filecmp.cmp(result_ql_ids, expected_ql_ids, shallow=False))
            self.assertTrue(filecmp.cmp(result_a6_ids, expected_a6_ids, shallow=False))

            result_all_homodimers = os.path.join(outinter, 'all_homodimers.txt')
            self.assertTrue(filecmp.cmp(result_all_homodimers, expected_all_homodimers, shallow=False))

            result_db_ql = os.path.join(outlib, 'contactsdb_partial', 'ql.tsv')
            result_db_a6 = os.path.join(outlib, 'contactsdb_partial', 'a6.tsv')
            self.assertTrue(filecmp.cmp(result_db_ql, expected_db_ql, shallow=False)) 
            self.assertTrue(filecmp.cmp(result_db_a6, expected_db_a6, shallow=False)) 


