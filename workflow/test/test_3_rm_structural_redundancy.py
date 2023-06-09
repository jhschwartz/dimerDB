import unittest
import os
import filecmp

from run_tmp_snakemake import run_tmp_snakemake

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = os.path.join(test_dir, 'data', '3_rm_structural_redundancy')

conda_env = os.environ['conda_env']

snakefile = os.path.join('..', '3_rm_structural_redundancy.smk')
snake_exe = os.path.join(conda_env, 'bin/snakemake')
lib = os.path.abspath(os.path.join(data_dir, 'lib'))
scripts = os.path.abspath('../../scripts')
bin = os.path.abspath('../../bin')



seq_clusters = ['1bb1-a1-m1-cD', '2cc2-a1-m1-cB', '3aa3-a1-m1-cD']

aa_expected = os.path.join(lib, 'tmdb_expected', 'aa.tsv')
bb_expected = os.path.join(lib, 'tmdb_expected', 'bb.tsv')
cc_expected = os.path.join(lib, 'tmdb_expected', 'cc.tsv')



class TestSubworkflow3(unittest.TestCase):

    def test_sub3_calc_all_tms(self):
        # in this case, no start tmdb is provided at all, forcing calc of all.
        config = os.path.abspath('config_for_test_3_tmdb_notexist.yaml')
       
        rule = 'all'
        intermediates = os.path.join(data_dir, 'intermediates')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outinter = os.path.join(os.path.dirname(tmpsnake), 'intermediates')
            outlib = os.path.join(os.path.dirname(tmpsnake), 'lib')
            outclusts = [os.path.join(outinter, 'cluster', 'seq_clusters', c) for c in seq_clusters]

            # check that lookup file is empty:
            for result_lookup in [os.path.join(c, 'dist_lookup.tsv') for c in outclusts]:
                self.assertEqual(os.stat(result_lookup).st_size, 0)

            # check that calc results are correct:
            with open(os.path.join(outclusts[0], 'dist_calc.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 91) # 14 choose 2
                expected = '1cc1-a1-m1-cC_1cc1-a1-m1-cD\t1cc1-a1-m1-cD_1cc1-a1-m1-cE\t0.96\t0.97'
                self.assertTrue(expected in lines)
            with open(os.path.join(outclusts[1], 'dist_calc.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 300) # 25 choose 2 
                expected = '2cc2-a1-m1-cB_2cc2-a1-m1-cC\t2cc2-a1-m1-cD_2cc2-a1-m1-cE\t0.95\t0.97'
                self.assertTrue(expected in lines)
            with open(os.path.join(outclusts[2], 'dist_calc.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 190) # 20 choose 2
                expected = '3cc3-a1-m1-cD_3cc3-a1-m1-cH\t3cc3-a1-m1-cG_3cc3-a1-m1-cH\t0.99\t0.99'
                self.assertTrue(expected in lines)
            
            # check resolu
            self.assertTrue(os.path.exists(os.path.join(outinter, 'resolu.idx')))

            # check methods
            self.assertTrue(os.path.exists(os.path.join(outinter, 'pdb_entry_type.txt')))

            # check cluster results
            with open(os.path.join(outclusts[0], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '1aa1-a1-m1-cB_1aa1-a1-m1-cF')
            with open(os.path.join(outclusts[0], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 14)
                self.assertTrue('1cc1-a1-m1-cB_1cc1-a1-m1-cF\t1aa1-a1-m1-cB_1aa1-a1-m1-cF'\
                                    in lines)
                self.assertTrue('1bb1-a1-m1-cD_1bb1-a1-m1-cE\t1aa1-a1-m1-cB_1aa1-a1-m1-cF'\
                                    in lines)
            with open(os.path.join(outclusts[1], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '2bb2-a1-m1-cB_2bb2-a1-m1-cC')
            with open(os.path.join(outclusts[1], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 25)
                self.assertTrue('2aa2-a2-m1-cE_2aa2-a2-m1-cG\t2bb2-a1-m1-cB_2bb2-a1-m1-cC'\
                                    in lines)
                self.assertTrue('2cc2-a1-m1-cC_2cc2-a1-m1-cD\t2bb2-a1-m1-cB_2bb2-a1-m1-cC'\
                                    in lines)
            with open(os.path.join(outclusts[2], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '3aa3-a1-m1-cD_3aa3-a1-m1-cH')
            with open(os.path.join(outclusts[2], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 20)
                self.assertTrue('3cc3-a1-m1-cD_3cc3-a1-m1-cE\t3aa3-a1-m1-cD_3aa3-a1-m1-cH'\
                                    in lines)
                self.assertTrue('3aa3-a1-m1-cG_3aa3-a1-m1-cH\t3aa3-a1-m1-cD_3aa3-a1-m1-cH'\
                                    in lines)

            result_tmdb = os.path.join(outlib, 'tmdb_notexist')
            aa = os.path.join(result_tmdb, 'aa.tsv')
            bb = os.path.join(result_tmdb, 'bb.tsv')
            cc = os.path.join(result_tmdb, 'cc.tsv')
            self.assertTrue(filecmp.cmp(aa, aa_expected))
            self.assertTrue(filecmp.cmp(bb, bb_expected))
            self.assertTrue(filecmp.cmp(cc, cc_expected))


    def test_sub3_lookup_all_tms(self):
        # in this case, the full tmdb is provided, forcing an entire lookup and no calc.
        config = os.path.abspath('config_for_test_3_tmdb_all.yaml')
       
        rule = 'all'
        intermediates = os.path.join(data_dir, 'intermediates')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outinter = os.path.join(os.path.dirname(tmpsnake), 'intermediates')
            outlib = os.path.join(os.path.dirname(tmpsnake), 'lib')
            outclusts = [os.path.join(outinter, 'cluster', 'seq_clusters', c) for c in seq_clusters]


            # check that lookup results are correct:
            with open(os.path.join(outclusts[0], 'dist_lookup.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 91) # 14 choose 2
                expected = '1cc1-a1-m1-cC_1cc1-a1-m1-cD\t1cc1-a1-m1-cD_1cc1-a1-m1-cE\t0.96\t0.97'
                self.assertTrue(expected in lines)
            with open(os.path.join(outclusts[1], 'dist_lookup.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 300) # 25 choose 2 
                expected = '2cc2-a1-m1-cB_2cc2-a1-m1-cC\t2cc2-a1-m1-cD_2cc2-a1-m1-cE\t0.95\t0.97'
                self.assertTrue(expected in lines)
            with open(os.path.join(outclusts[2], 'dist_lookup.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 190) # 20 choose 2
                expected = '3cc3-a1-m1-cD_3cc3-a1-m1-cH\t3cc3-a1-m1-cG_3cc3-a1-m1-cH\t0.99\t0.99'
                self.assertTrue(expected in lines)
            
            # check that calc files empty
            for result_calc in [os.path.join(c, 'dist_calc.tsv') for c in outclusts]:
                self.assertEqual(os.stat(result_calc).st_size, 0)
            
            # check resolu
            self.assertTrue(os.path.exists(os.path.join(outinter, 'resolu.idx')))

            # check methods
            self.assertTrue(os.path.exists(os.path.join(outinter, 'pdb_entry_type.txt')))

            # check cluster results
            with open(os.path.join(outclusts[0], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '1aa1-a1-m1-cB_1aa1-a1-m1-cF')
            with open(os.path.join(outclusts[0], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 14)
                self.assertTrue('1cc1-a1-m1-cB_1cc1-a1-m1-cF\t1aa1-a1-m1-cB_1aa1-a1-m1-cF'\
                                    in lines)
                self.assertTrue('1bb1-a1-m1-cD_1bb1-a1-m1-cE\t1aa1-a1-m1-cB_1aa1-a1-m1-cF'\
                                    in lines)
            with open(os.path.join(outclusts[1], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '2bb2-a1-m1-cB_2bb2-a1-m1-cC')
            with open(os.path.join(outclusts[1], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 25)
                self.assertTrue('2aa2-a2-m1-cE_2aa2-a2-m1-cG\t2bb2-a1-m1-cB_2bb2-a1-m1-cC'\
                                    in lines)
                self.assertTrue('2cc2-a1-m1-cC_2cc2-a1-m1-cD\t2bb2-a1-m1-cB_2bb2-a1-m1-cC'\
                                    in lines)
            with open(os.path.join(outclusts[2], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '3aa3-a1-m1-cD_3aa3-a1-m1-cH')
            with open(os.path.join(outclusts[2], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 20)
                self.assertTrue('3cc3-a1-m1-cD_3cc3-a1-m1-cE\t3aa3-a1-m1-cD_3aa3-a1-m1-cH'\
                                    in lines)
                self.assertTrue('3aa3-a1-m1-cG_3aa3-a1-m1-cH\t3aa3-a1-m1-cD_3aa3-a1-m1-cH'\
                                    in lines)

            result_tmdb = os.path.join(outlib, 'tmdb_all')
            aa = os.path.join(result_tmdb, 'aa.tsv')
            bb = os.path.join(result_tmdb, 'bb.tsv')
            cc = os.path.join(result_tmdb, 'cc.tsv')
            self.assertTrue(filecmp.cmp(aa, aa_expected))
            self.assertTrue(filecmp.cmp(bb, bb_expected))
            self.assertTrue(filecmp.cmp(cc, cc_expected))






    def test_sub3_lookup_some_calc_some_tms(self):
        # in this case, a half tmdb is provided, forcing some lookup, some calc
        config = os.path.abspath('config_for_test_3_tmdb_half.yaml')
       
        rule = 'all'
        intermediates = os.path.join(data_dir, 'intermediates')
        with run_tmp_snakemake(snakefile, config, snake_exe, scripts, bin, rule, lib, intermediates) as tmpsnake:
            outinter = os.path.join(os.path.dirname(tmpsnake), 'intermediates')
            outlib = os.path.join(os.path.dirname(tmpsnake), 'lib')
            outclusts = [os.path.join(outinter, 'cluster', 'seq_clusters', c) for c in seq_clusters]

            # check that lookup files not empty
            for result_lookup in [os.path.join(c, 'dist_lookup.tsv') for c in outclusts]:
                self.assertNotEqual(os.stat(result_lookup).st_size, 0)
            
            # check that calc files not empty
            for result_calc in [os.path.join(c, 'dist_calc.tsv') for c in outclusts]:
                self.assertNotEqual(os.stat(result_calc).st_size, 0)

            # check aa items were lookup
            with open(os.path.join(outclusts[0], 'dist_lookup.tsv'), 'r') as f:
                lines = f.read().splitlines()
                expected = '1aa1-a1-m1-cD_1aa1-a1-m1-cE\t1cc1-a1-m1-cD_1cc1-a1-m1-cE\t0.99\t1.0'
                self.assertTrue(expected in lines)
            with open(os.path.join(outclusts[1], 'dist_lookup.tsv'), 'r') as f:
                lines = f.read().splitlines()
                expected = '2aa2-a2-m1-cG_2aa2-a2-m1-cH\t2bb2-a2-m1-cF_2bb2-a2-m1-cJ\t0.96\t0.98'
                self.assertTrue(expected in lines)
            with open(os.path.join(outclusts[2], 'dist_lookup.tsv'), 'r') as f:
                lines = f.read().splitlines()
                expected = '3aa3-a1-m1-cD_3aa3-a1-m1-cE\t3aa3-a1-m1-cE_3aa3-a1-m1-cF\t0.99\t0.99'
                self.assertTrue(expected in lines)
            

            # check cc items were calc
            with open(os.path.join(outclusts[0], 'dist_calc.tsv'), 'r') as f:
                lines = f.read().splitlines()
                expected = '1cc1-a1-m1-cB_1cc1-a1-m1-cF\t1cc1-a1-m1-cC_1cc1-a1-m1-cD\t0.99\t0.98'
                self.assertTrue(expected in lines)
            with open(os.path.join(outclusts[1], 'dist_calc.tsv'), 'r') as f:
                lines = f.read().splitlines()
                expected = '2cc2-a1-m1-cA_2cc2-a1-m1-cB\t2cc2-a1-m1-cC_2cc2-a1-m1-cD\t0.95\t0.97'
                self.assertTrue(expected in lines)
            with open(os.path.join(outclusts[2], 'dist_calc.tsv'), 'r') as f:
                lines = f.read().splitlines()
                expected = '3cc3-a1-m1-cD_3cc3-a1-m1-cH\t3cc3-a1-m1-cF_3cc3-a1-m1-cG\t0.99\t0.99'
                self.assertTrue(expected in lines)


            
            # check resolu
            self.assertTrue(os.path.exists(os.path.join(outinter, 'resolu.idx')))

            # check methods
            self.assertTrue(os.path.exists(os.path.join(outinter, 'pdb_entry_type.txt')))

            # check cluster results
            with open(os.path.join(outclusts[0], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '1aa1-a1-m1-cB_1aa1-a1-m1-cF')
            with open(os.path.join(outclusts[0], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 14)
                self.assertTrue('1cc1-a1-m1-cB_1cc1-a1-m1-cF\t1aa1-a1-m1-cB_1aa1-a1-m1-cF'\
                                    in lines)
                self.assertTrue('1bb1-a1-m1-cD_1bb1-a1-m1-cE\t1aa1-a1-m1-cB_1aa1-a1-m1-cF'\
                                    in lines)
            with open(os.path.join(outclusts[1], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '2bb2-a1-m1-cB_2bb2-a1-m1-cC')
            with open(os.path.join(outclusts[1], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 25)
                self.assertTrue('2aa2-a2-m1-cE_2aa2-a2-m1-cG\t2bb2-a1-m1-cB_2bb2-a1-m1-cC'\
                                    in lines)
                self.assertTrue('2cc2-a1-m1-cC_2cc2-a1-m1-cD\t2bb2-a1-m1-cB_2bb2-a1-m1-cC'\
                                    in lines)
            with open(os.path.join(outclusts[2], 'representatives.txt'), 'r') as f:
                line = f.read().strip()
                self.assertEqual(line, '3aa3-a1-m1-cD_3aa3-a1-m1-cH')
            with open(os.path.join(outclusts[2], 'representation.tsv'), 'r') as f:
                lines = f.read().splitlines()
                self.assertEqual(len(lines), 20)
                self.assertTrue('3cc3-a1-m1-cD_3cc3-a1-m1-cE\t3aa3-a1-m1-cD_3aa3-a1-m1-cH'\
                                    in lines)
                self.assertTrue('3aa3-a1-m1-cG_3aa3-a1-m1-cH\t3aa3-a1-m1-cD_3aa3-a1-m1-cH'\
                                    in lines)

            result_tmdb = os.path.join(outlib, 'tmdb_half')
            aa = os.path.join(result_tmdb, 'aa.tsv')
            bb = os.path.join(result_tmdb, 'bb.tsv')
            cc = os.path.join(result_tmdb, 'cc.tsv')
            self.assertTrue(filecmp.cmp(aa, aa_expected))
            self.assertTrue(filecmp.cmp(bb, bb_expected))
            self.assertTrue(filecmp.cmp(cc, cc_expected))

