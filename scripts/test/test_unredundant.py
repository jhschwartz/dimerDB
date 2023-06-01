import os
import sys
sys.path.append('..')

import unittest
from unredundant import *
import numpy as np

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
data_dir = f'{test_dir}/data/unredundant'
lib_path = os.path.join(data_dir, 'lib')
resolu_path = os.path.join(lib_path, 'resolu.idx')
methods_file = os.path.join(lib_path, 'pdb_entry_type.txt')


class TestUnredundantClusterFromMat(unittest.TestCase):

    def test_cluster_simple_fake_mat(self):
        mat = np.array([ [0., 0.9, 0.1, 0.9],
                         [0.9, 0., 0.9, 0.1],
                         [0.1, 0.9, 0., 0.9],
                         [0.9, 0.1, 0.9, 0.], ])
        result_labels = run_cluster_from_matrix(npmat=mat, thresh=0.5)
        self.assertEqual(result_labels[0], result_labels[2])
        self.assertEqual(result_labels[1], result_labels[3])
        self.assertNotEqual(result_labels[0], result_labels[1])
        self.assertEqual(len(set(result_labels)), 2)


    def test_cluster_one_element_mat(self):
        mat = np.array([0])
        result_labels = run_cluster_from_matrix(npmat=mat, thresh=0.5)
        self.assertEqual(mat, [0])


    def test_cluster_thresh_high(self):
        mat = np.array([ [0., 0.9, 0.1, 0.9],
                         [0.9, 0., 0.9, 0.1],
                         [0.1, 0.9, 0., 0.9],
                         [0.9, 0.1, 0.9, 0.], ])
        result_labels = run_cluster_from_matrix(npmat=mat, thresh=0.99)
        self.assertEqual(len(set(result_labels)), 1)


    def test_cluister_thresh_low(self):
        mat = np.array([ [0., 0.9, 0.1, 0.9],
                         [0.9, 0., 0.9, 0.1],
                         [0.1, 0.9, 0., 0.9],
                         [0.9, 0.1, 0.9, 0.], ])
        result_labels = run_cluster_from_matrix(npmat=mat, thresh=0.01)
        self.assertEqual(len(set(result_labels)), 4)




class TestUnredundantGetCluster(unittest.TestCase):

    def test_get_single_item(self):
        num = 0
        labs = [0]
        all_d = ['A']
        result = get_cluster_dimers(cluster_num=num, cluster_labels=labs, all_dimers=all_d)
        result = list(result)
        self.assertEqual(result, ['A'])


    def test_get_multi_items_one_cluster(self):
        num = 0
        labs = [0, 0, 0]
        all_d = ['A', 'B', 'C']
        result = get_cluster_dimers(cluster_num=num, cluster_labels=labs, all_dimers=all_d)
        result = list(result)
        self.assertEqual(result, ['A', 'B', 'C'])


    def test_get_multi_items_multi_cluster(self):
        num = 0
        labs = [0, 1, 0, 1]
        all_d = ['A', 'B', 'C', 'D']
        result = get_cluster_dimers(cluster_num=num, cluster_labels=labs, all_dimers=all_d)
        result = list(result)
        self.assertEqual(result, ['A', 'C'])




#   choose_rep(dimers, all_dimers, dist_mat, lib_path, resolu_path)
class TestUnredundantChooseRep(unittest.TestCase):
    
    def setUp(self):

        self.dimer_small_1 = '1sak-a1-m1-cC_1sak-a1-m1-cD' # 42 res/chain
        self.dimer_small_2 = '6rdj-a1-m1-cA_6rdj-a1-m1-cJ' # 72 res/chain
        
        # 327 res/chain, 1.301A
        self.dimer_large_xray_best_res = '5vej-a2-m1-cA_5vej-a2-m2-cA' 
        # 378 res/chain, 2.8A
        self.dimer_large_xray_second_res_1 = '4irn-a1-m1-cA_4irn-a1-m1-cD' 
        # 300 res/chain, 2.8A
        self.dimer_large_xray_second_res_2 = '4jrv-a1-m1-cA_4jrv-a1-m2-cA' 
        
        # 392 res/chain, EM res (ignored) 3.1A
        self.dimer_large_EM_1 = '7phk-a1-m1-cA_7phk-a1-m1-cD' 
        # 443 res/chain, EM res (ignored) 2.13A
        self.dimer_large_EM_2 = '7jfo-a1-m1-cC_7jfo-a1-m1-cG' 
        
        # 314 res/chain (double 3jc8-a1-m1-cA3_3jc8-a1-m1-cA6)
        self.dimer_large_NMR_1 = '3jc8-a1-m1-cA3_3jc8-a1-m1-cA6'
        # 314 res/chain (double 3jc8-a1-m1-cA2_3jc8-a1-m1-cAy)
        self.dimer_large_NMR_2 = '3jc8-a1-m1-cA2_3jc8-a1-m1-cAy' 

        self.all_dimers = [ self.dimer_small_1, 
                            self.dimer_small_2,
                            self.dimer_large_xray_best_res,
                            self.dimer_large_xray_second_res_1,
                            self.dimer_large_xray_second_res_2,
                            self.dimer_large_EM_1,
                            self.dimer_large_EM_2,
                            self.dimer_large_NMR_1, 
                            self.dimer_large_NMR_2                ]

        self.dist_mat = np.random.rand(len(self.all_dimers), len(self.all_dimers))



    def set_mat(self, dimer1, dimer2, val):
        i = self.all_dimers.index(dimer1)
        j = self.all_dimers.index(dimer2)
        self.dist_mat[i, j] = val
        self.dist_mat[j, i] = val



    # pick only dimers with sufficient coverage
    def test_rep_first_criterion(self):

        group = [self.dimer_small_1, self.dimer_small_2, self.dimer_large_EM_1]

        self.dist_mat = np.ones_like(self.dist_mat)

        result_rep = choose_rep(  dimers=group,
                                  all_dimers=self.all_dimers,
                                  dist_mat=self.dist_mat, 
                                  lib_path=lib_path, 
                                  resolu_path=resolu_path,
                                  methods_file=methods_file     )

        expected_rep = self.dimer_large_EM_1 # only not-too-small choice

        self.assertEqual(result_rep, expected_rep)



    # pick dimers past criterion 1 that are xray only
    def test_rep_second_criterion(self):
        
        group = [self.dimer_small_1, self.dimer_large_xray_second_res_1, 
                            self.dimer_large_EM_2, self.dimer_large_NMR_1]
        
        self.dist_mat = np.ones_like(self.dist_mat)

        result_rep = choose_rep(  dimers=group,
                                  all_dimers=self.all_dimers,
                                  dist_mat=self.dist_mat, 
                                  lib_path=lib_path, 
                                  resolu_path=resolu_path,
                                  methods_file=methods_file     )

        expected_rep = self.dimer_large_xray_second_res_1 # only xray choice

        self.assertEqual(result_rep, expected_rep)
        # '7jfo-a1-m1-cC_7jfo-a1-m1-cG' != '4irn-a1-m1-cA_4irn-a1-m1-cD'


    # pick only dimers that are most like the leftover, with 
    # xray dimers making it past criterion 2
    def test_rep_third_criterion_xray_exists(self):
        
        group = [self.dimer_large_xray_best_res, self.dimer_large_xray_second_res_1,
                                            self.dimer_large_xray_second_res_2        ]

        self.set_mat(  self.dimer_large_xray_best_res,
                       self.dimer_large_xray_second_res_1,
                       0.9                                  )
        self.set_mat(  self.dimer_large_xray_best_res,
                       self.dimer_large_xray_second_res_2,
                       0.1                                  )
        self.set_mat(  self.dimer_large_xray_second_res_1,
                       self.dimer_large_xray_second_res_2,
                       0.1                                  )
        
        result_rep = choose_rep(  dimers=group,
                                  all_dimers=self.all_dimers,
                                  dist_mat=self.dist_mat, 
                                  lib_path=lib_path, 
                                  resolu_path=resolu_path,
                                  methods_file=methods_file     )

        expected_rep = self.dimer_large_xray_second_res_2 # artificially most like xray choices

        self.assertEqual(result_rep, expected_rep)



    # pick only dimers that are most like the leftover, without 
    # any xray dimers making it past criterion 2
    def test_rep_third_criterion_no_xray(self):
        
        group = [self.dimer_large_NMR_1, self.dimer_large_NMR_2, self.dimer_large_EM_1, 
                                                                self.dimer_large_EM_2   ]   
        
        self.set_mat(  self.dimer_large_NMR_1,
                       self.dimer_large_NMR_2,
                       0.1                                  )
        self.set_mat(  self.dimer_large_NMR_1,
                       self.dimer_large_EM_1,
                       0.9                                  )
        self.set_mat(  self.dimer_large_NMR_1,
                       self.dimer_large_EM_2,
                       0.9                                  )
        self.set_mat(  self.dimer_large_NMR_2,
                       self.dimer_large_EM_1,
                       0.1                                  )
        self.set_mat(  self.dimer_large_NMR_2,
                       self.dimer_large_EM_2,
                       0.1                                  )
        self.set_mat(  self.dimer_large_EM_1,
                       self.dimer_large_EM_2,
                       0.9                                  )


        result_rep = choose_rep(  dimers=group,
                                  all_dimers=self.all_dimers,
                                  dist_mat=self.dist_mat, 
                                  lib_path=lib_path, 
                                  resolu_path=resolu_path,
                                  methods_file=methods_file     )

        expected_rep = self.dimer_large_NMR_2 # artifically most like others

        self.assertEqual(result_rep, expected_rep)


    # pick only dimers past third that are best resolution (xray only)
    def test_rep_fourth_criterion(self):
        
        group = [ self.dimer_large_NMR_1, self.dimer_large_NMR_2, self.dimer_large_EM_1, 
                  self.dimer_large_EM_2, self.dimer_large_xray_best_res, 
                  self.dimer_large_xray_second_res_1, self.dimer_large_xray_second_res_2  ]
        
        self.dist_mat = np.ones_like(self.dist_mat)
       

        result_rep = choose_rep(  dimers=group,
                                  all_dimers=self.all_dimers,
                                  dist_mat=self.dist_mat, 
                                  lib_path=lib_path, 
                                  resolu_path=resolu_path,
                                  methods_file=methods_file     )

        expected_rep = self.dimer_large_xray_best_res # best xray choice resolution

        self.assertEqual(result_rep, expected_rep)


    # pick dimers past fourth by first alphabetically 
    def test_rep_fifth_criterion_xray(self):
       
        group = [ self.dimer_large_NMR_1, self.dimer_large_NMR_2, self.dimer_large_EM_1, 
                            self.dimer_large_EM_2,
                  self.dimer_large_xray_second_res_1, self.dimer_large_xray_second_res_2  ]
        
        self.dist_mat = np.ones_like(self.dist_mat)

        
        result_rep = choose_rep(  dimers=group,
                                  all_dimers=self.all_dimers,
                                  dist_mat=self.dist_mat, 
                                  lib_path=lib_path, 
                                  resolu_path=resolu_path,
                                  methods_file=methods_file     )

        expected_rep = self.dimer_large_xray_second_res_2 # last alpha vs other xray choice

        self.assertEqual(result_rep, expected_rep)


    # pick dimers that had no fourth criterion because not 
    # xray, pick by alphabetically first
    def test_rep_fifth_criterion_no_xray(self):
       
        group = [ self.dimer_large_NMR_1, self.dimer_large_NMR_2, self.dimer_large_EM_1, 
                            self.dimer_large_EM_2                                        ]
        
        self.dist_mat = np.ones_like(self.dist_mat)

        result_rep = choose_rep(  dimers=group,
                                  all_dimers=self.all_dimers,
                                  dist_mat=self.dist_mat, 
                                  lib_path=lib_path, 
                                  resolu_path=resolu_path,
                                  methods_file=methods_file     )

        expected_rep = self.dimer_large_EM_1 # last alpha of group 

        self.assertEqual(result_rep, expected_rep)





