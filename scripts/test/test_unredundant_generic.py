import unittest
import sys
import numpy as np
import tempfile

sys.path.append('..')
from unredundant import RedundantThings


class RedundantGeneric(RedundantThings):
    def __init__(self, numbers, threshold):
        super().__init__(things=numbers, threshold=threshold)


    @classmethod 
    def distance(cls, number1, number2):
        diff = number1 - number2
        if diff < 0:
            return -1 * diff
        return diff

    
    @staticmethod
    def representative(cluster):
        return sorted(cluster)[int(len(cluster)/2)]



class TestRedundantGeneric(unittest.TestCase):
    def test_init(self):
        rg = RedundantGeneric(numbers=[1,2,3], threshold=0.1)
        self.assertEqual([1,2,3], rg.things)
        self.assertEqual(0.1, rg.threshold)
        

    def test_distance(self):
        dist = RedundantGeneric.distance(5, 6)
        self.assertEqual(dist, 1)

    def test_cluster(self):
        n1 = list(range(5))
        n2 = list(range(500,505))
        numbers = n1 + n2
        rg = RedundantGeneric(numbers=numbers, threshold=10)
        rg.initiate_distance_matrix()
        num_clusters = rg.initiate_clusters()
        self.assertEqual(num_clusters, 2)
        cluster1_index = rg.cluster_index_of_thing(0)
        cluster2_index = rg.cluster_index_of_thing(500)
        cluster1 = rg.retrieve_cluster(cluster1_index)
        cluster2 = rg.retrieve_cluster(cluster2_index)
        self.assertEqual(cluster1, list(range(5)))
        self.assertEqual(cluster2, list(range(500,505)))


    def test_distance_matrix_mono(self):
        rg = RedundantGeneric(numbers=[1,2,3], threshold=0.1)
        rg.initiate_distance_matrix()
        self.assertTrue(np.array_equal(rg.distance_matrix[0,:], np.array([0,1,2])))
        self.assertTrue(np.array_equal(rg.distance_matrix[1,:], np.array([1,0,1])))
        self.assertTrue(np.array_equal(rg.distance_matrix[2,:], np.array([2,1,0])))


    def test_distance_matrix_multi(self):
        rg = RedundantGeneric(numbers=[1,2,3], threshold=0.1)
        rg.initiate_distance_matrix(num_workers=2)
        self.assertTrue(np.array_equal(rg.distance_matrix[0,:], np.array([0,1,2])))
        self.assertTrue(np.array_equal(rg.distance_matrix[1,:], np.array([1,0,1])))
        self.assertTrue(np.array_equal(rg.distance_matrix[2,:], np.array([2,1,0])))


    def test_things_of_lowest_distance_to_others(self):
        rg = RedundantGeneric(numbers=[], threshold=0.1)
        result = rg._things_of_lowest_distance_to_others([1,2,3])
        self.assertEqual(result, [2])
        result = rg._things_of_lowest_distance_to_others([1,2,3,4])
        self.assertEqual(result, [2,3])
    
    
    def test_representative(self):
        result = RedundantGeneric.representative([3,4,1])
        self.assertEqual(result, 3)
        result = RedundantGeneric.representative([3,4,1,2])
        self.assertEqual(result, 3)
        

    def test_prune(self):
        n1 = list(range(5))
        n2 = list(range(500,505))
        numbers = n1 + n2
        rg = RedundantGeneric(numbers=numbers, threshold=10)
        result = rg.prune_redundancy()
        self.assertEqual(sorted(result), [2, 502])

    
    def test_cluster_index_of_thing(self):
        n1 = list(range(5))
        n2 = list(range(500,505))
        numbers = n1 + n2
        rg = RedundantGeneric(numbers=numbers, threshold=10)
        rg.initiate_distance_matrix()
        rg.initiate_clusters()
        i1 = rg.cluster_index_of_thing(4)
        i2 = rg.cluster_index_of_thing(3)
        i3 = rg.cluster_index_of_thing(501)
        i4 = rg.cluster_index_of_thing(502)
        self.assertEqual(i1, i2)
        self.assertEqual(i3, i4)
        self.assertNotEqual(i1, i3)


    def test_save_load_distance_matrix(self):
        n1 = list(range(5))
        n2 = list(range(500,505))
        numbers = n1 + n2
        rg = RedundantGeneric(numbers=numbers, threshold=10)
        rg.initiate_distance_matrix()
        old_mat = rg.distance_matrix.copy()
        with tempfile.TemporaryDirectory() as td:
            matfile = f'{td}/mat.npy'
            rg.save_dist_matrix(savefile=matfile)
            rg.distance_matrix = None
            rg.load_dist_matrix(loadfile=matfile)
        self.assertTrue(np.array_equal(rg.distance_matrix, old_mat))


    def test_prune_from_saved_distance_matrix(self):
        n1 = list(range(5))
        n2 = list(range(500,505))
        numbers = n1 + n2
        rg = RedundantGeneric(numbers=numbers, threshold=10)
        rg.initiate_distance_matrix()
        with tempfile.TemporaryDirectory() as td:
            matfile = f'{td}/mat.npy'
            rg.save_dist_matrix(savefile=matfile)
            rg.distance_matrix = None
            rg.load_dist_matrix(loadfile=matfile)
        result = rg.prune_redundancy(calc_dist_matrix=False)
        self.assertEqual(sorted(result), [2, 502])


    def test_noprune_without_distance_matrix(self):
        n1 = list(range(5))
        n2 = list(range(500,505))
        numbers = n1 + n2
        rg = RedundantGeneric(numbers=numbers, threshold=10)
        with self.assertRaises(ValueError):
            rg.prune_redundancy(calc_dist_matrix=False)


    def test_initiate_clusters(self):
        n1 = list(range(5))
        n2 = list(range(500,505))
        numbers = n1 + n2
        rg = RedundantGeneric(numbers=numbers, threshold=10)
        rg.initiate_distance_matrix()
        num_clusters = rg.initiate_clusters()
        self.assertEqual(num_clusters, 2)
        left_label = rg.things_cluster_labels[0]
        right_label = rg.things_cluster_labels[-1]
        l, r = left_label, right_label
        expected = [l,l,l,l,l,r,r,r,r,r]
        self.assertTrue(all([e==l for e,l in zip(expected, rg.things_cluster_labels)]))
         

    def test_initiate_cluster_skip_one_dimer(self):
        numbers = [100]
        rg = RedundantGeneric(numbers=numbers, threshold=10)
        rg.initiate_distance_matrix()
        num_clusters = rg.initiate_clusters()
        self.assertEqual(num_clusters, 1)
        self.assertEqual([0], rg.things_cluster_labels)
        self.assertTrue(all([e==l for e,l in zip([0], rg.things_cluster_labels)]))


if __name__ == '__main__':
    unittest.main()
