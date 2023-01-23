# '''
# prune_dimers.py - functions to remove structural redundancy from lists of dimers.

# Written by Jacob Schwartz (jaschwa@umich.edu) in November, 2022.
# Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
# https://freddolino-lab.med.umich.edu

# These functions are untested as of 11/18/2022.
# This work requires python >= 3.8
# '''
# import subprocess
# import random
# import re
# import os
# #import glob
# #import requests
# from prune_redundancy import prune_redundancy

# # todo : make Dimer compatible with generatic redundancy pruner

# class Dimer:
#     def __init__(self, chain1, chain2):
#         self.chain1 = chain1
#         self.chain2 = chain2

#     def download(self):
#         # TODO
#         raise NotImplementedError
#         self.chain1_file = None
#         self.chain2_file = None
#         self.complex_file = None

#     def clear(self):
#         os.remove(self.chain1_file)
#         os.remove(self.chain2_file)
#         os.remove(self.complex_file)
#         self.chain1_file = None
#         self.chain2_file = None
#         self.complex_file = None

#     def coverage(self):
#         # TODO
#         raise NotImplementedError




# #def dimer2data(dimer):
# #    '''
# #    Retrives the pdb files of a dimer.
# #
# #    :param dimer: tuple(str, str), the strings defining the dimer pair
# #    :returns: tuple(str, str), the data, each string is an entire
# #              pdb file.
# #    '''
# #    raise NotImplementedError



# def MMalign_distance(dimer1, dimer2, MMalign='../bin/MMalign'):
#     '''
#     Compute the distance between two dimers, which is equal to 
#     1 - (MMalign score)

#     :param dimer1: Dimer, the first thing to compare
#     :param dimer2: Dimer, the second thing to compare
#     :returns: float, the distance in structural similarity
#     '''
# #    rand = str(random.randint(1,99999999999999))
# #    complex1_tmp = f'/tmp/dimer1assembly-{rand}.pdb.tmp'
# #    complex2_tmp = f'/tmp/dimer2assembly-{rand}.pdb.tmp'
# #
# #    with open(complex1_tmp, 'w') as f1, open(complex2_tmp, 'w') as f2:
# #        f1.write(dimer1_A)
# #        f1.write('\n')
# #        f1.write(dimer1_B)
# #        f2.write(dimer2_A)
# #        f2.write('\n')
# #        f2.write(dimer2_B)

#     cmd = f'{MMalign} {dimer1.complex_file} {dimer2.complex_file}'
#     result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
#     if result.stderr != '':
#         raise RuntimeError(result.stderr)

#     TMscores = re.findall('TM-score=\s([\d\.]+)', result.stdout)
#     if len(TMscores) != 2:
#         raise ValueError(f'was unable to retrieve both TMscores, instead found {len(TMscores)} scores.')
#     score = max([float(tm) for tm in TMscores])

#     return 1 - score

        


# def dimer_cluster(dimers, distance_matrix, threshold):
#     '''
#     Wrapper for hierarchical agglomerative clustering using sklearn.

#     :param dimers: list[Dimer], our dimers to cluster
#     :param distance_matrix: np.arr[float, float], the pairwise 
#                             distance between the dimers
#     :param threshold: float, the maximum distance by which we consider 
#                       dimers similar (and of the same cluster)
#     :returns: list[list[Dimer]], a list of clusters in which 
#               each member list is one cluster list of Dimers 
#     '''
#     raise NotImplementedError



# def dimer_representative(cluster):
#     '''
#     Picks the representaive dimer of a cluster of dimers. Does so by
#     picking the dimer with the most sequence coverage. In a tiebreaker,
#     picks the dimer with the newest publication date.

#     :param cluster: list[Dimer], a list of dimers all in one cluster 
#     :returns dimer: Dimer, the dimer most representative of the cluster
#     '''
#     max_cover = -1
#     rep = None
#     for dimer in cluster:
#         cover = dimer.coverage()
#         if cover > max_cover:
#             max_cover = cover
#             rep = dimer
#     return rep




# def prune_dimers(dimers, threshold=0.5):
#     '''
#     Removes redundancy from a list of dimers.
#         - uses MMalign to compute distance between dimers.
#         - clusters using dimer_cluster, an wrapper for scikit learn 
#           hierarchical agglomerative clustering
#         - picks cluster represenatives by dimer_representative, 
#           which does so by primarily by most sequence coverage by
#           the model and second by most recent publication
#         - retrieves data by get_dimer, which checks a data folder
#           for the structure and downloads it if it is not there.
#           This data folder is meant to be temporary for only one
#           run of the prune_dimers function.
#     '''
#     for d in dimers:
#         d.download()
#     non_redundant =  prune_redundancy(things=dimers, 
#                             distance=MMalign_distance, 
#                             cluster=dimer_cluster, 
#                             threshold=threshold, 
#                             representative=dimer_representative, 
#                             thing2data=get_dimer)
#     for d in dimers:
#         d.clear()
#     return non_redundant
