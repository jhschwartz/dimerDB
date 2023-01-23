'''
check_chains_contact.py - functions to check that two chains of one PDB assembly, each in their own file, make contact

Written by Jacob Schwartz (jaschwa@umich.edu) in December 2022.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

These functions are unitested by test/test_check_chains_contact.py and passing as of 1/12/2023.
This work requires python >= 3.8
'''

import numpy as np


def _check_chains_contact_pairwise(chain1_pdb, chain2_pdb, label_funcs, contact_distance_threshold, contact_count_threshold):
    res_dic_1, nums_1, _  = label_funcs.read_chain(chain1_pdb)
    res_dic_2, nums_2, _  = label_funcs.read_chain(chain2_pdb)
    _, _, labels, _, _ = label_funcs.extract_labels_dimer(res_dic_1, nums_1, res_dic_2, nums_2, True)
    dist = labels['ca_dis']
    num_contacts = np.count_nonzero(dist < contact_distance_threshold)
    return num_contacts >= contact_count_threshold



def _dim_impossible(max_x_1, min_x_1, max_x_2, min_x_2, threshold):
    return min_x_2 - max_x_1 > 8 or min_x_1 - max_x_2 > 8
        


def _check_chains_contact_impossible(chain1_pdb, chain2_pdb, label_funcs, contact_distance_threshold):
    res_dic_1, _, _  = label_funcs.read_chain(chain1_pdb)
    res_dic_2, _, _  = label_funcs.read_chain(chain2_pdb)

    chain1_ca_locs = np.array([coordinate for residue, coordinate in res_dic_1.items() if 'CA' in residue])
    chain2_ca_locs = np.array([coordinate for residue, coordinate in res_dic_2.items() if 'CA' in residue])

    chain1_x_max, chain1_y_max, chain1_z_max = np.max(chain1_ca_locs.squeeze(), axis=0)
    chain1_x_min, chain1_y_min, chain1_z_min = np.min(chain1_ca_locs.squeeze(), axis=0)
    chain2_x_max, chain2_y_max, chain2_z_max = np.max(chain2_ca_locs.squeeze(), axis=0)
    chain2_x_min, chain2_y_min, chain2_z_min = np.min(chain2_ca_locs.squeeze(), axis=0)

    if _dim_impossible(chain1_x_max, chain1_x_min, chain2_x_max, chain2_x_min, contact_distance_threshold):
        return True

    if _dim_impossible(chain1_y_max, chain1_y_min, chain2_y_max, chain2_y_min, contact_distance_threshold):
        return True
    
    if _dim_impossible(chain1_z_max, chain1_z_min, chain2_z_max, chain2_z_min, contact_distance_threshold):
        return True

    return False



def check_chains_contact(chain1_pdb, chain2_pdb, label_funcs, contact_distance_threshold, contact_count_threshold):
    if _check_chains_contact_impossible(chain1_pdb, chain2_pdb, label_funcs, contact_distance_threshold):
        return False
    return _check_chains_contact_pairwise(chain1_pdb, chain2_pdb, label_funcs, contact_distance_threshold, contact_count_threshold)


