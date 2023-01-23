'''
name_pdb.py - functions to give the path to pdb files of chains in the rcsb pdb, assuming that the pdb
              has been downloaded according to the organization of the dimerDB pipeline.
              Works whether a pdb chain belongs to the normal PDB or to an oversized structure that
              was split and renamed. Assumes a lib_path, which should be defined by the config file
              and passed to here.

Written by Jacob Schwartz (jaschwa@umich.edu) in January 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

This function is unittested by test/test_name_pdb.py and passing as of 1/12/2023.
This work requires python >= 3.8
'''
import os


def name_pdb_file(pdb_base, chain, lib_path):
    div = pdb_base[1:3]
    path = f'{lib_path}/rcsb_pdb/{div}/{pdb_base}{chain}.pdb'
    if not os.path.exists(path):
        # try oversized
        path = f'{lib_path}/rcsb_oversized/{div}/{pdb_base}/split_renamed/{pdb_base}{chain}.pdb'
        if not os.path.exists(path):
            raise FileNotFoundError(f'was unable to find pdb file {path} for {pdb_base}/{chain}.') 
    return path



def dimer2pdbs(dimer_name, lib_path):
    pdb0, chain0 = dimer_name.split('-')[0].split('_')
    pdb1, chain1 = dimer_name.split('-')[1].split('_')
    return name_pdb_file(pdb0, chain0, lib_path), name_pdb_file(pdb1, chain1, lib_path)
    

