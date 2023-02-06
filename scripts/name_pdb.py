'''
name_pdb.py - functions to give the path to pdb files of chains in the rcsb pdb, assuming that the pdb
              has been downloaded according to the organization of the dimerDB pipeline.
              Assumes the chains in question have been appropriately converted from their mmCIF assembly
              source and converted into individual chain PDB files for each model of the assembly. 
              Additionnaly assumes a lib_path, which should be defined by the config file
              and passed to here.

Written by Jacob Schwartz (jaschwa@umich.edu) in January 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

This function is unittested by test/test_name_pdb.py and passing as of 2/1/2023.
This work requires python >= 3.8
'''
import os
import glob
import re


# individual pdb files take the naming scheme "1abc-X-1.pdb"
#       where "1abc" is the pdb entry, X is the chain, and 1 is the model number
def name_pdb_file(pdb_base, chain, model, lib_path, allow_nonexist=False):
    div = pdb_base[1:3]
    path = f'{lib_path}/rcsb_pdb/{div}/{pdb_base}-{chain}-{model}.pdb'
    if not allow_nonexist and not os.path.exists(path):
        raise FileNotFoundError(f'was unable to find pdb file {path} for pdb:{pdb_base}/chain:{chain}/model:{model}.') 
    return path


def dimer2pdbs(dimer_name, lib_path):
    pdb0, chain0, model0 = dimer_name.split('-')[0].split('_')
    pdb1, chain1, model1 = dimer_name.split('-')[1].split('_')
    return name_pdb_file(pdb0, chain0, model0, lib_path), name_pdb_file(pdb1, chain1, model1, lib_path)
    

def get_models_of_chain(pdb_base, chain, lib_path):
    wild_model = "*"
    globber = name_pdb_file(pdb_base, chain, wild_model, lib_path, allow_nonexist=True)
    model_files = glob.glob(globber)
    model_nums = [re.findall(r'-(\d+?).pdb', fn)[0] for fn in model_files]
    return model_nums


