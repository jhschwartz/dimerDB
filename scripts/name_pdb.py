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



def read_chain_names(name):
    # filename
    if re.match(r'.*[0-9a-z]{4}-a[0-9]+-m[0-9]+-c[A-Za-z0-9\-]+\.pdb', name):
        pdb = re.findall(r'.*([0-9a-z]{4})-a[0-9]+-m[0-9]+-c[A-Za-z0-9\-]+\.pdb', name)[0]
        assembly = re.findall(r'.*[0-9a-z]{4}-a([0-9]+)-m[0-9]+-c[A-Za-z0-9\-]+\.pdb', name)[0] 
        model = re.findall(r'.*[0-9a-z]{4}-a[0-9]+-m([0-9]+)-c[A-Za-z0-9\-]+\.pdb', name)[0]
        chain = re.findall(r'.*[0-9a-z]{4}-a[0-9]+-m[0-9]+-c([A-Za-z0-9\-]+)\.pdb', name)[0] 
        return pdb, assembly, model, chain

    # just a name
    elif re.match(r'[0-9a-z]{4}-a[0-9]+-m[0-9]+-c[A-Za-z0-9\-]+$', name):
        pdb = re.findall(r'([0-9a-z]{4})-a[0-9]+-m[0-9]+-c[A-Za-z0-9\-]+', name)[0]
        assembly = re.findall(r'[0-9a-z]{4}-a([0-9]+)-m[0-9]+-c[A-Za-z0-9\-]+', name)[0] 
        model = re.findall(r'[0-9a-z]{4}-a[0-9]+-m([0-9]+)-c[A-Za-z0-9\-]+', name)[0]
        chain = re.findall(r'[0-9a-z]{4}-a[0-9]+-m[0-9]+-c([A-Za-z0-9\-]+)', name)[0]
        return pdb, assembly, model, chain

    # invalid
    else:
        raise ValueError(f'invalid chain name encoutnered, neither a filename nor a simple name: {name}')


def get_div(name):
    if len(name) != 4:
        pdb, _, _, _ = read_chain_names(name)
    else:
        pdb = name
    return pdb[1:3]


def name_chain_from_filename(filename):
    pdb_base, assembly, model, chain = read_chain_names(filename)
    return f'{pdb_base}-a{assembly}-m{model}-c{chain}'



def name_pdb_file(pdb_base, assembly, model, chain, lib_path=None, allow_nonexist=False):
    div = get_div(pdb_base)
    basename = f'{pdb_base}-a{assembly}-m{model}-c{chain}.pdb'
    if lib_path:
        path = os.path.join(lib_path, 'rcsb', div, basename)
        if not allow_nonexist and not os.path.exists(path):
            raise FileNotFoundError(f'was unable to find pdb file {path} for pdb:{pdb_base}, assembly:{assembly}, model:{model}, chain:{chain}') 
        return path
    return basename


def dimer2pdbs(dimer_name, lib_path):
    pdb0, assembly0, model0, chain0 = read_chain_names(dimer_name.split('_')[0])
    pdb1, assembly1, model1, chain1 = read_chain_names(dimer_name.split('_')[1])
    return name_pdb_file(pdb0, assembly0, model0, chain0, lib_path), name_pdb_file(pdb1, assembly1, model1, chain1, lib_path)


def dimer2chains(dimer_name):
    return dimer_name.split('_')


def name_dimer(chain1name, chain2name):
    c1 = chain1name
    c2 = chain2name
    if chain2name < chain1name:
        c1 = chain2name
        c2 = chain1name
    return f'{c1}_{c2}'






