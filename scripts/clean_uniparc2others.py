'''
clean_uniparc2others.py - functions to filter the uniparc2others_expanded.yaml file, clearing out PDBs of a given uniparc if the pdb sequence does not appear to match the uniparc sequence. 

Written by Jacob Schwartz (jaschwa@umich.edu) in February 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu
'''

import yaml
from multiprocessing import Pool
import itertools

from align_tools import nw_fasta_to_pdb
from name_fasta import uniparc_fasta
from name_pdb import read_chain_names, name_pdb_file




def check_uniparc_matches_chain(uniparc_name, chain_name, config):
    lib = config['paths']['lib']
    nw = config['paths']['nwalign']
    pdb2fasta = config['paths']['pdb2fasta']
    
    uniparc_fastafile = uniparc_fasta(uniparc_name)
    pdb_file = name_pdb_file(*read_chain_names(chain_name), lib)
    
    score = nw_fasta_to_pdb(uniparc_fastafile, pdb_file, nw, pdb2fasta)
    return (uniparc_name, chain_name, score)



def check_structures_of_uniparc(args):
    uniparc, values, config = args
    chains = values['pdb']
    nw = config['paths']['nwalign']
    results = [check_uniparc_matches_chain(uniparc, chain, config) for chain in chains]
    return results


def clean_uniparc2others(inyaml, outyaml, config):
    with open(inyaml, 'r') as f:
        uniparc2others_data = yaml.safe_load(f)

    with Pool(num_processes=config['runtime']['max_threads']):
        args = ((k, v, config) for k, v in uniparc2others_data.items())
        results = Pool.imap(check_structures_of_uniparc, args)

    for (uniparc, pdb, score) in itertools.chain.from_iterable(results):
        if score < 0.95:
             print(f'Uniparc {uniparc} does not match chain {pdb}, nw = {score} < 0.95')
        else:
            if uniparc not in outdata:
                outdata[uniparc] = []
            if not pdb in outdata[uniparc]:
                outdata[uniparc].append(pdb)

    with open(outyaml, 'w') as f:
        yaml.dump(outdata, f)

