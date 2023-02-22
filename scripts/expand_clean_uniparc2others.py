'''
expand_uniparc2others.py - given a uniparc2others.yaml file, both expands chain
                           names to include assembly and model numbering and drops
                           any chains which do not match their uniparc seq and drops
                           uniparc seqs (and respective chains) that are too short

Written by Jacob Schwartz (jaschwa@umich.edu) in February 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu
'''

import yaml
import pickle
from read_fasta import read_prot_from_fasta
from align_tools import nw_fasta_to_pdb
from name_fasta import uniparc_fasta
from name_pdb import read_chain_names, name_pdb_file, name_chain_from_filename


def _expand_chains_across_assemblies_models(chains, pdb_index):
    '''
    Given a list of chains, each a string of the form "<pdb code>_<chain id>",
    returns a list of those chains with model numbers and assembly identifies, 
    duplicating chain names that occur across different models and/or assemblies.

    So if we have chains ['1abc_A', '1abc_B', '7cba_B', '1abc_C'] we might get
    something like ['1abc_a1_m1_cA', '1abc_a1_m2_cB', ...] 
    
    :param chains: list[str], a list of chains in pdb assemblies that correspond to a uniprot sequence
    :param pdb_index: dict, wherein keys are pdb entry codes, second-depth keys are chain names, and
                        third-depth members are lists of which each member is a pdb filename.
                            E.g.:
                                {
                                    '1abc': {
                                        'A': ['1abc-a1-m1-cA.pdb', ...'],
                                        ...
                                    },
                                    ...
                                }
    '''
    all_chains = []
    for chain in chains:
        pdb_base, chain_ID = chain.split('_')
        try:
            chain_files = pdb_index[pdb_base][chain_ID]
        except KeyError:
            print(f'Warning: the chain {pdb_base}:{chain_ID} does not exist in the local pdb and has been skipped while creating possible homodimers list') 
            continue
        for cf in chain_files:
            all_chains.append(name_chain_from_filename(cf))
    return all_chains





def _compare_uniparc_to_chain(uniparc_id, chain_name, config):
    lib = config['paths']['lib']
    nw = config['paths']['nwalign']
    pdb2fasta = config['paths']['pdb2fasta']

    uniparc_fastafile = uniparc_fasta(uniparc_id)
    pdb_file = name_pdb_file(*read_chain_names(chain_name), lib)

    score = nw_fasta_to_pdb(uniparc_fastafile, pdb_file, nw, pdb2fasta)
    return score




def expand_clean_uniparc2others(inyaml, outyaml, config):
    lib_path = config['paths']['lib']
    seqmatch_thresh = config['database_settings']['uniparc_chain_seqmatch_id_thresh']


    # open the uniparc2others yaml as a dict
    with open(inyaml, 'r') as f:
        uniparc2others = yaml.safe_load(f)

    # open the pdb index
    indexfile = f'{lib_path}/rcsb_index.pkl'
    with open(indexfile, 'rb') as f:
        pdb_index = pickle.load(f)

    # the homodimers dict we are going to output
    homodimers = {}

    for uniparc, entry in uniparc2others.items():
        div = uniparc[-2:]
        seqs_loc = config['paths']['uniparc_seqs']
        fasta = f'{seqs_loc}/{div}/{uniparc}.fasta'
        _, seq = next(read_prot_from_fasta(fasta))
        if len(seq) < config['database_settings']['chain_min_seq_len']:
            continue

        chains = entry['pdb']
        expanded_chains = _expand_chains_across_assemblies_models(chains, pdb_index)

        matching_chains = []
        for ec in expanded_chains:
            score = _compare_uniparc_to_chain(uniparc, ec, config)
            if score >= seqmatch_thresh:
                matching_chains.append(ec)
            else:
                print(f'encountered non-matching chain: uniparc={uniparc}, chain={ec}, id={score}')

        if len(matching_chains) > 0:
            homodimers[uniparc] = matching_chains

    with open(outyaml, 'w') as f:
        yaml.dump(homodimers, f)
