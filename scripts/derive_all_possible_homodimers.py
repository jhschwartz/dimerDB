'''
derive_all_possible_homodimers.py - from a preprocessed yaml file that matches uniparc sequences to pdb structures
                             generates a yaml file of all possible homodimers in the pdb.
                             Does not check that chains are in contact. Merely checks that two chains of
                             the same sequence sequence are in the same overall pdb assembly.

Written by Jacob Schwartz (jaschwa@umich.edu) in Winter 2022-23.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

These functions are unittested by test/test_derive_all_possible_homodimers.py and passing as of 1/12/2023. NNED UPDATE 2/5/23
This work requires python >= 3.8
'''

import yaml
import itertools
import pickle
import re
from name_pdb import name_chain_from_filename
from read_fasta import read_prot_from_fasta


def expand_chains_across_assemblies_models(chains, pdb_index):
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


def group_chains(chains):
    '''
    Given a list of chains, each a string of the form "<pdb code>_a<assembly>_m<model>_c<chain>", 
    returns a list of lists, where each list is the chains of one pdb_code and one assembly.

    So if we have chains ['1abc_a1_m1_cA', '1abc_a1_m1_cB', '1abc_a2_m1_cA', '7cba_a1_m1_cB', '1abc_a2_m1_cC'],
    we should return [['1abc_a1_m1_cA', '1abc_a1_m1_cB'], ['1abc_a2_m1_cA', '1abc_a2_m1_cC'], ['7cba_B']]

    :param chains: list[str] a list of chains in pdb assemblies that correspond to a uniprot sequence
    '''
    
    chains = sorted(chains)

    # 1abc_a1_m1_cA -> 1abca1
    pdb_and_assembly = lambda fullname: ''.join(fullname.split('_')[:2])
    
    group_result = itertools.groupby(chains, key=pdb_and_assembly)
    result = []
    for k, v in group_result:
        result.append(list(v))
    return result


def derive_homodimers_from_groups(grouped_chains):
    '''
    Given a list of grouped chains (grouped by gSo if we have chains ['1abc_A', '1abc_B', '7cba_B', '1abc_C']roup_chains(..) above), returns an unnested list
    of combinations limited to combinations within individual groups.

    So if we have [['1abc_A', '1abc_B', '1abc_C'], ['7cba_B']]
    we would return  [('1abc_A', '1abc_B'), ('1abc_A', '1abc_C'), ('1abc_B', '1abc_C')]

    :param grouped_chains: list[list[str]], grouped pdb code and chain according to pdb code
    '''

    homodimers = []

    # remove groups of 1 (no homodimer in pdb assembly)
    for group in grouped_chains:
        if len(group) < 2:
            grouped_chains.remove(group)

    # nothing left to do if no groups left
    if len(grouped_chains) == 0:
        return homodimers 

    # init list of homodimers for the current prot
    for group in grouped_chains:
        combinations = itertools.combinations(group, 2)
        homodimers.append(combinations)

    # unnest the homodimers list
    homodimers = list(itertools.chain(*homodimers))

    return homodimers




def homodimers(infile, outfile, config):
    '''
    This function makes a dict, which is later saved in a yaml file, of 
    uniparc IDs matched to their mutliple homodimers. For example, if I 
    have a listing in a uniparc2others.yaml file that looks like 
        { 'UPXXX': {'pdb': ['1abc_A', '1abc_B', '7cba_B', '1abc_C'] ...} }
    then the result should be all intra-assemply unique combinations:
        {'UPXXX': [('1abc_A', '1bc_B'), ('1abc_A', '1bc_C'), ('1abc_B', '1bc_C')]}
    The result only includes keys which have any homodimer combinations 
    and the value lists only inlude PDB codes that have more than one
    identical chain. So, above, 7cba did not make the cut because it 
    had only one chain of UPXXX. Say there were another uniparc listing,
    UPZZZ, which only had entries like 7cba (only one chain of the seq).
    In this case, UPZZZ would be entirely omitted from the result.

    So if we have chains ['1abc_A', '1abc_B', '7cba_B', '1abc_C'],
    the result would be [('1abc_A', '1abc_B'), ('1abc_A', '1abc_C'), ('1abc_B', '1abc_C')]
    
    :param infile: str, the path to the uniparc2others.yaml file made by uniparc_to_uniprot_and_pdb.py
    :param outfile: str, the path where we are putting the output dict as a yaml
    '''

    lib_path = config['paths']['lib']

    # open the uniparc2others yaml as a dict
    with open(infile, 'r') as f:
        uniparc2others = yaml.safe_load(f)

    # open the pdb index
    indexfile = f'{lib_path}/rcsb_index.pkl'
    with open(indexfile, 'rb') as f:
        pdb_index = pickle.load(f)

    # the homodimers dict we are going to output
    homodimers = {}
    
    # for each uniparc sequence and its matching chains
    for uniparc, chains in [(uniparc, subdict['pdb']) for uniparc, subdict in uniparc2others.items()]:

        # read fasta and skip if too short
        div = uniparc[-2:]
        seqs_loc = config['paths']['uniparc_seqs']
        fasta = f'{seqs_loc}/{div}/{uniparc}.fasta'
        _, seq = next(read_prot_from_fasta(fasta))
        if len(seq) < config['database_settings']['chain_min_seq_len']:
            continue

        # expand list of chains to include all models
        all_chains = expand_chains_across_assemblies_models(chains, pdb_index)

        # group the chains by uniprot code, making a nested list
        grouped_chains = group_chains(all_chains)

        # make an unnested list of tuples which describe intra-assembly (same pdb code)
        # combinations of chains, e.g. ('1abc_A', '1abc_C')
        prot_homodimers = derive_homodimers_from_groups(grouped_chains)

        # convert each tuple to a combined name like '1abc_A-1abc_C'
        prot_homodimers = [f'{ph[0]}-{ph[1]}' for ph in prot_homodimers]

        # save the homodimers to the dict if there are any
        if prot_homodimers != []:
            homodimers[uniparc] = prot_homodimers

    # save to yaml
    with open(outfile, 'w') as f:
        yaml.dump(homodimers, f, default_flow_style=None)



