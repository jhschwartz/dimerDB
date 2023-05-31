'''
derive_all_possible_homodimers.py - from a preprocessed pickle file that matches uniparc sequences to pdb structures
                             generates a pickle file of all possible homodimers in the pdb.
                             Does not check that chains are in contact. Merely checks that two chains of
                             the same sequence sequence are in the same overall pdb assembly.

Written by Jacob Schwartz (jaschwa@umich.edu) in Winter 2022-23.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu
'''

import pickle
import itertools


def _group_chains(chains):
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


def _derive_homodimers_from_groups(grouped_chains):
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




def derive_homodimers(infile, outfile):
    '''
    This function makes a dict, which is later saved in a pickle file, of 
    uniparc IDs matched to their mutliple homodimers. For example, if I 
    have a listing in a uniparc2others.pkl file that looks like 
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
    
    :param infile: str, the path to the uniparc2others.pkl file made by uniparc_to_uniprot_and_pdb.py
    :param outfile: str, the path where we are putting the output dict as a pickle
    '''

    # open the uniparc2others pickle as a dict
    with open(infile, 'rb') as f:
        uniparc2others = pickle.load(f)
 
    # the homodimers dict we are going to output
    homodimers = {}
    
    # for each uniparc sequence and its matching chains
    for uniparc, chains in uniparc2others.items():

        # group the chains by uniprot code, making a nested list
        grouped_chains = _group_chains(chains)

        # make an unnested list of tuples which describe intra-assembly (same pdb code)
        # combinations of chains, e.g. ('1abc_A', '1abc_C')
        prot_homodimers = _derive_homodimers_from_groups(grouped_chains)

        # convert each tuple to a combined name like '1abc_A-1abc_C'
        prot_homodimers = [f'{ph[0]}-{ph[1]}' for ph in prot_homodimers]

        # save the homodimers to the dict if there are any
        if prot_homodimers != []:
            homodimers[uniparc] = prot_homodimers

    # save to pickle
    with open(outfile, 'wb') as f:
        pickle.dump(homodimers, f)
