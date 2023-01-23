'''
derive_all_possible_heterodimers.py - from a preprocessed yaml file that matches uniparc sequences to pdb structures
                             generates a yaml file of all possible heterodimers in the pdb.
                             Does not check that chains are in contact. Merely checks that two chains of
                             different sequences are in the same overall pdb assembly.

Written by Jacob Schwartz (jaschwa@umich.edu) in December 2022.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

These functions are unittested by test/test_derive_all_possible_heterodimers.py and passing as of 1/12/2023.
This work requires python >= 3.8
'''

import itertools
import yaml



def group_by_assembly(monomers):
    grouper = lambda t: t[1] # the pdb base name
    groups = itertools.groupby(
                    sorted(monomers, key=grouper),
                    key=grouper
                )
    result = []
    for k, v in groups:
        result.append(list(v))
    return result




def chains_as_tuples(uniparc2others):
    # a list in which we are going to place every uniparc/chain pair
    monomers = []

    # for each uniparc sequence and its matching chains
    for uniparc, values in uniparc2others.items():
        pdbs = values['pdb']

        for pdb in pdbs:
            pdb_base, chain = pdb.split('_')
            chain_tuple = (uniparc, pdb_base, chain)
            monomers.append(chain_tuple)

    return monomers


def derive_heterodimers_from_assembly_groups(groups):
    heterodimers = {}

    for monomers_of_one_assembly in groups:

        pairs = itertools.combinations(monomers_of_one_assembly, 2)

        for p in pairs:
            m0_uniparc = p[0][0]
            m1_uniparc = p[1][0]
            m0_pdb = p[0][1] + '_' + p[0][2]
            m1_pdb = p[1][1] + '_' + p[1][2]

            # only consider heterodimers, i.e. where uniparcs aren't equal
            if m0_uniparc != m1_uniparc:
                # naming convention is alphabetical by uniparc
                if m0_uniparc < m1_uniparc:
                    uniparc_pair = f'{m0_uniparc}-{m1_uniparc}'
                    pdb_pair = f'{m0_pdb}-{m1_pdb}'
                else:
                    uniparc_pair = f'{m1_uniparc}-{m0_uniparc}'
                    pdb_pair = f'{m1_pdb}-{m0_pdb}'
                
                if uniparc_pair in heterodimers:
                    heterodimers[uniparc_pair].add(pdb_pair)
                else:
                    heterodimers[uniparc_pair] = set([pdb_pair])

    for uniparc_pair, chain_pairs in heterodimers.items():
        heterodimers[uniparc_pair] = sorted(list(chain_pairs))

    return heterodimers





def heterodimers(infile, outfile):

    # open the uniparc2others yaml as a dict
    with open(infile, 'r') as f:
        uniparc2others = yaml.safe_load(f)

    monomers = chains_as_tuples(uniparc2others)

    monomers_grouped_by_pdb = group_by_assembly(monomers)

    heterodimers = derive_heterodimers_from_assembly_groups(monomers_grouped_by_pdb)

    # save to yaml
    with open(outfile, 'w') as f:
        yaml.dump(heterodimers, f, default_flow_style=None)
