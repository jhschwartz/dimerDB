import os
import sys 
from datetime import datetime

configfile: 'config.yaml'

import tempfile
import itertools
import math

sys.path.append('scripts')
import name_pdb

sys.path.append('bin/check_contact')
from check_contact import check_contact_many_parallel
from align_tools import parallel_calc_nwalign_glocal

subworkflow_done = config['subworkflow_done']['1_derive_dimers']
intermediates = config['paths']['intermediates_dir']
lib_path = config['paths']['lib'] 

pdb_index = os.path.join(lib_path, 'pdb_index.txt')


# workflow parameters
contact_max_angstroms = config['workflow_params']['define_contact_max_dist_angstroms']
contact_min_residue_pairs = config['workflow_params']['define_contact_min_num_residue_pairs']
max_threads = config['workflow_params']['max_threads']

nw_exe = config['exe']['nwalign']


outfile = {
    'all': {
        'homodimers': os.path.join(intermediates, 'all_homodimers.txt'),
        #'heterodimers': os.path.join(intermediates, 'all_heterodimers.txt'),
#        'dimer_seq_ids': os.path.join(lib_path, 'pairs_info', 'dimer_seq_ids.tsv'),
#        'contacting_pairs': os.path.join(lib_path, 'pairs_info', 'contacting_pairs.txt'),
#        'non_contacting_pairs': os.path.join(lib_path, 'pairs_info', 'non_contacting_pairs.txt'),
    },
    'div': {
        'chains': os.path.join(intermediates, 'check_pairs', '{div}', 'chains.txt'),
        'all_chain_pairs': os.path.join(intermediates, 'check_pairs', '{div}', 'intra_assembly_chain_pairs.txt'),
        'dimers': os.path.join(intermediates, 'check_pairs', '{div}', 'dimers.txt'),
        'dimer_seq_ids': os.path.join(intermediates, 'check_pairs', '{div}', 'dimer_seq_ids.tsv'),
    }
}




divs = set()
with open(pdb_index, 'r') as f:
    for line in f:
        entry, _, _, _ = name_pdb.read_chain_names(line.rstrip())
        div = name_pdb.get_div(entry)
        divs.add(div)
divs = list(divs)




localrules: all
rule all:
    input:
        all_homodimers = outfile['all']['homodimers']
    output:
        done = subworkflow_done
    run:
        with open(output.done, 'w') as f:
            f.write('Subworkflow 1 done at ')
            f.write(str(datetime.utcnow()))
            f.write('\n')





# shadow input is pdb_index
rule write_chains:
    '''
    This rule, which is run once for each div of pdbs,
    writes a chains.txt listing all chains discovered
    from assemblies of entries of one div. This is used
    downstream for contact-checking and dimer
    categorization upon the chains of all entries of
    each div.
    '''
    output:
        chainsfile_div = outfile['div']['chains']
    #group: 'div_group'
    run:
        # pdb_index is pre-sorted
        with open(pdb_index, 'r') as indexfile, open(output.chainsfile_div, 'w') as chainsfile:
            for line in indexfile:
                chain = name_pdb.name_chain_from_filename(line.rstrip())
                entry, _, _, _ = name_pdb.read_chain_names(chain)
                div = name_pdb.get_div(entry)
                if div == wildcards.div:
                    chainsfile.write(f'{chain}\n')





rule pair_possible_contacting_chains:
    '''
    This rule, which operates separately on samples that are each
    the chains of all pdb entries of a div, computes the combinations
    of two chains within each assembly of a div. These pairs are
    the pairs to check for spatial contact, and are therefore
    possibly, but not definitely, homodimers.
    '''
    input:
        chainsfile_div = outfile['div']['chains']
    output:
        pairsfile_div = outfile['div']['all_chain_pairs']
    #group: 'div_group'
    run:
        with open(input.chainsfile_div, 'r') as f:
            chains = [line.rstrip() for line in f]

        # 1abc_a1_m1_cA -> 1abca1
        entry_assembly = lambda fullname: ''.join(name_pdb.read_chain_names(fullname)[:2])
        chains.sort(key=entry_assembly)        
        
        with open(output.pairsfile_div, 'w') as f:
            for key, group in itertools.groupby(chains, key=entry_assembly):
                chain_pairs = itertools.combinations(group, 2)
                for cp in chain_pairs:
                    potential_dimer = name_pdb.name_dimer(*cp)
                    f.write(f'{potential_dimer}\n')






# shadow input is the pdb library at {lib_path}/rcsb
rule check_contacts:
    '''
    This rule, which operates separately on samples that are each
    the set of same-assembly chain pairs of the pdb entries of
    one div, checks chain pairs for spatial contact, where contact
    is defined by at least 10 (non-exclusive) C-beta (or C-alpha if 
    C-beta not existing for a residue) atom pairs between chains 
    being less than 8 angstroms apart. 
    '''
    input:
        pairsfile_div = outfile['div']['all_chain_pairs']
    output:
        contactsfile_div = outfile['div']['dimers']
    threads: lambda wildcards, attempt: 2**attempt
    resources:
        time = lambda wildcards, attempt: '{hrs}:00:00'.format(hrs=6*attempt),
        mem_mb = lambda wildcards, attempt: str(2000 * 2**attempt)
    run:
        with open(input.pairsfile_div, 'r') as f:
            dimers = [line.rstrip() for line in f] 
          
        results = []
        if len(dimers) > 0: 
            path_pairs = []
            for d in dimers:
                c1p, c2p = name_pdb.dimer2pdbs(dimer_name=d, lib_path=lib_path)
                path_pairs.append( (c1p, c2p) ) 
            results = check_contact_many_parallel(  pairs_list=path_pairs, 
                                                    thresh_max_dist=contact_max_angstroms, 
                                                    thresh_min_pairs=contact_min_residue_pairs, 
                                                    cores=threads, 
                                                    num_series=5000                                  )
       
        with open(output.contactsfile_div, 'w') as f:
            for in_contact, dimer in zip(results, dimers):
                if in_contact: 
                    f.write(f'{dimer}\n')




rule calc_dimer_seq_ids:
    '''
    This rule, which operates separately on samples that are each
    the in-contact chains of all pdb entries of one div, 
    calcualtes the sequence identity between chains of a dimer.
    Specifically, it calculates glocal sequence identity with
    reference to the shorter chain of each dimer. 
    '''
    input:
        contactsfile_div = outfile['div']['dimers']
    output:
        dimer_seq_ids_file = outfile['div']['dimer_seq_ids']
    threads: lambda wildcards, attempt: 2**attempt
    resources:
        time = lambda wildcards, attempt: '{hrs}:00:00'.format(hrs=6*attempt),
        mem_mb = lambda wildcards, attempt: str(2000 * 2**attempt)
    run:
        with open(input.contactsfile_div, 'r') as fi:
            d2p = lambda dimer_name: name_pdb.dimer2pdbs(dimer_name, lib_path)
            pdb_pairs = ( d2p(line.rstrip()) for line in fi )
            scores = parallel_calc_nwalign_glocal(  USnw=nw_exe,
                                                    pdb_pairs=pdb_pairs,
                                                    cores=threads        )
            fi.seek(0)
            with open(output.dimer_seq_ids_file, 'w') as fo:
                for line, score in zip(fi, scores):
                    dimer = line.rstrip()
                    score = round(score, 4)
                    fo.write(f'{dimer}\t{score}\n')





rule categorize_all_dimers:
    '''
    This rule, which is run only once, uses dimer_seq_ids files
    from all divs to create a consolidated list of homodimers.
    A homodimer is here defined as an in-contact pair of chains
    with glocal sequence identity, with reference to the shorter
    chain, of >= 98%.
    '''
    input:
        dimer_seq_ids_files = expand(outfile['div']['dimer_seq_ids'], div=divs)
    output:
        all_homodimers = outfile['all']['homodimers']
    run:
        with open(output.all_homodimers, 'w') as fo:

            for dsi_file in input.dimer_seq_ids_files:

                with open(dsi_file, 'r') as fi:
                    for line in fi:
                        dimer, score = line.split()
                        if float(score) >= 0.98:
                            fo.write(f'{dimer}\n')


        shell('''

            sort -s --key=1.1,1.4 -o {output.all_homodimers} {output.all_homodimers};

                ''')
         


