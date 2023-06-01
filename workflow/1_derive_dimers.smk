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
from align_tools import calc_nwalign_glocal

subworkflow_done = config['workflow']['1_derive_dimers_done']
intermediates = config['paths']['intermediates_dir']
lib_path = config['paths']['lib'] # TODO put in config, match to sub0, must contain folder "rcsb"

pdb_index = os.path.join(lib_path, 'pdb_index.txt')


# workflow parameters
contact_max_angstroms = config['workflow_params']['define_contact_max_dist_angstroms']
contact_min_residue_pairs = config['workflow_params']['define_contact_min_num_residue_pairs']
max_threads = config['workflow_params']['max_threads']

outfile = {
    'all': {
        'homodimers': os.path.join(intermediates, 'check_pairs', 'all_homodimers.txt')
        #'heterodimers': os.path.join(intermediates, 'check_pairs', 'all_heterodimers.txt') 
    },
    'entry_template': {
        'chains': intermediates+'/check_pairs/div/{div}/{entry}/chains.txt',
        'all_chain_pairs': intermediates+'/check_pairs/div/{div}/{entry}/all_intra_assembly_chain_pairs.txt',
        'contacting_pairs': intermediates+'/check_pairs/div/{div}/{entry}/contacting_chain_pairs.txt',
        'homodimers': intermediates+'/check_pairs/div/{div}/{entry}/homodimers.txt'
        #'heterodimers': intermediates+'/check_pairs/div/{div}/{entry}/heterodimers.txt' 
    }
}



### define samples via entries and divs
entries_counter = {}
with open(pdb_index, 'r') as f:
    for line in f:
        entry, _, _, _ = name_pdb.read_chain_names(line.rstrip())
        if entry in entries_counter:
            entries_counter[entry] += 1
        else:
            entries_counter[entry] = 1
entries = list(entries_counter.keys())
divs = [name_pdb.get_div(e) for e in entries]




rule all:
    input:
        all_homodimers = outfile['all']['homodimers']
    output:
        done = subworkflow_done
    run:
        with open(output.done, 'w') as f:
            f.write('Subworkflow 1 done at ')
            f.write(str(datetime.utcnow()))





# shadow input is pdb_index
rule write_chains:
    '''
    This rule, which is run only once, writes a chain.txt
    for each pdb entry in pdb_index, listing all chains
    discovered from pdb assemblies of that entry. This
    is used for downstream contact-checking and dimer 
    categorization upon the chains of each entry.
    '''
    output:
        chainsfile = expand(outfile['entry_template']['chains'], zip, div=divs, entry=entries)
    run:
        # pdb_index is pre-sorted
        with open(pdb_index, 'r') as f: 
            last_entry = None
            ofile = None
            for line in f:
                chain = name_pdb.name_chain_from_filename(line.rstrip())
                entry, _, _, _ = name_pdb.read_chain_names(chain)
                div = name_pdb.get_div(entry)
                if entry != last_entry:
                    last_entry = entry
                    if ofile and not ofile.closed:
                        ofile.close()
                    ofile = open(outfile['entry_template']['chains'].format(div=div, entry=entry), 'w')
                ofile.write(f'{chain}\n')
            if ofile and not ofile.closed:
                ofile.close()
                 





rule pair_possible_contacting_chains:
    '''
    This rule, which operates separately on samples that are each
    the chains of one pdb entry, computes the combinations of
    two chains within each assembly of a pdb entry. These pairs
    are the pairs to check for spatial contact, and are therefore
    possibly, but not definitely, homodimers.
    '''
    input:
        chainsfile = outfile['entry_template']['chains']
    output:
        pairsfile = outfile['entry_template']['all_chain_pairs']
    run:
        with open(input.chainsfile, 'r') as f:
            chains = [line.rstrip() for line in f]

        # 1abc_a1_m1_cA -> 1abca1
        entry_assembly = lambda fullname: ''.join(name_pdb.read_chain_names(fullname)[:2])
        chains.sort(key=entry_assembly)        
        
        with open(output.pairsfile, 'w') as f:
            for key, group in itertools.groupby(chains, key=entry_assembly):
                chain_pairs = itertools.combinations(group, 2)
                for cp in chain_pairs:
                    potential_dimer = name_pdb.name_dimer(*cp)
                    f.write(f'{potential_dimer}\n')





# shadow input is the pdb library at {lib_path}/rcsb
rule check_contacts:
    '''
    This rule, which operates separately on samples that are each
    the set of same-assembly chain pairs of one pdb entry, checks
    chain pairs for spatial contact, where contact is defined by
    at least 10 (non-exclusive) C-beta (or C-alpha if C-beta not
    existing for a residue) atom pairs between chains being less
    than 8 angstroms apart. 
    '''
    input:
        pairsfile = outfile['entry_template']['all_chain_pairs']
    output:
        contactsfile = outfile['entry_template']['contacting_pairs']
    threads:
        # set threads based on number of chains in entry
        lambda wildcards: min(math.ceil(entries_counter[wildcards.entry] / 1000), max_threads)
    run:
        with open(input.pairsfile, 'r') as f:
            dimers = [line.rstrip() for line in f] 
           
        path_pairs = []
        for d in dimers:
            c1p, c2p = name_pdb.dimer2pdbs(dimer_name=d, lib_path=lib_path)
            path_pairs.append( (c1p, c2p) ) 
        results = check_contact_many_parallel(  pairs_list=path_pairs, 
                                                thresh_max_dist=contact_max_angstroms, 
                                                thresh_min_pairs=contact_min_residue_pairs, 
                                                cores=threads, 
                                                num_series=1000                                  )
       
        with open(output.contactsfile, 'w') as f:
            for in_contact, dimer in zip(results, dimers):
                if in_contact: 
                    f.write(f'{dimer}\n')




rule categorize_dimers:
    '''
    This rule, which operates separately on samples that are each
    the in-contact chains of one pdb entry, checks sequence 
    identity between chains of a dimer to categorize each dimer
    as either a homodimer or heterodimer. A homodimer is defined
    as a dimer whose sequences have >= 98% glocal (both) sequence
    identity with reference to the shorter of the two chains.
    '''
    input:
        contactsfile = outfile['entry_template']['contacting_pairs']
    output:
        homodimersfile = outfile['entry_template']['homodimers']
    threads: 1
    run:
        with open(input.contactsfile, 'r') as fi, \
                open(output.homodimersfile, 'w') as fo:
            for line in fi:
                dimer = line.rstrip()
                pdb1, pdb2 = name_pdb.dimer2pdbs(dimer_name=dimer, lib_path=lib_path)
                score = calc_nwalign_glocal('bin/USalign/NWalign', pdb1, pdb2)
                if score >= 0.98:
                    fo.write(line)




rule unify_homodimers:
    '''
    This rule, which is run only once, combines the homodimers found
    across pdb entries into one unified homodimers file.
    '''
    input:
        homodimer_files = expand(outfile['entry_template']['homodimers'], zip, div=divs, entry=entries)
    output:
        all_homodimers = outfile['all']['homodimers']
    shell:
        '''
        for hf in {input.homodimer_files};
        do
            cat $hf >> {output.all_homodimers};
        done

        sort -s --key=1.1,1.4 -o {output.all_homodimers} {output.all_homodimers};
        '''

