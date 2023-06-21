import os
import sys 
from datetime import datetime

configfile: 'config.yaml'


import tempfile
import itertools
import math
import shutil
from sortedcontainers import SortedList, SortedSet
import subprocess

sys.path.append('scripts')
import name_pdb
from contact_database import ContactDB

sys.path.append('bin/check_contact')
from check_contact import check_contact_many_parallel
from align_tools import parallel_calc_nwalign_glocal

subworkflow_done = config['subworkflow_done']['1_derive_dimers']
intermediates = config['paths']['intermediates_dir']
lib_path = config['paths']['lib'] 

contacts_db_path = config['paths']['contacts_db']
seqid_db_path = config['paths']['seqid_db']


pdb_index = os.path.join(lib_path, 'pdb_index.txt')


# workflow parameters
contact_max_angstroms = config['workflow_params']['define_contact_max_dist_angstroms']
contact_min_residue_pairs = config['workflow_params']['define_contact_min_num_residue_pairs']
max_threads = config['workflow_params']['max_threads']

nw_exe = config['exe']['nwalign']
check_contact_exe = config['exe']['check_contact']

outfile = {
    'all': {
        'homodimers': os.path.join(intermediates, 'all_homodimers.txt'),
        'contacts_saved': os.path.join(intermediates, 'sub1', 'contacts_saved.done')
        #'heterodimers': os.path.join(intermediates, 'all_heterodimers.txt'),
#        'dimer_seq_ids': os.path.join(lib_path, 'pairs_info', 'dimer_seq_ids.tsv'),
#        'contacting_pairs': os.path.join(lib_path, 'pairs_info', 'contacting_pairs.txt'),
#        'non_contacting_pairs': os.path.join(lib_path, 'pairs_info', 'non_contacting_pairs.txt'),
    },
    'div': {
        'chains': os.path.join(intermediates, 'check_pairs', '{div}', 'chains.txt'),
        'all_chain_pairs': os.path.join(intermediates, 'check_pairs', '{div}', 'intra_assembly_chain_pairs.txt'),
        'contacts_lookup': os.path.join(intermediates, 'check_pairs', '{div}', 'contacts_lookup.txt'),
        'contacts_needcalc': os.path.join(intermediates, 'check_pairs', '{div}', 'contacts_needcalc.txt'),
        'contacts_calc': os.path.join(intermediates, 'check_pairs', '{div}', 'contacts_calc.txt'),
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
        all_homodimers = outfile['all']['homodimers'],
        saved_contacts_done = outfile['all']['contacts_saved']
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
        chains = []

        with open(pdb_index, 'r') as indexfile:
            for line in indexfile:
                chain = name_pdb.name_chain_from_filename(line.rstrip())
                entry, _, _, _ = name_pdb.read_chain_names(chain)
                div = name_pdb.get_div(entry)
                if div == wildcards.div:
                    chains.append(chain)

        chains.sort(key=name_pdb.name_chain_from_filename)

        with open(output.chainsfile_div, 'w') as chainsfile:
            for chain in chains:
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



rule lookup_contacts:
    input:
        pairsfile_div = outfile['div']['all_chain_pairs']
    output:
        contactsfound_div = outfile['div']['contacts_lookup']
    threads: 1
    resources:
        time = '3:00:00',
        mem_mb = lambda wildcards, attempt: str(2000 * 2**attempt)
    retries: 4
    run:
        db = ContactDB(contacts_db_path)
        
        with open(input.pairsfile_div, 'r') as fi:
            pairs = [ name_pdb.dimer2chains(line.strip()) for line in fi ]
            
            with open(output.contactsfound_div, 'w') as fo:
                for contact_found in db.get(pairs):
                    c1, c2, is_contact = contact_found
                    dimer = name_pdb.name_dimer(c1, c2)
                    fo.write(f'{dimer}\t{is_contact}\n')
        


rule pairs_needing_contact_calc:
    input:
        pairsfile_div = outfile['div']['all_chain_pairs'],
        contactsfound_div = outfile['div']['contacts_lookup']
    output:
        contactsneedcalc_div = outfile['div']['contacts_needcalc']
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: str(1000 * 2**attempt)
    retries: 4
    shell:
        # ok to use GNU sort here because it's just for generating the need_calc file,
        # which can be used without sorting and which gets its results sorted later
        '''
        tmp_not_need_calc=$(mktemp /tmp/XXXXXXXX);
        tmp_all_pairs=$(mktemp /tmp/XXXXXXXX);
        awk '{{print $1}}' {input.contactsfound_div} > $tmp_not_need_calc; 
        sort -s  -o $tmp_not_need_calc $tmp_not_need_calc;
        sort -s  -o $tmp_all_pairs {input.pairsfile_div}; 
        comm -23 --check-order $tmp_all_pairs $tmp_not_need_calc > {output.contactsneedcalc_div};
        '''



# shadow input is the pdb library at {lib_path}/rcsb
rule calc_contacts:
    '''
    This rule, which operates separately on samples that are each
    the set of same-assembly chain pairs of the pdb entries of
    one div, checks chain pairs for spatial contact, where contact
    is defined by at least 10 (non-exclusive) C-beta (or C-alpha if 
    C-beta not existing for a residue) atom pairs between chains 
    being less than 8 angstroms apart. 
    '''
    input:
        contactsneedcalc_div = outfile['div']['contacts_needcalc']
    output:
        contactscalc_div = outfile['div']['contacts_calc']
    threads: 1
    resources:
        time = '2-00:00:00',
        mem_mb = '2000' 
    run:
        # special TMPDIR settings for SDSC Expanse, if that is the current machine
        jobid = os.environ.get('SLURM_JOB_ID', None)
        user = os.environ['USER']
        potentialtmp = f'/scratch/{user}/job_{jobid}'
        use_tmpdir_base = None
        if jobid and os.path.exists(potentialtmp):
            use_tmpdir_base = potentialtmp
            print(f'using TMPDIR {potentialtmp} for python tempfile')
        else:
            print('unable to use custom TMPDIR for python tempfile')


        
        with tempfile.TemporaryDirectory(dir=use_tmpdir_base) as tempdir:

            templib = os.path.join(tempdir, 'lib')
            temprcsb = os.path.join(templib, 'rcsb')
            os.makedirs(temprcsb)

            # copy div to node
            source_div = os.path.join(lib_path, 'rcsb', wildcards.div)
            dest_div = os.path.join(temprcsb, wildcards.div)
            shutil.copytree(source_div, dest_div)

            check_contact_infile = os.path.join(tempdir, 'pairs.txt')
            check_contact_outfile = os.path.join(tempdir, 'contact_info.txt')
            
            with open(check_contact_infile, 'w') as fo, \
                                        open(input.contactsneedcalc_div, 'r') as fi:
                
                for line in fi:
                    if line == '':
                        continue
                    chain1name, chain2name = name_pdb.dimer2chains(line.strip())
                    fo.write(f'{chain1name}\t{chain2name}\n')
    

            # if pairsfile is not empty
            if os.stat(check_contact_infile).st_size > 0:
                cmd = f'{check_contact_exe} {check_contact_infile} {check_contact_outfile} 8 10 {temprcsb}'
                exe_result = subprocess.run(cmd.split())
                if exe_result.returncode != 0:
                    raise RuntimeError(f'encountered error from {check_contact_exe}: {exe_result.stderr} {exe_result.stdout}') 

                
                with open(output.contactscalc_div, 'w') as fout, \
                                        open(check_contact_outfile, 'r') as fin_result, \
                                        open(input.contactsneedcalc_div, 'r') as fin_name:
                    for dimer, result in zip(fin_name, fin_result):
                        if result.split()[0] != '0':
                            fout.write(f'{dimer.strip()}\n')
    

            else:
                shell(''' touch {output.contactscalc_div} ''')
            


rule write_dimers:
    input:
        contactsfound_div = outfile['div']['contacts_lookup'],
        contactscalc_div = outfile['div']['contacts_calc']
    output:
        dimersfile_div = outfile['div']['dimers']
    threads: 1
    run:
        sort_by_first_chain = lambda dimer: name_pdb.dimer2chains(dimer)
        dimers = SortedList(key=sort_by_first_chain)

        with open(input.contactsfound_div, 'r') as f_found:
            for line in f_found:
                dimer, is_contact = line.split()
                is_contact = is_contact == 'True'
                if is_contact:
                    dimers.add(dimer)

        with open(input.contactscalc_div, 'r') as f_calc:
            for line in f_calc:
                dimer = line.strip()
                dimers.add(dimer)

        with open(output.dimersfile_div, 'w') as fo:
            for dimer in dimers:
                fo.write(f'{dimer}\n')




rule calc_dimer_seq_ids:
    '''
    This rule, which operates separately on samples that are each
    the in-contact chains of all pdb entries of one div, 
    calcualtes the sequence identity between chains of a dimer.
    Specifically, it calculates glocal sequence identity with
    reference to the shorter chain of each dimer. 
    '''
    input:
        dimersfile_div = outfile['div']['dimers']
    output:
        dimer_seq_ids_file = outfile['div']['dimer_seq_ids']
    threads: lambda wildcards, attempt: 2**attempt
    retries: 5
    resources:
        time = lambda wildcards, attempt: '{hrs}:00:00'.format(hrs=6*attempt),
        mem_mb = lambda wildcards, attempt: str(2000 * 2**attempt)
    run:
        if os.stat(input.dimersfile_div).st_size == 0:
            shell(''' touch {output.dimer_seq_ids_file} ''')

        else:
            with open(input.dimersfile_div, 'r') as fi:
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



rule save_contacts:
    input:
        need_calc_files = expand(outfile['div']['contacts_needcalc'], div=divs),
        calc_contacing_files = expand(outfile['div']['contacts_calc'], div=divs)
    output:
        done = outfile['all']['contacts_saved']
    run:
        db = ContactDB(contacts_db_path)

        for need_calc_file, calc_file in zip(input.need_calc_files, input.calc_contacing_files):
            is_contact = SortedSet()
            data_to_save = []
            with open(calc_file, 'r') as f:
                for line in f:
                    is_contact.add(line.strip())
            with open(need_calc_file, 'r') as f:
                for line in f:
                    dimer = line.strip()
                    c1, c2 = name_pdb.dimer2chains(dimer)
                    contacting = dimer in is_contact
                    db.auto_buffer( (c1, c2, contacting) )

        db.update_clear_buffer()

        shell(''' touch {output.done} ''')



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


