import os
import sys
import tempfile
from sortedcontainers import SortedSet, SortedList
from datetime import datetime

configfile: 'config.yaml'

sys.path.append('scripts')
import name_pdb
import read_fasta 
from wrap_mmseqs import run_mmseqs_cluster, derive_mmseqs_reps_from_tsv

subworkflow_done = config['subworkflow_done']['2_seq_cluster']
intermediates = config['paths']['intermediates_dir']
lib_path = config['paths']['lib'] 

pdb_index = os.path.join(lib_path, 'pdb_index.txt')
all_homodimers = os.path.join(intermediates, 'all_homodimers.txt')
subflow1_div_dimer_seq_ids = os.path.join(intermediates, 'check_pairs', '{div}', 'dimer_seq_ids.tsv')

if config.get('test'):
    all_homodimers = os.path.join(lib_path, 'test_in_homodimers.txt')
    
max_threads = config['workflow_params']['max_threads']
seq_cluster_id = config['workflow_params']['seq_cluster_id']

mmseqs_exe = config['exe']['mmseqs']

outfile = {
    'fasta': {
        'index': os.path.join(lib_path, 'fasta_index.tsv'),
        'homodimers_seqs': os.path.join(intermediates, 'homodimer_seqs.fasta'),
        #'all_seqs': os.path.join(lib_path, 'all_seqs.fasta'),
        'single_seq_template': os.path.join(lib_path, 'fasta', '{div}', '{chain}.fasta'),
        'div_template': os.path.join(lib_path, 'fasta', 'homodimers_{div}.fasta')
    },
    'cluster': {
        'index': os.path.join(intermediates, 'cluster', 'cluster_index.txt'),
        'cluster_members_dimers_template': os.path.join(intermediates, 'cluster', 
                                            'seq_clusters', '{cluster_name}', 'members_dimers.txt'),
        'cluster_members_chains_template': os.path.join(intermediates, 'cluster', 
                                            'seq_clusters', '{cluster_name}', 'members_chains.txt'),
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
#### BEGIN RULE ALL TARGETS
    input:
        clusterindex = outfile['cluster']['index'],
        fastaindex = outfile['fasta']['index']
    output:
        done = subworkflow_done
    run:
        with open(output.done, 'w') as f:
            f.write('Subworkflow 2 done at ')
            f.write(str(datetime.utcnow()))
            f.write('\n')
#### END RULE ALL TARGETS
# defining begin/end of rule all in comments is necessary for unittesting





rule write_div_fasta:
    input:
        dsi_file = subflow1_div_dimer_seq_ids
    output:
        fasta = outfile['fasta']['div_template']
    run:
        with tempfile.NamedTemporaryFile('w+t') as tf:
        
            div_homodimerizing_chains = SortedSet()
            with open(input.dsi_file, 'r') as fi:
                for line in fi:
                    dimer, score = line.split()
                    if float(score) >= 0.98:
                        c1, _ = name_pdb.dimer2chains(dimer)
                        div_homodimerizing_chains.add(c1)

            for chain in div_homodimerizing_chains:
                chainfile = name_pdb.name_pdb_file(*name_pdb.read_chain_names(chain))
                line = f'{wildcards.div}/{chainfile}\n'
                tf.write(line)
            
            tf.seek(0)

            rcsb_dir = os.path.join(lib_path, 'rcsb/')

            shell(
                ''' 
                set +e;
                bin/USalign/pdb2fasta -mol prot -het 2 -dir {rcsb_dir} {tf.name} > {output.fasta}; 
                sed -i -r "s/>.*\//>/g" {output.fasta};
                sed -i -r "s/\.pdb\S*//g" {output.fasta};
                '''
            )  # resulting headers like ">ql/6qle-a1-m1-c1.fasta 123", where 123 is seq length
            # TODO: move seds to new rule! 
            # instead of second sed:
                # cut -d. -f1 input.fasta > output.fasta 
            # instead of first sed: cd into the div directory and avoid it entirely




rule combine_all_fastas:
    input:
        fastas = expand(outfile['fasta']['div_template'], div=divs)
    output:
        fasta = outfile['fasta']['homodimers_seqs']
    shell:
        '''
        for fastafile in {input.fastas};
        do
            cat $fastafile;
        done > {output.fasta};
        '''



rule write_index_individual_fastas:
    input:
        div_fastas = expand(outfile['fasta']['div_template'], div=divs)
    output:
        index = outfile['fasta']['index']
    run:
        with open(output.index, 'w') as outindex:

            for fasta in input.div_fastas:

                for header, seq in read_fasta.read_prot_from_fasta(fasta):
                    if not header:
                        print('WTF', fasta, header, seq)
                        break
                    chain = header.split()[0].lstrip('>')
                    div = name_pdb.get_div(chain)
                    seq_fasta = outfile['fasta']['single_seq_template'].format(div=div, chain=chain)
                    dir_ = os.path.dirname(seq_fasta)
                    os.makedirs(dir_, exist_ok=True)
                    L = len(seq)
                    with open(seq_fasta, 'w') as f:
                        f.write(f'{header}\n{seq}\n')

                    outindex.write(f'{div}/{chain}.fasta\t{L}\n')

        # sort the index file - the "-s", "-t" and "-k" sort only by letters before ".pdb"
        shell(''' sort -s -t "." -k 1,1 -o {output.index} {output.index} ''')


#
## shadow input pdb_index
#rule process_fastas:
#    '''
#    This rule, which runs only once without samples, runs
#    pdb2fasta to generate one fasta for the entire set of
#    pdb chains found to homodimerize by subflow1.
#    This large fasta then has its headers reformatted,
#    then individual fasta files for each sequence are 
#    written to {lib_path}/fasta. Additionally, an index
#    of all these fastas is written to {fasta_index}.
#    '''
#    output:
#        all_seqs_fasta = outfile['fasta']['all_seqs'],
#        indexfile_tsv = outfile['fasta']['index']
#    resources:
#        time = '48:00:00',
#        mem_mb = '10000'
#    run:
#        rcsb_dir = os.path.join(lib_path, 'rcsb/')
#        
#        # write all pdbs to one fasta and edit its headers
#        shell(
#            ''' 
#            set +e;
#            bin/USalign/pdb2fasta -mol prot -het 2 -dir {rcsb_dir} {pdb_index} > {output.all_seqs_fasta}; 
#            sed -i -r "s/>.*\//>/g" {output.all_seqs_fasta};
#            sed -i -r "s/\.pdb\S*//g" {output.all_seqs_fasta};
#            '''
#        )  # resulting headers like ">ql/6qle-a1-m1-c1.fasta 123", where 123 is seq length
#
#        # write each seq to its own fasta file and keep index
#        with open(output.indexfile_tsv, 'w') as indexf:
#            for header, seq in read_fasta.read_prot_from_fasta(output.all_seqs_fasta):
#                chain = header.lstrip('>').split()[0]
#                L = header.rstrip().split()[1]
#                div = name_pdb.get_div(chain)
#                outfasta = outfile['fasta']['single_seq_template'].format(div=div, chain=chain) 
#                os.makedirs(os.path.dirname(outfasta), exist_ok=True)
#                with open(outfasta, 'w') as fo:
#                    fo.write(f'{header}\n{seq}\n')
#                indexf.write(f'{div}/{chain}.fasta\t{L}\n')
#
#
#        # sort the index file - the "-s", "-t" and "-k" sort only by letters before ".pdb"
#        shell(''' sort -s -t "." -k 1,1 -o {output.indexfile_tsv} {output.indexfile_tsv} ''')
#
#
#
## shadow input is {all_homodimers}
#rule write_homodimers_chains_fastas:
#    input:
#        all_seqs_fasta = outfile['fasta']['all_seqs']
#    output:
#        fasta = outfile['fasta']['homodimers_seqs']
#    resources:
#        time = '48:00:00',
#        mem_mb = '5000'
#    run:
#        with open(all_homodimers, 'r') as f:
#            homodimers = [line.strip() for line in f]
#        
#        homodimer_chains = SortedSet()
#        for dimer in homodimers:
#            c1, _ = name_pdb.dimer2chains(dimer)
#            homodimer_chains.add(c1)
#
#        with open(output.fasta, 'w') as f:
#            for header, seq in read_fasta.read_prot_from_fasta(input.all_seqs_fasta):
#                chain = header.split()[0].lstrip('>')
#                if chain in homodimer_chains:
#                    f.write(f'{header}\n{seq}\n')
#
#


# shadow input is {all_homodimers} 
rule run_seq_cluster:
    input:
        fasta = outfile['fasta']['homodimers_seqs'],
    output:
        cluster_index = outfile['cluster']['index']
    threads: 1
    resources:
        mem_mb = '2000',
        time = '48:00:00'
    run:
        with open(all_homodimers, 'r') as f:
            homodimers = [line.strip() for line in f]

        with tempfile.TemporaryDirectory() as td:

            # run mmseqs and write its output files to tempdir
            outpref = os.path.join(td, 'mmseqs_out')
            mmseqs_cluster_tsv, _ = run_mmseqs_cluster(  infasta=input.fasta, 
                                                         outprefix=outpref,
                                                         seq_id=seq_cluster_id, 
                                                         mmseqs_exe=mmseqs_exe, 
                                                         cores=threads          )

            print('1: mmseqs done', flush=True)

            # read membership relationships from mmseqs output
            membership = {}
            with open(mmseqs_cluster_tsv, 'r') as f:
                for line in f:
                    chain_rep, chain_mem = line.rstrip().split()
                    if not chain_rep in membership:
                        membership[chain_rep] = SortedList()
                    membership[chain_rep].add(chain_mem)
           
            print('2: membership read done', flush=True)

            # write seqs cluster index
            with open(output.cluster_index, 'w') as f:
                for cluster_name in sorted(list(membership.keys())):
                    f.write(f'{cluster_name}\n')
            
            print('3: write cluster index done', flush=True)

            i = 0
            I = len(membership.keys())
            
            # write dimers and chains membership files for each cluster
            for cluster_name, member_chains in membership.items():
                template = outfile['cluster']['cluster_members_chains_template']
                chains_file = template.format(cluster_name=cluster_name)
                template = outfile['cluster']['cluster_members_dimers_template'] 
                dimers_file = template.format(cluster_name=cluster_name)

                # write member chains of one cluster
                os.makedirs(os.path.dirname(chains_file), exist_ok=True)
                with open(chains_file, 'w') as f:
                    for chain in sorted(member_chains):
                        f.write(f'{chain}\n')

                # derive member dimers
                dimer_in_clust = lambda dimer: name_pdb.dimer2chains(dimer)[0] in member_chains
                member_dimers = filter(dimer_in_clust, homodimers)
                
                # write member dimers of one cluster
                with open(dimers_file, 'w') as f: 
                    for dimer in member_dimers:
                        f.write(f'{dimer}\n')

                i += 1
                if i % 100 == 0:
                    print(f'4: {i}/{I} steps done', flush=True)



