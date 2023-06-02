import os
import sys
import tempfile
from sortedcontainers import SortedSet
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

if config.get('test'):
    all_homodimers = os.path.join(lib_path, 'test_in_homodimers.txt')
    
max_threads = config['workflow_params']['max_threads']
seq_cluster_id = config['workflow_params']['seq_cluster_id']

mmseqs_exe = config['exe']['mmseqs']

outfile = {
    'fasta': {
        'index': os.path.join(lib_path, 'fasta_index.tsv'),
        'homodimers_seqs': os.path.join(intermediates, 'homodimer_seqs.fasta'),
        'all_seqs': os.path.join(lib_path, 'all_seqs.fasta'),
        'single_seq_template': os.path.join(lib_path, 'fasta', '{div}', '{chain}.fasta')
    },
    'cluster': {
        'index': os.path.join(intermediates, 'cluster', 'cluster_index.txt'),
        'cluster_members_dimers_template': os.path.join(intermediates, 'cluster', 
                                            'seq_clusters', '{cluster_name}', 'members_dimers.txt'),
        'cluster_members_chains_template': os.path.join(intermediates, 'cluster', 
                                            'seq_clusters', '{cluster_name}', 'members_chains.txt')
    }
}




localrules: all
rule all:
#### BEGIN RULE ALL TARGETS
    input:
        indexfile = outfile['cluster']['index']
    output:
        done = subworkflow_done
    run:
        with open(output.done, 'w') as f:
            f.write('Subworkflow 2 done at ')
            f.write(str(datetime.utcnow()))
            f.write('\n')
#### END RULE ALL TARGETS
# defining begin/end of rule all in comments is necessary for unittesting






# shadow input pdb_index
rule process_fastas:
    '''
    This rule, which runs only once without samples, runs
    pdb2fasta to generate one fasta for the entire set of
    pdb chains found to homodimerize by subflow1.
    This large fasta then has its headers reformatted,
    then individual fasta files for each sequence are 
    written to {lib_path}/fasta. Additionally, an index
    of all these fastas is written to {fasta_index}.
    '''
    output:
        all_seqs_fasta = outfile['fasta']['all_seqs'],
        indexfile_tsv = outfile['fasta']['index']
    run:
        rcsb_dir = os.path.join(lib_path, 'rcsb/')
        
        # write all pdbs to one fasta and edit its headers
        shell(
            ''' 
            set +e;
            bin/USalign/pdb2fasta -mol prot -het 2 -dir {rcsb_dir} {pdb_index} > {output.all_seqs_fasta}; 
            sed -i -r "s/>.*\//>/g" {output.all_seqs_fasta};
            sed -i -r "s/\.pdb\S*//g" {output.all_seqs_fasta};
            '''
        )  # resulting headers like ">ql/6qle-a1-m1-c1.fasta 123", where 123 is seq length

        # write each seq to its own fasta file and keep index
        with open(output.indexfile_tsv, 'w') as indexf:
            for header, seq in read_fasta.read_prot_from_fasta(output.all_seqs_fasta):
                chain = header.lstrip('>').split()[0]
                L = header.rstrip().split()[1]
                div = name_pdb.get_div(chain)
                outfasta = outfile['fasta']['single_seq_template'].format(div=div, chain=chain) 
                os.makedirs(os.path.dirname(outfasta), exist_ok=True)
                with open(outfasta, 'w') as fo:
                    fo.write(f'{header}\n{seq}\n')
                indexf.write(f'{div}/{chain}.fasta\t{L}\n')


        # sort the index file - the "-s", "-t" and "-k" sort only by letters before ".pdb"
        shell(''' sort -s -t "." -k 1,1 -o {output.indexfile_tsv} {output.indexfile_tsv} ''')



# shadow input is {all_homodimers}
rule write_homodimers_chains_fastas:
    input:
        all_seqs_fasta = outfile['fasta']['all_seqs']
    output:
        fasta = outfile['fasta']['homodimers_seqs']
    run:
        with open(all_homodimers, 'r') as f:
            homodimers = [line.strip() for line in f]
        
        homodimer_chains = SortedSet()
        for dimer in homodimers:
            c1, _ = name_pdb.dimer2chains(dimer)
            homodimer_chains.add(c1)

        with open(output.fasta, 'w') as f:
            for header, seq in read_fasta.read_prot_from_fasta(input.all_seqs_fasta):
                chain = header.split()[0].lstrip('>')
                if chain in homodimer_chains:
                    f.write(f'{header}\n{seq}\n')




# shadow input is {all_homodimers} 
rule run_seq_cluster:
    input:
        fasta = outfile['fasta']['homodimers_seqs'],
    output:
        cluster_index = outfile['cluster']['index']
    threads: max_threads
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

            # read membership relationships from mmseqs output
            membership = {}
            with open(mmseqs_cluster_tsv, 'r') as f:
                for line in f:
                    chain_rep, chain_mem = line.rstrip().split()
                    if chain_rep in membership:
                        membership[chain_rep].append(chain_mem)
                    else:
                        membership[chain_rep] = [chain_mem]
            
            # write seqs cluster index
            with open(output.cluster_index, 'w') as f:
                for cluster_name in sorted(list(membership.keys())):
                    f.write(f'{cluster_name}\n')
            
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





