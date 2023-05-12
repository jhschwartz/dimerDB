import os
import sys

configfile: 'config.yaml'

if config.get('test'):
   outdata = 'tmpout'
else:
    outdata = config['paths']['out_dir']


from datetime import datetime


sys.path.append('scripts')
import name_pdb
from tar_pdb_chains import tar_pdb_chains


subworkflow_done = config['workflow']['4_prep_webserver_data_done']

# data locations set by config
intermediates = config['paths']['intermediates_dir']
lib_path = config['paths']['lib'] # TODO put in config, match to sub0, must contain folder "rcsb"

# input data, given data locations
pdb_index = os.path.join(lib_path, 'pdb_index.txt') # TODO match sub0
fasta_index = os.path.join(lib_path, 'fasta_index.tsv')
all_homodimers_file = os.path.join(intermediates, 'check_pairs', 'all_homodimers.txt')

cluster_index = os.path.join(intermediates, 'cluster', 'cluster_index')
cluster_data = os.path.join(intermediates, 'cluster', 'cluster_index', 'seq_clusters')



# outfile naming
outfile = {
    'pdb': {
        'seqs': os.path.join(outdata, 'pdb', 'seqs.fasta'),
        'tar_template': outdata+'/pdb/div/{div}.tar.gz'
    },
    'all': {
        'chains': os.path.join(outdata, 'all', 'chains.txt'), 
        'dimers': os.path.join(outdata, 'all', 'dimers.txt'),
        'seqs_shorter': os.path.join(outdata, 'all', 'seqs_shorter.fasta'),
        'seqs_longer': os.path.join(outdata, 'all', 'seqs_longer.fasta')
    },
    'cluster': {
        'clusters': os.path.join(outdata, 'cluster', 'clusters.txt'),
        'membership': os.path.join(outdata, 'cluster', 'membership.tsv')
    },
    'nonredundant': {
        'chains': os.path.join(outdata, 'nonredundant', 'chains.txt'), 
        'dimers': os.path.join(outdata, 'nonredundant', 'dimers.txt'),
        'seqs_shorter': os.path.join(outdata, 'nonredundant', 'seqs_shorter.fasta'),
        'seqs_longer': os.path.join(outdata, 'nonredundant', 'seqs_longer.fasta')
    },
    'extra': {
        'dimers_template': outdata+'/extra/cluster{num}_dimers.txt',
        'seqs_longer_template': outdata+'/extra/cluster{num}_seqs_longer.fasta',
        'seqs_shorter_template': outdata+'/extra/cluster{num}_seqs_shorter.fasta'
    },
    'summary': os.path.join(outdata, 'summary.txt') 
}









### retrieve list of divs used in pdb_index  ###
### and use for rule tar_pdb_chain_lib later ###
chain_divs = set()
with open(pdb_index, 'r') as f:
    for line in f:
        chain_name = line.strip()
        div = name_pdb.get_div(chain_name)
        chain_divs.add(div)


# set ids for running samples of rule out_extra_clusters
ids_for_extra = [30, 40, 50, 60, 70, 80]



# rule all and the two rules that use wildcard samples, 
# out_pdb_div_tar, and out_extra_clusters are defined 
# first because their samples are defined just above
rule all:
    input:
        # rule out_pdb_div_tar
        expand(outfile['pdb']['tar_template'], div=chain_divs),

#        # rule out_extra_clusters
#        expand(outfile['extra']['dimers_template'], num=ids_for_extra), 
#        expand(outfile['extra']['seqs_shorter_template'], num=ids_for_extra), 
#        expand(outfile['extra']['seqs_longer_template'], num=ids_for_extra)
    




rule out_pdb_div_tar:
    '''
    This rule archives all single chain pdb files that 
    are found to be involved in homodimerization. That
    is, chains of all pairs found to be of the same
    sequence and to be in contact, before any redundancy
    removal. PDB chains are divided according to the
    middle two characters of their entry ID, and each
    division is saved in its own archive, e.g. 1aa1.tar.gz
    These are written to output files at pdb/div

    This rule utilizes samples, where each sample is a
    division of pdb entries, to parallelize this process.
    '''
    output:
        tarf = outfile['pdb']['tar_template']
    threads: 1
    resources:
        time = '2:00:00',
        memory = '5000'
    run:
        tar_pdb_chains(indexfile=pdb_index,
                       div=wildcards.div,
                       chains_lib=lib_path,
                       target_dir=os.path.dirname(output.tarf)
                       )




rule out_extra_clusters:
    '''
    Given the resulting dimers and seqs written to nonredundant/,
    this rule further clusters the nonredundant dimers by 
    sequence identity according to the shorter chains' seqs.
    This clustering is done to a variety of sequence identities
    30%-80% to create sub datasets for non-ML purposes where a 
    user might want a nonredundant, but not comprehensive, 
    dimer dataset according to some sequence identity.

    This rule utilizes samples defined above in "ids_for_extra"
    to parallelize the rule across each id 30%-80%.
    
    Uses all/dimers.txt to preserve length-based naming.
    '''
    input:
        dimers_txt = outfile['nonredundant']['dimers'],
        seqs_shorter = outfile['nonredundant']['seqs_shorter'],
        seqs_longer = outfile['nonredundant']['seqs_longer']
    output:
        dimers_txt = outfile['extra']['dimers_template'], 
        seqs_shorter = outfile['extra']['seqs_shorter_template'],
        seqs_longer = outfile['extra']['seqs_longer_template']
    run:
        raise NotImplementedError








# defacto input is the fasta library at lib_path and fasta_index
rule out_pdb_one_fasta:
    '''
    This rule writes the sequences of all chains involved
    in homodimerization to a single output file pdb/seqs.fasta
    These sequences are simply collected from lib fasta_index.txt
    and the fasta files in <lib_path>/fasta/...
    '''
    output:
        fasta = outfile['pdb']['seqs']
    run:
        raise NotImplementedError





# defacto input is all_homodimers_file 
rule out_all_chains_dimers:
    '''
    This rule writes all homodimer pairs discovered, without 
    consideration for redundancy removal, to output file
    all/dimers.txt, as well as all the chains which make up
    these dimers to all/chains.txt. Both output files are 
    sorted, and dimers.txt is named using the shorter chain
    name first. This lengths are collected from fasta_index.
    
    The dimers are collected from subworkflow2's intermediate
    file that lists all homodimers, all_homodimers_file
    '''
    output:
        chains_txt = outfile['all']['chains'],
        dimers_txt = outfile['all']['dimers']
    run:
        raise NotImplementedError



# defacto inputs data in lib/fasta
rule out_all_seqs:
    '''
    This rule writes the sequences corresponding to all
    dimers discovered, before redundancy removal. It 
    writes both all/seqs_shorter.fasta and 
    all/seqs_longer.fasta, corresponding respectively to
    the sequences which make up every dimer's shorter
    chain and longer chain.

    Uses all/dimers.txt to preserve length-based naming.
    '''
    input:
        dimers_txt = outfile['all']['dimers']
    output:
        fasta_shorter = outfile['all']['seqs_shorter'],
        fasta_longer  = outfile['all']['seqs_longer']
    run:
        raise NotImplementedError



# defacto inputs cluster_index and many representative.tsv files
rule out_cluster_info:
    '''
    This rule collects and writes sequence clustering and
    structural grouping information. Specifically, uses
    cluster_index and many representative.tsv files across
    cluster data to create a consolidated cluster/cluster.txt
    and cluster/membership.tsv

    cluster.txt gives the names of all sequence clusters, 
    given by the name of the dimer whose shorter chain
    is the center of the cluster.

    membership.tsv defines the membership of every dimer
    (in all/dimers.txt) to a sequence cluster and to a
    structural group, where the sequence cluster is 
    defined the same way as above, and the structure 
    cluster is defined by the name of the dimer that was
    chosen as the representative structure of the group.

    Uses all/dimers.txt to preserve length-based naming.
    '''
    input:
        dimers_txt = outfile['all']['dimers']
    output:
        clusters_txt = outfile['cluster']['clusters'],
        membership_tsv = outfile['cluster']['membership']
    run:
        raise NotImplementedError




rule out_nonredundant_dimers_chains:
    '''
    This rule writes nonredundant homodimer pairs to 
    output file nonredundant/dimers.txt, as well as all the 
    chains which make up these dimers to 
    nonredundant/chains.txt. Both output files are sorted.
   
    This data is collected from the third column of
    cluster/membership.tsv, generated by the rule above.
    '''
    input:
        membership_tsv = outfile['cluster']['membership']
    output:
        dimers_txt = outfile['nonredundant']['dimers'], 
        chains_txt = outfile['nonredundant']['chains']
    run:
        raise NotImplementedError





# defacto input sequences in lib/fasta
rule out_nonredundant_seqs:
    '''
    This rule writes the sequences corresponding to nonredundant
    dimers. It writes both nonredundant/seqs_shorter.fasta and 
    nonredundant/seqs_longer.fasta, corresponding respectively to
    the sequences which make up each nonredndant dimer's shorter
    chain and longer chain.
    '''
    input:
        dimers_txt = outfile['nonredundant']['dimers']
    output:
        seqs_shorter_fasta = outfile['nonredundant']['seqs_shorter'],
        seqs_longer_fasta = outfile['nonredundant']['seqs_longer']
    run:
        raise NotImplementedError










rule out_summary:
    '''
    This rule writes statistics on the output folder
    and timestamps it.
    '''
    input:
        all_chains = outfile['all']['chains'],
        all_dimers = outfile['all']['dimers'], 
        nonred_chains = outfile['nonredundant']['chains'],
        nonred_dimers = outfile['nonredundant']['dimers']
    run:
        now = str(datetime.utcnow()) 

        with open(input.all_chains, 'r') as f:
            all_num_chains = len(f.readlines())
        with open(input.all_dimers, 'r') as f:
            all_num_dimers = len(f.readlines())
        with open(input.nonred_chains, 'r') as f:
            nr_num_chains = len(f.readlines())
        with open(input.nonred_dimers, 'r') as f:
            nr_num_dimers = len(f.readlines())

        with open(outfile['summary'], 'w') as f:
            f.write(f'LAST_UPDATE_UTC\t{now}\n')
            f.write(f'FULL_SET_NUM_CHAIN\t{all_num_chains}\n')
            f.write(f'FULL_SET_NUM_DIMER\t{all_num_dimers}\n')
            f.write(f'NONRED_NUM_CHAIN\t{nr_num_chains}\n')
            f.write(f'NONRED_NUM_DIMER\t{nr_num_dimers}\n')




