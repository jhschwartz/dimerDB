import os
import sys
from datetime import datetime

configfile: 'config.yaml'

if config.get('test'):
   outdata = 'tmpout'
else:
    outdata = config['paths']['out_dir']


from datetime import datetime
import tempfile
from sortedcontainers import SortedSet
import gzip
import shutil

sys.path.append('scripts')
import name_pdb
from tar_pdb_chains import tar_pdb_chains
from wrap_mmseqs import run_mmseqs_cluster, derive_mmseqs_reps_from_tsv
from align_tools import calc_nwalign_glocal

sys.path.append('bin/check_contact')
from check_contact import count_contact_many_parallel

subworkflow_done = config['subworkflow_done']['4_prep_final_data']

# data locations set by config
intermediates = config['paths']['intermediates_dir']
lib_path = config['paths']['lib_full'] 

# input data, given data locations
pdb_index = os.path.join(lib_path, 'pdb_index.txt') 
fasta_index = os.path.join(lib_path, 'fasta_index.tsv')
fasta_lib = os.path.join(lib_path, 'fasta')
all_homodimers_file = os.path.join(intermediates, 'check_pairs', 'all_homodimers.txt')

cluster_index = os.path.join(intermediates, 'cluster', 'cluster_index.txt')
cluster_data = os.path.join(intermediates, 'cluster', 'seq_clusters')

mmseqs_exe = config['exe']['mmseqs']


# outfile naming
outfile = {
    'pdb': {
        'tar_template': outdata+'/pdb/div/{div}.tar.gz'
    },
    'all': {
        'chains': os.path.join(outdata, 'all', 'chains.txt'), 
        'dimers': os.path.join(outdata, 'all', 'dimers.txt'),
        'seqs': os.path.join(outdata, 'all', 'seqs.fasta')
    },
    'cluster': {
        'clusters': os.path.join(outdata, 'cluster', 'clusters.txt'),
        'membership': os.path.join(outdata, 'cluster', 'membership.tsv')
    },
    'nonredundant': {
        'chains': os.path.join(outdata, 'nonredundant', 'chains.txt'), 
        'dimers': os.path.join(outdata, 'nonredundant', 'dimers.txt'),
        'seqs': os.path.join(outdata, 'nonredundant', 'seqs.fasta'),
        'info': os.path.join(outdata, 'nonredundant', 'dimers_info.tsv')
    },
    'extra': {
        'dimers_template': outdata+'/extra/cluster{num}_dimers.txt',
        'seqs_template': outdata+'/extra/cluster{num}_seqs.fasta',
        'chains_template': outdata+'/extra/cluster{num}_chains.txt'
    },
    'summary': os.path.join(outdata, 'summary.txt'),
    'info': {
        'seqids': os.path.join(outdata, 'info', 'seq_ids.tsv'),
        'contact_counts': os.path.join(outdata, 'info', 'contact_counts.tsv'),
        'entries': os.path.join(outdata, 'info', 'entries.idx'),
        'pdb_chain_uniprot': os.path.join(outdata, 'info', 'pdb_chain_uniprot.tsv'),
        'pdb_chain_taxonomy': os.path.join(outdata, 'info', 'pdb_chain_taxonomy.tsv'),
        'ncbi_taxonomy_archive': os.path.join(outdata, 'info', 'taxdump.tar.gz'),
        'ncbi_taxonomy_names': os.path.join(outdata, 'info', 'names.dmp'),
        'taxid2name': os.path.join(outdata, 'info', 'names_processed.tsv')
    }
}




## read and save all chain lengths into memory ##
chain_lengths = {}
with open(fasta_index, 'r') as f:
    
    # fasta_index line is like "ql/6qle-a1-m1-cH.fasta  126"
    for line in f:
        chain_fasta = line.split()[0]
        chain_name = os.path.basename(chain_fasta).split('.fasta')[0]
        chain_length = int(line.split()[1])
        chain_lengths[chain_name] = chain_length


## func name_dimer_by_len uses chain_lengths to name
## every dimer shorter-chain_longer-chain
def name_dimer_by_len(dimer_name):
    c1, c2 = name_pdb.dimer2chains(dimer_name)
    if chain_lengths[c2] < chain_lengths[c1]:
        dimer_name = f'{c2}_{c1}'
    return dimer_name



## func write_fastas_of_chains writes the seqs
## that correspond to a list of chain names to
## one fasta file; seqs come from the specified
## fasta_lib. This is reused a lot.
def write_fastas_of_chains(chains_list, fasta_lib, outfasta):
    fasta_files = []
    for chain in sorted(list(chains_list)): 
        div = name_pdb.get_div(chain)
        fasta_path = os.path.join(fasta_lib, div, f'{chain}.fasta')
        fasta_files.append(fasta_path)
    with open(outfasta, 'w') as fo:
        for chain_fasta in fasta_files:
            with open(chain_fasta, 'r') as fi:
                fo.write(fi.read())
    




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
localrules: all
rule all:
    input:
#### BEGIN RULE ALL TARGETS
        # target of out_pdb_div_tar
        expand(outfile['pdb']['tar_template'], div=chain_divs),

        ## targets of out_extra_clusters
        expand(outfile['extra']['dimers_template'], num=ids_for_extra), 
        expand(outfile['extra']['seqs_template'], num=ids_for_extra),
        expand(outfile['extra']['chains_template'], num=ids_for_extra),

	# target of everything
	outfile['nonredundant']['info']
    output:
        done = config['subworkflow_done']['4_prep_final_data']
    run:
        # cleanup
        shutil.rmtree(os.path.dirname(outfile['info']['seqids']), ignore_errors=True)
        with open(output.done, 'w') as f:
            f.write('Subworkflow 4 done at ')
            f.write(str(datetime.utcnow()))
            f.write('\n')
#### END RULE ALL TARGETS
# defining begin/end of rule all in comments is necessary for unittesting



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
        tar_pdb_chains(  indexfile=pdb_index,
                         div=wildcards.div,
                         chains_lib=lib_path,
                         target_dir=os.path.dirname(output.tarf)  )




# shadow input is the fasta library at lib_path and fasta_index
rule out_extra_clusters:
    '''
    Given the resulting dimers and seqs written to nonredundant/,
    this rule further clusters the nonredundant dimers by 
    sequence identity according to the shorter chains' seqs.
    This clustering is done to a variety of sequence identities
    30%-80% to create sub datasets for non-ML purposes where a 
    user might want a nonredundant, but not comprehensive, 
    dimer dataset according to some sequence identity.

    Note that output file extra/seqs<num>_seqs.fasta has the seqs
    of all dimers (both chains' seqs) in the named subset, while 
    extra/seqs<num>_dimers.txt names the dimes in the subset.

    This rule utilizes samples defined above in "ids_for_extra"
    to parallelize the rule across each id 30%-80%.
    '''
    input:
        dimers_txt = outfile['nonredundant']['dimers'],
    output:
        dimers_txt = outfile['extra']['dimers_template'], 
        chains_txt = outfile['extra']['chains_template'],
        seqs_fasta = outfile['extra']['seqs_template'] 
    threads: 8
    run:
        # read nonred dimers
        dimers = []
        with open(input.dimers_txt, 'r') as f:
            for line in f:
                dimers.append(line.strip())

        # get the shorter chains of nonred dimers
        shorter_chains = set()        
        for d in dimers:
            shorter_chain, _ = name_pdb.dimer2chains(d)
            shorter_chains.add(shorter_chain)

        # do seq clustering in a tempdir
        with tempfile.TemporaryDirectory() as td:
            tmp_inseqs = os.path.join(td, 'in.fasta')
            tmp_outprefix = os.path.join(td, 'out')
            
            # write temp fasta of shorter_chains' seqs for mmseqs input
            write_fastas_of_chains(  chains_list=sorted(list(shorter_chains)), 
                                     fasta_lib=fasta_lib,
                                     outfasta=tmp_inseqs                       )

            # cluster the shorter_chains' sequences; outputs put in td/out*
            mmseqs_membership_tsv, _ = \
                run_mmseqs_cluster(  infasta=tmp_inseqs, 
                                     outprefix=tmp_outprefix, 
                                     seq_id=wildcards.num,
                                     mmseqs_exe=mmseqs_exe,
                                     cores=threads           )
            mmseqs_reps = derive_mmseqs_reps_from_tsv(mmseqs_membership_tsv) 


        # derive the resulting dimers
        resulting_dimers = []
        for dimer in dimers:
            shorter_chain, _ = name_pdb.dimer2chains(dimer)
            if shorter_chain in mmseqs_reps:
                resulting_dimers.append(dimer)

        # write resulting dimers to file
        with open(output.dimers_txt, 'w') as f:
            for dimer in sorted(resulting_dimers):
                f.write(f'{dimer}\n')

        # use resulting dimers to resolve all chains used in this subset
        resulting_chains = set()
        for dimer in resulting_dimers:
            c1, c2 = name_pdb.dimer2chains(dimer)
            resulting_chains.add(c1)
            resulting_chains.add(c2)

        # write resulting chains to file
        with open(output.chains_txt, 'w') as f:
            for chain in sorted(list(resulting_chains)):
                f.write(f'{chain}\n')

        # write seqs of resulting chains to file
        write_fastas_of_chains(  chains_list=sorted(list(resulting_chains)), 
                                 fasta_lib=fasta_lib,
                                 outfasta=output.seqs_fasta                  )
        
             
# CONSIDERING REMOVAL
#
## shadow input is the fasta library at lib_path and fasta_index
#rule out_all_seqs:
#    '''
#    This rule writes the sequences of all chains involved
#    in homodimerization to a single output file all/seqs.fasta
#    These sequences are simply collected from lib fasta_index.txt
#    and the fasta files in <lib_path>/fasta/...
#    '''
#    output:
#        fasta = outfile['all']['seqs']
#    run:
#        chains_list = set() 
#        with open(fasta_index, 'r') as f:
#            for line in f:
#                path = line.split()[0]
#                chain_name = os.path.basename(path).split('.fasta')[0]
#                chains_list.add(chain_name)
#        write_fastas_of_chains(  chains_list=sorted(list(chains_list)), 
#                                 fasta_lib=fasta_lib, 
#                                 outfasta=output.fasta                  )
#
#



# shadow inputs are all_homodimers_file and chain_lengths
rule out_all_dimers_chains:
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
        dimers = []
        with open(all_homodimers_file, 'r') as f:
            for line in f:
                dimer = line.strip()
                dimer = name_dimer_by_len(dimer)
                dimers.append(dimer)
        
        dimers.sort()
        with open(output.dimers_txt, 'w') as f:
            for dimer in dimers:
                f.write(f'{dimer}\n')

        chains = set()
        for dimer in dimers:
            c1, c2 = name_pdb.dimer2chains(dimer)
            chains.add(c1)
            chains.add(c2)

        chains = sorted(list(chains))
        with open(output.chains_txt, 'w') as f:
            for chain in chains:
                f.write(f'{chain}\n')



# shadow inputs cluster_index and each cluster's representation.tsv
rule out_cluster_info:
    '''
    This rule collects and writes sequence clustering and
    structural grouping information. Specifically, uses
    cluster_index and many representation.tsv files across
    cluster data to create a consolidated cluster/cluster.txt
    and cluster/membership.tsv

    cluster.txt gives the names of all sequence clusters, 
    given by the name of the chain whose sequence is the
    center of the cluster.

    membership.tsv defines the membership of every dimer
    (in all/dimers.txt) to a sequence cluster and to a
    structural group, where the sequence cluster is 
    defined the same way as above, and the structure 
    cluster is defined by the name of the dimer that was
    chosen as the representative structure of the group.
    '''
    output:
        clusters_txt = outfile['cluster']['clusters'],
        membership_tsv = outfile['cluster']['membership']
    run:
        # 1: create paths to all clusters' representation files
        cluster_repfiles = []
        with open(cluster_index, 'r') as f:
            for line in f:
                cluster_name = line.strip()
                repfile = os.path.join(cluster_data, 
                                       cluster_name, 
                                       'representation.tsv'
                                       )
                cluster_repfiles.append( (cluster_name, repfile) )

        # 2: read relationships between sequence clusters and
        #    structural groups from each representation file
        membership = {}
        for cluster_name, repfile in cluster_repfiles:
            with open(repfile, 'r') as f:
                for line in f:
                    member_dimer, rep_dimer = line.split()
                    member_dimer = name_dimer_by_len(member_dimer)
                    rep_dimer = name_dimer_by_len(rep_dimer) 
                    membership[member_dimer] = {
                        'cluster': cluster_name,
                        'rep': rep_dimer
                    }

        # 3: write list of all clusters and sort
        with open(output.clusters_txt, 'w') as f:
            for cluster, _ in cluster_repfiles:
                f.write(f'{cluster}\n')
        shell(''' sort -o {output.clusters_txt} {output.clusters_txt} ''') 

        # 4: write membership table and sort
        with open(output.membership_tsv, 'w') as f:
            for member_dimer in membership.keys():
                cluster = membership[member_dimer]['cluster']
                rep = membership[member_dimer]['rep']
                f.write(f'{member_dimer}\t{cluster}\t{rep}\n')
        shell(''' sort -o {output.membership_tsv} {output.membership_tsv} ''') 






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
        dimers = set()
        with open(input.membership_tsv, 'r') as f:
            for line in f:
                dimers.add(line.split()[2])

        with open(output.dimers_txt, 'w') as f:
            for d in sorted(list(dimers)):
                f.write(f'{d}\n')
        
        chains = set()
        for d in dimers:
            c1, c2 = name_pdb.dimer2chains(d)
            chains.add(c1)
            chains.add(c2)

        with open(output.chains_txt, 'w') as f:
            for c in sorted(list(chains)):
                f.write(f'{c}\n')






# shadow input sequences in lib/fasta
rule out_nonredundant_seqs:
    '''
    This rule writes the sequences corresponding to nonredundant
    dimers. It writes nonredundant/seqs.fasta, corresponding to
    the sequences of chains that make up nonredundat/chains.txt
    '''
    input:
        chains_txt = outfile['nonredundant']['chains']
    output:
        fasta = outfile['nonredundant']['seqs']
    run:
        chains_list = set()
        with open(input.chains_txt, 'r') as f:
            for line in f:
                chains_list.add(line.strip())
        write_fastas_of_chains(  chains_list=sorted(list(chains_list)),
                                 fasta_lib=fasta_lib, 
                                 outfasta=output.fasta                  )






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






rule download_info_flatfiles:
    output:
       entries_file = outfile['info']['entries'],
       pdb2uniprot_tsv = outfile['info']['pdb_chain_uniprot'],
       pdb2species_tsv = outfile['info']['pdb_chain_taxonomy'],
       ncbi_taxonomy_names = outfile['info']['ncbi_taxonomy_names'],
       ncbi_taxonomy_targz = outfile['info']['ncbi_taxonomy_archive']
    threads: 1
    shell:
        '''
        wget -O {output.entries_file} https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx;
        wget -O {output.pdb2uniprot_tsv} ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz;
        wget -O {output.pdb2species_tsv} ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_taxonomy.tsv.gz;
        wget -O {output.ncbi_taxonomy_targz} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz;
        tar -zxvf {output.ncbi_taxonomy_targz} names.dmp --to-stdout > {output.ncbi_taxonomy_names}; 
        '''



rule process_ncbi_taxonomy_file:
    input:
        ncbi_taxonomy_names_raw = outfile['info']['ncbi_taxonomy_names']
    output:
        ncbi_taxonomy_names_processed = outfile['info']['taxid2name']
    run:
        # from ncbi, for input names.dmp file:
        #   Field terminator is "\t|\t"
        #   Row terminator is "\t|\n"

        with open(input.ncbi_taxonomy_names_raw, 'r') as fi, \
                    open(output.ncbi_taxonomy_names_processed, 'w') as fo:
           
            is_sci_name = lambda row: 'scientific name' in row.lower()

            for row in filter(is_sci_name, fi):
                tax_id, sci_name, _, _ = row.rstrip('\t|\n').split('\t|\t')
                fo.write(f'{tax_id}\t{sci_name}\n')



rule bulk_count_contacts:
    input:
        nonred_dimers = outfile['nonredundant']['dimers']
    output:
        counts_file = outfile['info']['contact_counts']
    threads: 8
    run:
        dimer_tuples = []
        dimer_names = []
        with open(input.nonred_dimers, 'r') as f:
            for line in f:
                dimer = line.strip()
                dimer_names.append(dimer)
                c1, c2 = name_pdb.dimer2pdbs(dimer_name=dimer, lib_path=lib_path)
                dimer_tuples.append( (c1, c2) )

        counts = count_contact_many_parallel(  pairs_list=dimer_tuples,
                                               thresh_max_dist=8,
                                               cores=threads             )

        with open(output.counts_file, 'w') as f:
            for count, dimer in zip(counts, dimer_names):
                f.write(f'{dimer}\t{count}\n')




rule bulk_calc_seqid:
    input:
        nonred_dimers = outfile['nonredundant']['dimers']
    output:
        seqids_file = outfile['info']['seqids']
    run:
        with open(input.nonred_dimers, 'r') as fi, \
                open(output.seqids_file, 'w') as fo:
            for line in fi:
                dimer = line.strip()
                pdb1, pdb2 = name_pdb.dimer2pdbs(dimer_name=dimer, lib_path=lib_path)
                score = calc_nwalign_glocal('bin/USalign/NWalign', pdb1, pdb2)
                fo.write(f'{dimer}\t{score}\n')




rule out_dimers_info:
    input:
        nonred_dimers = outfile['nonredundant']['dimers'],
        contact_counts_file = outfile['info']['contact_counts'],
        seqids_file = outfile['info']['seqids'],
        entries_file = outfile['info']['entries'],
        pdb2uniprot_tsv = outfile['info']['pdb_chain_uniprot'],
        pdb2species_tsv = outfile['info']['pdb_chain_taxonomy'],
        taxid2name_tsv = outfile['info']['taxid2name']
    output:
        infofile = outfile['nonredundant']['info']
    run:
        dimers = []
        chains_nonred = SortedSet()
        entries_nonred = SortedSet()
        with open(input.nonred_dimers, 'r') as f:
            for line in f:
                dimer = line.strip()
                dimers.append(dimer)
                c1, c2 = name_pdb.dimer2chains(dimer)
                entry, _, _ , c1 = name_pdb.read_chain_names(c1)
                _, _, _, c2 = name_pdb.read_chain_names(c2)
                c1, c2 = f'{entry}_{c1}', f'{entry}_{c2}'
                chains_nonred.add(c1)
                chains_nonred.add(c2)
                entries_nonred.add(entry)

        with open(input.contact_counts_file, 'r') as f:
            contact_counts = [line.split()[1] for line in f]

        with open(input.seqids_file, 'r') as f:
            seqids = [line.split()[1] for line in f]

        chain2uniprot = {}
        with gzip.open(input.pdb2uniprot_tsv, 'rt') as f:
            f.readline() # skip line 1 - timestamp
            f.readline() # skip line 2 - header
            for line in f:
                entry, chain, uniprot = line.split('\t')[:3]
                entry = entry.lower()
                chain_name = f'{entry}_{chain}'
                if chain_name in chains_nonred:
                    chain2uniprot[chain_name] = uniprot
        
        chain2taxid = {}
        taxids_of_interest = SortedSet()
        with gzip.open(input.pdb2species_tsv, 'rt') as f:
            f.readline() # skip line 1 - timestamp
            f.readline() # skip line 2 - header
            for line in f:
                entry, chain, species_id, species_name = line.rstrip().split('\t')
                chain_name = f'{entry}_{chain}'
                if chain_name in chains_nonred:
                    chain2taxid[chain_name] = species_id
                    taxids_of_interest.add(species_id)

        taxid2species = {}
        with open(input.taxid2name_tsv, 'r') as f:
            for line in f:
                taxid, taxname = line.rstrip().split('\t')
                if taxid in taxids_of_interest:
                    taxid2species[taxid] = taxname

        entry2title = {}
        entry2res = {}
        entry2method = {}
        entry2date = {}
        with open(input.entries_file, 'r') as f:
            for line in f:
                spl = line.split('\t')
                entry = spl[0].lower()
                if entry in entries_nonred:
                    entry2date[entry] = spl[2]
                    entry2title[entry] = spl[3]
                    entry2res[entry] = spl[6]
                    entry2method[entry] = spl[7].rstrip()
                    if spl[7] == '':
                        entry2method[entry] = 'X-RAY DIFFRACTION'
        
        with open(output.infofile, 'w') as f:
            for dimer, num_contacts, seqid in zip(dimers, contact_counts, seqids):
                c1, c2 = name_pdb.dimer2chains(dimer)
                entry, _, _ , c1 = name_pdb.read_chain_names(c1)
                _, _, _, c2 = name_pdb.read_chain_names(c2)
                c1, c2 = f'{entry}_{c1}', f'{entry}_{c2}'
                
                title = entry2title[entry]
                up1 = chain2uniprot.get(c1, '')
                up2 = chain2uniprot.get(c2, '')
                res = entry2res[entry]
                method = entry2method[entry]
                
                c1_taxid = chain2taxid.get(c1, '')
                c2_taxid = chain2taxid.get(c2, '')
                c1_taxname = taxid2species.get(c1_taxid, '')
                c2_taxname = taxid2species.get(c2_taxid, '')
                c1_species = c1_taxid
                c2_species = c2_taxid
                if c1_species != '':
                    c1_species = f'{c1_taxid} ({c1_taxname})'
                if c2_species != '':
                    c2_species = f'{c2_taxid} ({c2_taxname})'

                f.write(f'{dimer}\t{title}\t{up1}\t{up2}\t{res}\t{method}')
                f.write(f'\t{num_contacts}\t{seqid}\t{c1_species}\t{c2_species}\n')

