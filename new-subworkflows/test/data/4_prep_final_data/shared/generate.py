import pickle
import itertools
import random
import tempfile
import subprocess
import sys
import os 
import shutil
from datetime import datetime


outdir = './expected'
#pdbdir = '/nfs/turbo/umms-petefred/jaschwa/explore_homodimers/truncate-lib/rcsb'
pdbdir = '/nfs/turbo/umms-petefred/jaschwa/dimerDB/lib/rcsb'
index = '/nfs/turbo/umms-petefred/jaschwa/explore_homodimers/truncate-lib/rcsb_index.pkl'
check_contact_dir = '/nfs/turbo/umms-petefred/jaschwa/check-dimer-contact'
pdb2fasta = '/nfs/turbo/umms-petefred/jaschwa/dimerDB-optimize-dev/bin/USalign/pdb2fasta'
read_fasta_dir = '/nfs/turbo/umms-petefred/jaschwa/dimerDB-optimize-dev/scripts'


sys.path.append(check_contact_dir)
from check_contact import check_contact_many

sys.path.append(read_fasta_dir)
from read_fasta import read_prot_from_fasta as read_fasta


def pdb_path(pdb):
    div = pdb[1:3]
    path = f'{pdbdir}/{div}/{pdb}'
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} not existing')
    return path 

def dimer_pdb_path(dimer_tuple):
    return (pdb_path(dimer_tuple[0]), pdb_path(dimer_tuple[1]))

def seq_of_pdb(pdbfile):
    with tempfile.NamedTemporaryFile('w+t') as f:
        subprocess.run(f'{pdb2fasta} -mol prot -het 1 {pdbfile} > {f.name}', shell=True)
        f.seek(0)
        header, seq = next(read_fasta(f.name))
        return header, seq


# read index
with open(index, 'rb') as f:
    index = pickle.load(f)

all_entries = list(index.keys())

########################################################################
os.makedirs(f'{outdir}/all', exist_ok=True)

# STEP 1a: randomly make dimers
n_dimers = 10
dimers_paths = []
while len(dimers_paths) < n_dimers:
    contacts_to_check = []
    while len(contacts_to_check) < 1000:
        picked_entry = random.choice(all_entries)
        pdbs_of_entry = []
        for chain_name, pdbs in index[picked_entry].items():
            for pdb in pdbs:
                pdbs_of_entry.append(pdb)
        if len(pdbs_of_entry) < 2:
            continue
        assembly_groups = []
        assembly_of_pdb = lambda pdb: pdb.split('-')[1] # '119l-a1-m1-cA.pdb' -> 'a1'
        for assembly_name, assembly_group in itertools.groupby(pdbs_of_entry, key=assembly_of_pdb):
            chains = list(assembly_group)
            if len(chains) > 1:
                assembly_groups.append(chains)

        if assembly_groups == []:
            continue

        picked_group = random.choice(assembly_groups)
        dimers_of_group = list(itertools.combinations(picked_group, 2))
        contacts_to_check += dimers_of_group
   
    contacts_to_check = contacts_to_check[:1000]
    tmp = []
    for dimer in contacts_to_check:
        if dimer[1] < dimer[0]:
            dimer = (dimer[1], dimer[0])
        tmp.append(dimer)
    contacts_to_check = tmp
    contacts_to_check = [dimer_pdb_path(p) for p in contacts_to_check]
    dimers_contacting = check_contact_many(pairs_list=contacts_to_check,
                                           thresh_max_dist=8,
                                           thresh_min_pairs=10)
    i = 0
    for dimer, is_contacting in zip(contacts_to_check, dimers_contacting):
        if is_contacting:
            L1 = len(seq_of_pdb(dimer[0])[1])
            L2 = len(seq_of_pdb(dimer[1])[1])
            #if L1 == L2:
            #    continue
            #print(dimer, L1, L2)
            if L1 > L2:
                dimer = (dimer[1], dimer[0])
            if not dimer in dimers_paths:
                dimers_paths.append(dimer)
                i += 1
    print(f'found additional {i} dimers...running total {len(dimers_paths)}...')

dimers_paths = sorted(dimers_paths)[:10]

# note: each dimer is a size-two tuple of paths to pdb chains

def dimer_name(path1, path2):
    c1 = path1.split('/')[-1].rstrip(".pdb")
    c2 = path2.split('/')[-1].rstrip(".pdb")
    return f'{c1}_{c2}'

dimers = []
# STEP 1b: write all/dimers.txt
with open(f'{outdir}/all/dimers.txt', 'w') as f:
    for c1, c2 in dimers_paths:
        name = dimer_name(c1, c2)
        f.write(name)
        f.write('\n')
        dimers.append(name)


########################################################################################




def write_dimer_fastas(outdir, dimers_paths):

    # write longer seq of each dimer to all/seqs_longer.fasta and shorter of each dimer to all/seqs_shorter.fasta
    with open(f'{outdir}/seqs_shorter.fasta', 'w') as short, open(f'{outdir}/seqs_longer.fasta', 'w') as longer:
        for dimer_path in dimers_paths:
            dimer_named = dimer_name(*dimer_path)
            c1_path, c2_path = dimer_path
            c1_name = dimer_named.split('_')[0]
            c2_name = dimer_named.split('_')[1]
            _, seq1 = seq_of_pdb(c1_path)
            _, seq2 = seq_of_pdb(c2_path)

            longer_seq = None
            shorter_seq = None
            longer_name = None
            shorter_name = None

            if len(seq1) < len(seq2):
                longer_seq = seq2
                longer_name = c2_name
                shorter_seq = seq1
                shorter_name = c1_name
            else:
                longer_seq = seq1
                longer_name = c1_name
                shorter_seq = seq2
                shorter_name = c2_name
            
            short.write(f'>{dimer_named};\t{shorter_name}\n')
            short.write(f'{shorter_seq}\n')
            longer.write(f'>{dimer_named};\t{longer_name}\n')
            longer.write(f'{longer_seq}\n')

write_dimer_fastas(outdir=f'{outdir}/all', dimers_paths=dimers_paths)

########################################################################################


def dimers_of_cluster(dimers, cluster_center, assignments):
    cluster = []
    for dimer, assignment in zip(dimers, assignments):
        if assignment == cluster_center:
            cluster.append(dimer)
    return cluster


# STEP 3a: generate fake clusters and cluster assignments
#N_clusters = random.randint(15, 20)
N_clusters = 500
cluster_centers = random.choices(dimers, k=N_clusters)
cluster_assigns = []
for dimer in dimers:
    assigned = random.choice(cluster_centers)
    if dimer in cluster_centers:
        assigned = dimer
    cluster_assigns.append(assigned)


membership = {}

# STEP 3b: generate fake groups and group assignments
for cluster_center in cluster_centers:
    cluster = dimers_of_cluster(dimers, cluster_center, cluster_assigns)
    N_groups = random.randint(1,5)
    if N_groups > len(cluster):
        N_groups = len(cluster)
    group_centers = random.choices(cluster, k=N_groups)
    #group_assigns = []
    for dimer in cluster:
        assigned = random.choice(group_centers)
        if dimer in group_centers:
            assigned = dimer
        #group_assigns.append(assigned)
        membership[dimer] = {'clust': cluster_center, 'group': assigned}


# STEP 3c: write membership and relate to files
os.makedirs(f'{outdir}/cluster', exist_ok=True)
groups = set()
clusters = cluster_centers
with open(f'{outdir}/cluster/membership.tsv', 'w') as f:
    for dimer in membership:
        cluster = membership[dimer]['clust']
        group = membership[dimer]['group']
        f.write(f'{dimer}\t{cluster}\t{group}\n')
        groups.add(group)

with open(f'{outdir}/cluster/clusters.txt', 'w') as f:
    for cluster in clusters:
        f.write(f'{cluster}\n')

with open(f'{outdir}/cluster/groups.txt', 'w') as f:
    for group in groups:
        f.write(f'{group}\n')


########################################################################################
os.makedirs(f'{outdir}/nonredundant', exist_ok=True)

# STEP 4a: pick fake reps of each group for nonredundant
reps = []
for group_name in groups:
    group_members = []
    for dimer, info in membership.items():
        if info['group'] == group_name:
            group_members.append(dimer)
    rep = random.choice(group_members)
    reps.append(rep)
reps = sorted(list(reps))

# STEP 4b: write reps
with open(f'{outdir}/nonredundant/dimers.txt', 'w') as f:
    for rep in reps:
        f.write(f'{rep}\n')

with open(f'{outdir}/nonredundant/chains.txt', 'w') as f:
    rep_chains = set()
    for rep in reps:
        rep_chains.add(rep.split('_')[0])
        rep_chains.add(rep.split('_')[1])
    rep_chains = list(sorted(rep_chains))
    for chain in rep_chains:
        f.write(f'{chain}\n')



# STEP 4d: write rep fasta
reps_paths = []
for rep in reps:
    i = dimers.index(rep)
    path = dimers_paths[i]
    reps_paths.append(path)
        
write_dimer_fastas(outdir=f'{outdir}/nonredundant', dimers_paths=reps_paths)


## get dimers, each named by shorter seq first
#dimers_named_by_shorter_chain = []
#with open(f'{outdir}/nonredundant/seqs_shorter.fasta', 'r') as f:
#    for line in f:
#        if line.startswith('>'):
#            init_dimer_name = line.split(';')[0].lstrip('>')
#            shorter_chain = line.split(';')[1].strip()
#            if init_dimer_name.startswith(shorter_chain):
#                dimers_named_by_shorter_chain.append(init_dimer_name)
#            else:
#                spl = init_dimer_name.split('_')
#                new_dimer_name = f'{spl[1]}_{spl[0]}'
#                dimers_named_by_shorter_chain.append(new_dimer_name)
##dimers_named_by_shorter_chain.sort()
#
#
#
#with open(f'{outdir}/nonredundant/dimers_named_by_shorter_chain.txt', 'w') as f:
#    for dimer in dimers_named_by_shorter_chain:
#        f.write(f'{dimer}\n')


########################################################################################
os.makedirs(f'{outdir}/pdb/div', exist_ok=True)

source_pdbs = set()
for dimer in dimers_paths:
    source_pdbs.add(dimer[0])
    source_pdbs.add(dimer[1])

getdiv = lambda pdbpath: pdbpath.split('/')[-1][1:3]

# pdb tars 
pdb_names = []
with tempfile.TemporaryDirectory() as td:
    for div_name, div_pdbs in itertools.groupby(source_pdbs, key=getdiv):
        divtmpdir = f'{td}/{div_name}' 
        os.makedirs(divtmpdir, exist_ok=True)
        for pdb in div_pdbs:
            shutil.copy(pdb, divtmpdir)
            pdb_name = pdb.split('/')[-1].rstrip('.pdb')
            pdb_names.append(pdb_name)
        subprocess.run(f'tar -czf {outdir}/pdb/div/{div_name}.tar.gz -C {td} {div_name}', shell=True)

# pdb index 
with open(f'{outdir}/all/chains.txt', 'w') as findex:
    for pdb_name in pdb_names:        
        findex.write(f'{pdb_name}\n')


seqs = []
for header, seq in read_fasta(f'{outdir}/all/seqs_shorter.fasta'):
    header = header.split()[1]
    seqs.append( (header, seq) )
for header, seq in read_fasta(f'{outdir}/all/seqs_longer.fasta'):
    header = header.split()[1]
    seqs.append( (header, seq) )

# sort by header
seqs.sort(key=lambda seq: seq[0])

with open(f'{outdir}/pdb/seqs.fasta', 'w') as f:
    for header, seq in seqs:
        f.write(f'>{header}\n')
        f.write(f'{seq}\n')

########################################################################################

os.makedirs(f'{outdir}/extra', exist_ok=True)

for num in [30, 40, 50, 60, 70, 80]:
    N = int(len(reps)*num/100)
    samples = random.choices(reps, k=N)
    with open(f'{outdir}/extra/cluster{num}_dimers.txt', 'w') as f:
        for s in samples:
            f.write(f'{s}\n')
    with open(f'{outdir}/extra/cluster{num}_dimers_longer.fasta', 'w') as fo:
        for header, seq in read_fasta(f'{outdir}/all/seqs_longer.fasta'):
            d = header.split(';')[0].lstrip('>')
            if d in samples:
                fo.write(f'{header}\n')
                fo.write(f'{seq}\n')
    with open(f'{outdir}/extra/cluster{num}_dimers_shorter.fasta', 'w') as fo:
        for header, seq in read_fasta(f'{outdir}/all/seqs_shorter.fasta'):
            d = header.split(';')[0].lstrip('>')
            if d in samples:
                fo.write(f'{header}\n')
                fo.write(f'{seq}\n')

########################################################################################

now = str(datetime.utcnow()) 

with open(f'{outdir}/all/chains.txt', 'r') as f:
    all_num_chains = len(f.readlines())
with open(f'{outdir}/all/dimers.txt', 'r') as f:
    all_num_dimers = len(f.readlines())
with open(f'{outdir}/nonredundant/chains.txt', 'r') as f:
    nr_num_chains = len(f.readlines())
with open(f'{outdir}/nonredundant/dimers.txt', 'r') as f:
    nr_num_dimers = len(f.readlines())




with open(f'{outdir}/summary.txt', 'w') as f:
    f.write(f'LAST_UPDATE_UTC\t{now}\n')
    f.write(f'FULL_SET_NUM_CHAIN\t{all_num_chains}\n')
    f.write(f'FULL_SET_NUM_DIMER\t{all_num_dimers}\n')
    f.write(f'NONRED_NUM_CHAIN\t{nr_num_chains}\n')
    f.write(f'NONRED_NUM_DIMER\t{nr_num_dimers}\n')




