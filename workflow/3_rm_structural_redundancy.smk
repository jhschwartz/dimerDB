import os
import sys
import itertools
import numpy as np
from sortedcontainers import SortedSet
from datetime import datetime

configfile: 'config.yaml'

sys.path.append('scripts')
import name_pdb
import tmscore_database as TMDB
import unredundant as unred
from wrap_tmscore import calculate_many_dimers_TM_score

subworkflow_done = config['workflow']['3_rm_structural_redundancy_done']
intermediates = config['paths']['intermediates_dir']
lib_path = config['paths']['lib'] # TODO put in config, match to sub0, must contain folder "rcsb"

tm_distances_database = config['paths']['tmscore_db']


max_threads = config['workflow_params']['max_threads']
cluster_thresh = config['workflow_params']['dimer_clustering_threshold']


cl_base = os.path.join(intermediates, 'cluster', 'seq_cluster')
outfile = {
    'cluster_index': os.path.join(intermediates, 'cluster', 'cluster_index.txt'), 
    'calc_dists_stored_done': os.path.join(intermediates, 'cluster', 'calcs_stored.done'),
    'seq_cluster': {
        'template_dists_lookup': os.path.join(cl_base, '{cluster_name}', 'dist_lookup.tsv'),
        'template_dists_calc': os.path.join(cl_base, '{cluster_name}', 'dist_calc.tsv'),
        'template_mem_dimers': os.path.join(cl_base, '{cluster_name}', 'members_dimers.txt'),
        'template_reps_list': os.path.join(cl_base, '{cluster_name}', 'representatives.txt'),
        'template_reps_table': os.path.join(cl_base, '{cluster_name}', 'representation.tsv'),
    },
    'resolu': os.path.join(intermediates, 'resolu.idx'),
    'methods': os.path.join(intermediates, 'pdb_entry_type.txt')
}



## define samples as cluster names
with open(outfile['cluster_index'], 'r') as f:
    cluster_names = [line.rstrip() for line in f]



## helper func to load dimers of seq cluster
def get_dimers(cluster_name):
    dimer_file = outfile['seq_cluster']['template_mem_dimers']\
                        .format(cluster_name=cluster_name)
    with open(dimer_file, 'r') as f:
        dimers = [line.rstrip() for line in f]
        return dimers
    


## helper func to scale threads by num_dimers^2/10000 for resource allocation
def threads_func(cluster_name):
    n = len(get_dimers(cluster_name))
    cores = min(int(n*n/10), max_threads)
    return max(cores, 1)




rule all:
    input:
        expand(outfile['seq_cluster']['template_reps_list'], cluster_name=cluster_names),
        expand(outfile['seq_cluster']['template_reps_table'], cluster_name=cluster_names),
        outfile['calc_dists_stored_done']
    output:
        done = subworkflow_done
    run:
        with open(output.done, 'w') as f:
            f.write('Subworkflow 3 done at ')
            f.write(str(datetime.utcnow()))



# shadow inputs {tm_distances_database} and members_dimers.txt
rule lookup_distances:
    '''
    This rule retrieves the already-computed and stored dimer-dimer 
    TM score distances from the tm_distances_database. This prevents
    unnecessary TM-score computation which is a major bottleneck in
    this pipeline. In an update run, it is expected that this rule would
    find all distances except those of the most recent PDB release; in
    an initial run, this would yield nothing. 
    This rule runs in a separate iteration for each sequence cluster. 
    '''
    output:
        dists_lookup = outfile['seq_cluster']['template_dists_lookup'] 
    run:
        # retrieve dimers
        dimers_of_seq_clust = sorted(get_dimers(wildcards.cluster_name))

        # construct all pairings
        dimer_pairs = itertools.combinations(dimers_of_seq_clust, 2)

        # lookup scores
        scores = []
        for found_dimer_pair in TMDB.lookup_scores(dimer_pairs, tm_distances_database):
            print(found_dimer_pair)
            dimer1, dimer2, score1to2, score2to1 = found_dimer_pair 
            scores.append( (dimer1, dimer2, score1to2, score2to1) )

        # write found scores
        scores = sorted(scores, key=lambda score: (score[0], score[1]) )
        with open(output.dists_lookup, 'w') as f:
            for dimer1, dimer2, score1, score2 in scores:
                f.write(f'{dimer1}\t{dimer2}\t{score1}\t{score2}\n')







# shadow input members_dimers.txt
rule compute_distances:
    '''
    This rule calculates dimer-dimer structural TM score distances 
    for all dimer-dimer pairings not found in the tm_distances_database.
    This rule runs in a separate iteration for each sequence cluster. 
    '''
    input:
        dists_lookup = outfile['seq_cluster']['template_dists_lookup']
    output:
        dists_calc = outfile['seq_cluster']['template_dists_calc']
    threads: 
        lambda wildcards: threads_func(wildcards.cluster_name) 
    run:
        ## STEP 1: figure out which dimers need calc
        dimers_of_clust = sorted(get_dimers(wildcards.cluster_name))
       
        dimer_pairs = itertools.combinations(dimers_of_clust, 2)
        
        with open(input.dists_lookup, 'r') as f:
            dists_found = SortedSet()
            for line in f:
                dimer1, dimer2 = line.split()[:2]
                dists_found.add( (dimer1, dimer2) )
        
        dimer_pairs_need_calc = []
        for dp in dimer_pairs:
            if not dp in dists_found:
                dimer_pairs_need_calc.append(dp) 


        ## STEP 2: perform calc
        scores = calculate_many_dimers_TM_score(  dimer_pairs=dimer_pairs_need_calc, 
                                                  usalign_exe='bin/USalign/USalign',
                                                  lib_path=lib_path,
                                                  cores=threads                      )

        ## STEP 3: write results
        with open(output.dists_calc, 'w') as f:
            for (dimer1, dimer2), (score1, score2) in zip(dimer_pairs_need_calc, scores):
                f.write(f'{dimer1}\t{dimer2}\t{score1}\t{score2}\n')



rule download_resolu_file:
    '''
    downloads resolu.idx from RCSB for getting xray structures'
    resolution in angstroms. Used by rule cluster_poses_and_choose_rep.
    '''
    output:
        resolu = outfile['resolu']
    run:
        shell(''' wget -O {output.resolu} https://files.wwpdb.org/pub/pdb/derived_data/index/resolu.idx ''')
        if config.get('test'):
            shell(''' echo "3CC3  ;    3.0" >> {output.resolu} ''')



rule download_methods_file:
    '''
    downloads pdb_entry_type from RCSB for getting structures' experiment
    type. Used by rule cluster_poses_and_choose_rep to preference xray structures.
    '''
    output:
        methods = outfile['methods']
    run:
        shell(''' wget -O {output.methods} https://files.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt ''')
        if config.get('test'):
            shell(''' echo "3cc3\tprot\tdiffraction" >> {output.methods} ''')



rule cluster_poses_and_choose_rep:
    '''
    This rule reads the looked-up and calculated dimer-dimer TM score
    structural distances from the two rules above, and from this follows
    the strucutral redundancy removal routine defined in 
    scripts/unredundant.py to perform higherarchial agglomerative 
    clutering and representative-choosing from each resutlant structural
    cluster of dimers. This yields a set of nonredundant dimer structures
    for each sequence cluster. (this rule is run once for each sequence
    cluster - these seq clusters are each a snakemake "sample")
    '''
    input:
        dists_lookup = outfile['seq_cluster']['template_dists_lookup'],
        dists_calc = outfile['seq_cluster']['template_dists_calc'],
        resolu = outfile['resolu'],
        methods = outfile['methods']
    output:
        reps_txt = outfile['seq_cluster']['template_reps_list'],
        reps_tsv = outfile['seq_cluster']['template_reps_table']
    threads: 
        lambda wildcards: threads_func(wildcards.cluster_name) 
    run:
        # retrieve all dimer names
        dimers_of_seq_clust = sorted(get_dimers(wildcards.cluster_name))
        
        # load all distances into dictionary
        pair2key = lambda pair_tuple: name_pdb.name_dimer(*pair_tuple)
        scores = {}
        with open(input.dists_lookup, 'r') as f:
            for line in f:
                dimer1, dimer2, score1, score2 = line.split()
                pair_key = pair2key( (dimer1, dimer2) )
                scores[pair_key] = max(score1, score2)
        with open(input.dists_calc, 'r') as f:
            for line in f:
                dimer1, dimer2, score1, score2 = line.split()
                pair_key = pair2key( (dimer1, dimer2) )
                scores[pair_key] = max(score1, score2)

        combs = len(list(itertools.combinations(dimers_of_seq_clust, 2)))
        if combs != len(scores.keys()):
            raise ValueError('was unable to load either a lookup or calc score '    \
                          + f'for all dimer pairings in {wildcards.cluster_name}. ' \
                          + f'There are {combs} dimer combinations but only'        \
                          + f' {len(scores.keys())} scores found'                    )
        
        # init dist matrix
        L = len(dimers_of_seq_clust)
        dist_mat = np.zeros([L, L], dtype=np.float16)

        # fill dist matrix
        for i in range(L):
            for j in range(i+1, L):
                dimer1 = dimers_of_seq_clust[i]
                dimer2 = dimers_of_seq_clust[j]
                pair_key = pair2key( (dimer1, dimer2) )
                dist = 1 - float(scores[pair_key]) # invert TMscore to dist
                dist_mat[i,j] = dist
                dist_mat[j,i] = dist


        # run clustering upon dist mat
        dimer_group_labels = unred.run_cluster_from_matrix(npmat=dist_mat, 
                                                       thresh=cluster_thresh)
        clusters = set(dimer_group_labels)

        # choose reps of each structural cluster (group)
        reps = {}
        for c_num in clusters:
            group_dimers = unred.get_cluster_dimers(  cluster_num=c_num,
                                                      cluster_labels=dimer_group_labels,
                                                      all_dimers=dimers_of_seq_clust     )
            rep = unred.choose_rep(  dimers=group_dimers,
                                     all_dimers=dimers_of_seq_clust, 
                                     dist_mat=dist_mat, 
                                     lib_path=lib_path, 
                                     resolu_path=input.resolu,
                                     methods_file=input.methods       )
            reps[c_num] = rep

        # write reps_txt
        with open(output.reps_txt, 'w') as f:
            for rep in reps.values():
                f.write(f'{rep}\n')

        # write reps_tsv
        with open(output.reps_tsv, 'w') as f:
            for dimer_name, dimer_label in zip(dimers_of_seq_clust, dimer_group_labels):
                #dimer_name = name_pdb.name_dimer(*dimer_t)
                rep = reps[dimer_label]
                f.write(f'{dimer_name}\t{rep}\n')





# shadow input {tm_distances_database}
rule store_distances:
    '''
    To expedite future updates of the database, this rule stores
    the distances that required calculation (not found via lookup 
    in the database) in the tm_distances_database. Runs once for
    all sequence clusters, without separate smk samples. 
    '''
    input:
        calc_files = expand(  outfile['seq_cluster']['template_dists_calc'],
                              cluster_name=cluster_names                     )
    output:
        done = outfile['calc_dists_stored_done']
    run:
        dimers_scores = []
        for calc_file in input.calc_files:
            with open(calc_file, 'r') as f:
                for line in f:
                    dimer1, dimer2, score1, score2 = line.split()
                    dimers_scores.append( (dimer1, dimer2, score1, score2) )

        TMDB.update_db(new_pairs_scores=dimers_scores, db_path=tm_distances_database)

        shell(''' touch {output.done} ''')







