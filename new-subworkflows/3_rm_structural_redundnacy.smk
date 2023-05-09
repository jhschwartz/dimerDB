import os
import sys 
import pickle
from helpers import list_seq_clusters

files = config['paths']['pipeline_files']
out = config['paths']['out']

clustersfile = files+'/seq_cluster/??'
seq_clusters = list_seq_clusters(clustersfile)


intracluster_distance_matrix_pkl = files+'/cluster_data/{seq_cluster}/USalign_distances.pkl'
intracluster_representative = files+'/cluster_data/{seq_cluster}/representative.txt'
intracluster_distances_stored_done = files+'/cluster_data/{seq_cluster}/distances_stored.done'
final_homodimers = os.path.join(out, 'homodimers.txt')





rule all:
    input:
        final_homodimers
    run:
        raise NotImplementedError #TODO


rule compute_distances:
    output:
        pklfile = intracluster_distance_matrix_pkl
    run:
        raise NotImplementedError #TODO


rule cluster_poses_and_choose_rep:
    input:
        pklfile = intracluster_distance_matrix_pkl
    output:
        txtfile = intracluster_representative
    run:
        raise NotImplementedError #TODO


rule store_distances:
    input:
        pklfile = intracluster_distance_matrix_pkl
    output:
        done = intracluster_distances_stored_done
    run:
        raise NotImplementedError #TODO


rule collate_final_db:
    input:
        rep_txtfiles = expand(intracluster_representative, seq_cluster=seq_clusters)
    output:
        final_homodimers
    run:
        raise NotImplementedError #TODO

        






