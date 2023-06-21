#!/bin/sh

cluster_dir="../intermediates/cluster/seq_clusters";
cluster_dir=$(realpath $cluster_dir);

all_pairs_tmp=$(mktemp);

for cluster in $cluster_dir/*;
do
    ../env/bin/python cluster_pairings.py $cluster/members_dimers.txt;
done > $all_pairs_tmp;

split_dir="splits$RANDOM";

split -l 1000 $all_pairs_tmp $splits_dir/pairs_;

for pairsfile in $splits_dir/pairs_*;
do
    sbatch compute_pairs_tms.sh $pairsfile;
done
