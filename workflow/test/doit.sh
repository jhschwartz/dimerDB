#!/bin/sh



for i in $(ls data/2_seq_cluster/expected_intermediates/cluster/seq_clusters/);
do
    grep $i tmpout_cluster.tsv | awk '{print $2}' > data/2_seq_cluster/expected_intermediates/cluster/seq_clusters/$i/member_chains.txt
done
