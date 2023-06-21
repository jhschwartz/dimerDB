#!/bin/bash

set -e;

split_dir="splitdir$RANDOM";
split_prefix="$split_dir/splitin_";

mkdir $split_dir;
echo $split_dir;


split -l 1000 /expanse/lustre/projects/mia174/jaschwa/dimerDB/lib/pdb_index.txt $split_prefix;

for inputfile in $split_dir/splitin_*;
do
    sbatch -A mia174 -p shared -n 1 -N 1 --mem=2G -t "1:00:00" truncategroup.sh $inputfile;
done



