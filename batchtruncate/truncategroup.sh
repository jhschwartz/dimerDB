#!/bin/sh

index=$1;



function truncatePDB () {
    input=$1;
    output=$2;
    if [[ -f $output ]]; then
        echo "skipping truncation of $input because $output already exists...";
    else
        mkdir -p $(dirname $output);
        grep -P "( CA )|( CB )" $input | cut -c1-54 > $output;
    fi
}

export -f truncatePDB;


for pdb in $(cat $index);
do
    truncatePDB /expanse/lustre/projects/mia174/jaschwa/dimerDB/lib/incl_rcsb_sidechains/rcsb/$pdb /expanse/lustre/projects/mia174/jaschwa/dimerDB/lib/rcsb/$pdb;
done

