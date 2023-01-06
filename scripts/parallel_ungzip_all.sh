#!/bin/bash

dir_=$1;
threads=$2;

# see comments below
parallel_batch_size=$3;

# read all *.gz files from recursive find into array
#   technique courtesy of stackoverflow user benjamin-w
#   https://stackoverflow.com/a/54561526/4176019
readarray -d '' gz_files < <(find $dir_ -name "*.gz" -print0);

num_files=${#gz_files[@]};

start_i=0;
while [ "$start_i" -lt "$num_files" ]
do
    # batch arguments to parallel to avoid maxing out bash command length  
    batch=("${gz_files[@]:start_i:parallel_batch_size}");
    
    # pass the batch as a string to parallel ungzip
    batch_as_string="${batch[@]}";
    parallel -j$threads gzip -d {} ::: $batch_as_string;
    
    # increment start for next batch
    batch_size=${#batch[@]};
    start_i=$(($start_i + $batch_size));

done

