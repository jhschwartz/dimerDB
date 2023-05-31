#!/bin/bash

# parallel_ungzip_all.pl - script to ungzip the files in a folder (searching recursively) using parallel
#                          processing. Given the directory to search in, and the number of threads to use, 
#                          queues the jobs and uses GNU parallel to perform the ungzips. This is necessary
#                          because in this pipeline there might be too many files to split to fit in a
#                          single shell command, thus we cannot use GNU parallel directly without a queue.
# 
# Written by Jacob Schwartz (jaschwa@umich.edu) in December 2022.
# Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
# https://freddolino-lab.med.umich.edu
# 
# This function is unittested by test/test_parallel_ungzip.py and passing as of 1/12/2023
# This work requires python >= 3.8

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

