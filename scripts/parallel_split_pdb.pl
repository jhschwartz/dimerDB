#!/usr/bin/perl

# parallel_split_pdb.pl - script to split the pdb files in a folder (searching recursively) using parallel
#                         processing. Given the directory to search in, the number of threads to use, and
#                         the path the the PDBParser split_chain executable, queues the jobs and uses GNU
#                         parallel to perform the splits. This is necessary because in this pipeline there
#                         might be too many files to split to fit in a single shell command, thus we cannot
#                         use GNU parallel directly without a queue. 
# 
# Written by Jacob Schwartz (jaschwa@umich.edu) in December 2022.
# Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
# https://freddolino-lab.med.umich.edu
# 
# This function is unittested by test/test_parallel_split_chain.py and passing as of 1/12/2023.
# This work requires python >= 3.8

use strict;
use warnings;

my $parent_dir = shift;
my $threads = shift;
my $split_chain_exe = shift;


my @pdbs = `find $parent_dir -name "*.pdb*"`;
my @pdbs_not_split = ();

foreach my $pdb (@pdbs) {
    if ($pdb =~ /\/.{4}\./ or $pdb =~ /\/.{4}\-pdb\-bundle\d{1}\.pdb/) {
        push @pdbs_not_split, $pdb;
    }
}


# do for all pdbs
while (scalar @pdbs_not_split > 0) {
    
    # build list of 1000 pdbs (or less)
    my @_1000_pdbs = ();
    while ( scalar @pdbs_not_split > 0 and scalar @_1000_pdbs < 1000 ) {
        my $pdb = pop @pdbs_not_split;
        push @_1000_pdbs, $pdb;
    }

    # perform parallel upon the 1000 pdbs
    chomp @_1000_pdbs;
    my $pdbstr = join(" ", @_1000_pdbs);
    print `parallel -j$threads $split_chain_exe {} ::: $pdbstr`;

}



