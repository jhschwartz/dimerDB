#!/usr/bin/perl
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



