#!/bin/bash

set -e;

if [[ $1 == "-n" ]]; then
    dryrun="-n";
fi

. ./set_vars.sh

env/bin/python -m snakemake -s workflow/0_local_pdb.smk --profile $smk_profile $dryrun;

# subflow1: run pre-contact calculation, so only until "needed calcs" 
env/bin/python -m snakemake -s workflow/1_derive_dimers.smk --until pairs_needing_contact_calc --profile $smk_profile $dryrun;

# subflow1: run large contact calculations as non-snakemake, more-efficient process
?

# subflow1: rerun lookup and through to end now that large contact calculations are in the database.
env/bin/python -m snakemake -s workflow/1_derive_dimers.smk --force lookup_contacts --profile $smk_profile $dryrun;


env/bin/python -m snakemake -s workflow/2_seq_cluster.smk --profile $smk_profile $dryrun;

env/bin/python -m snakemake -s workflow/3_rm_structural_redundancy.smk --profile $smk_profile $dryrun;

env/bin/python -m snakemake -s workflow/4_prep_final_data.smk --profile $smk_profile $dryrun;




