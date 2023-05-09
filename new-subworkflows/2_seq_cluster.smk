import os
import sys 

files = config['paths']['pipeline_files']

entries = resolve_entries_list(f'{files}/div')
divs = divs_from_entries(entries)

entry_fasta_done = files+'/div/{div}/{entry}/fastas.done'
all_chains_fasta = files'/no_sample/all_chains.fasta'
clustersfile = files+'/seq_cluster/??'



rule all:
    input:
        clustersfile
    run:
        raise NotImplementedError # TODO


rule write_fastas:
    output:
        done = entry_fasta_done
    run:
        raise NotImplementedError # TODO


rule one_fasta:
    input:
        expand(entry_fasta_done, zip, divs, entries)
    output:
        fasta = all_chains_fasta
    run:
        raise NotImplementedError # TODO


rule seq_cluster:
    input:
        fasta = all_chains_fasta
    output:
        clustersfile
    run:
        raise NotImplementedError # TODO




