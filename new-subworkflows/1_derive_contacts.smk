import os
import sys 
import pickle
from helpers import resolve_entries_list, divs_from_entries

configfile: 'config.yaml'

files = config['paths']['pipeline_files']

entries = resolve_entries_list(f'{files}/div')
divs = divs_from_entries(entries)


entry_potential_pairs = files+'/div/{div}/{entry}/potential_pairs.txt'
entry_contacting_pairs = files+'/div/{div}/{entry}/contacting_pairs.txt'



rule all:
    input:
        expand(entry_contacting_pairs, zip, div, entry)
    run:
        raise NotImplementedError # TODO




rule pair_possible_contacting_chains:
    output:
        pairsfile = entry_potential_pairs
    run:
        raise NotImplementedError # TODO



rule check_contacts:
    input:
        pairsfile = entry_potential_pairs
    output:
        contactsfile = entry_contacting_pairs
    run:
        raise NotImplementedError # TODO












