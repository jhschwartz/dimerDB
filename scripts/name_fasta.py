'''
name_fasta.py - given a homodimer uniparc id function gives the path
                to a fasta file(s) that represent(s) the sequence(s). Sequences are sourced from a prior 
                intermediate passed as filtering_dir. This dir should be defined in the config file of the 
                pipeline. Examples (from the initial pipeline config in development) are: 
                    - dimerDB/intermediates/homodimer_filtering

Written by Jacob Schwartz (jaschwa@umich.edu) in January 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

This function is unittested by test/test_name_fasta.py and passing as of 1/12/2023.
This work requires python >= 3.8
'''
import os



def get_homodimer_fasta(uniparc_id, filtering_dir):
    div = uniparc_id[-2:]
    path = os.path.join(filtering_dir, div, uniparc_id, 'seq.fasta')
    if not os.path.exists(path):
        raise FileNotFoundError(f'unable to find fasta for homodimer with uniparc id {uniparc_id} at path {path}')
    return os.path.realpath(path)
