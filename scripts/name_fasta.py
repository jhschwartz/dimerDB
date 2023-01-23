'''
name_fasta.py - given a homodimer uniparc id or heterodimer uniparc id pair, these functions give the path
                to a fasta file(s) that represent(s) the sequence(s). Sequences are sourced from a prior 
                intermediate passed as filtering_dir. This dir should be defined in the config file of the 
                pipeline. Examples (from the initial pipeline config in development) are: 
                    - dimerDB/intermediates/homodimer_filtering
                    - dimerDB/intermediates/heterodimer_filtering

Written by Jacob Schwartz (jaschwa@umich.edu) in January 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

This function is unittested by test/test_name_fasta.py and passing as of 1/12/2023.
This work requires python >= 3.8
'''
import os



def get_homodimer_fasta(uniparc_id, filtering_dir):
    div = uniparc_id[-2:]
    path = f'{filtering_dir}/{div}/{uniparc_id}/seq.fasta'
    if not os.path.exists(path):
        raise FileNotFoundError(f'unable to find fasta for homodimer with uniparc id {uniparc_id} at path {path}')
    return os.path.realpath(path)



def get_heterodimer_fasta(dimer_name, filtering_dir):
    uniparc1, uniparc2 = dimer_name.split('-')
    div = uniparc1[-1] + '-' + uniparc2[-1]

    path1 = f'{filtering_dir}/{div}/{dimer_name}/seq1.fasta'
    path2 = f'{filtering_dir}/{div}/{dimer_name}/seq2.fasta'


    if not os.path.exists(path1):
        raise FileNotFoundError(f'unable to find fasta for homodimer {dimer_name} with chain uniparc id {uniparc1} at path {path1}')
    if not os.path.exists(path2):
        raise FileNotFoundError(f'unable to find fasta for homodimer {dimer_name} with chain uniparc id {uniparc2} at path {path2}')
   
    return os.path.realpath(path1), os.path.realpath(path2)

