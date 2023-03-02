'''
read_fasta.py - function to retrieve individual protein entries in a FASTA file as a tuple generator

Written by Jacob Schwartz (jaschwa@umich.edu) in November, 2022.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

These functions are unittested in test/test_read_fasta.py and found to be passing as of 11/17/22.
This work requires python >= 3.8
'''

def read_prot_from_fasta(fasta: str) -> str:
    '''
    From a fasta file, makes a generator of tuples, for which each has the first
    element the fasta header and the second the fasta sequence.
    We use a generator so that we don't have to load huge fasta files into memory.

    :param fasta: str, path to a valid fasta file

    :yields: tuple of size 2, each is a single protein's entry in the fasta file,
                    with the first element of the tuple the header and the 
                    second element the fasta sequence.
    '''
    seq = ''
    with open(fasta, 'r') as f:
        header = f.readline().strip()
        for line in f:
            if seq != '' and line.startswith('>'):
                yield header, seq
                header = line.strip()
                seq = ''
            else:
                seq += line.strip()
        yield header, seq
