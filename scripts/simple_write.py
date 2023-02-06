'''
simple_write.py - a function to write to a file in one line. Used throughout this
    pipeline for information flow between parts of snakemake subworkflows.

Written by Jacob Schwartz (jaschwa@umich.edu) in November, 2022.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu
'''

def simple_write(filename, content=None):
    '''
    Writes a string to a file, handling opening and closing.
    If parameter conent is not supplied, acts as a touch.
    This function overwrites.

    :param filename: str, the path to the file we are creating and writing.
    :param content: str (optional), what we are writing
    '''
    with open(filename, 'w') as f:
        if content is not None:
            f.write(content)
