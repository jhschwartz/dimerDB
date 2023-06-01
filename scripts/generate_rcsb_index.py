import os
from name_pdb import read_chain_names, get_div
import itertools
import pickle


def generate_rcsb_index(rcsb_path, index_file):
    names = []
    emptyfiles = []
    for root, dirs, files in os.walk(rcsb_path):
        for file in files:
            if file.endswith('.pdb'):
                names.append(f'{get_div(file)}/{file}')
                path = os.path.join(root, file)
                if os.stat(path).st_size == 0: # empty file
                    emptyfiles.append(path)

    with open(index_file, 'w') as f:
        for name in sorted(names):
            f.write(f'{name}\n')

    return emptyfiles

