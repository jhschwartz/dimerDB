import os
from name_pdb import read_chain_names
import itertools
import pickle


def generate_rcsb_index(rcsb_path, index_file):
    names = []
    emptyfiles = []
    for root, dirs, files in os.walk(rcsb_path):
        for file in files:
            if file.endswith('.pdb'):
                names.append(file)
                path = os.path.join(root, file)
                if os.stat(path).st_size == 0: # empty file
                    emptyfiles.append(path)

    with open(index_file, 'w') as f:
        for name in sorted(names):
            f.write(name+'\n')

    return emptyfiles


def rcsb_index_to_pkl(index_file_txt, index_file_pkl):
    get_pdb = lambda filename: read_chain_names(filename)[0]
    get_chain = lambda filename: read_chain_names(filename)[3]

    whole_pdb = {}
    with open(index_file_txt, 'r') as f: 
        for pdb_name, pdb_group in itertools.groupby(f, get_pdb):
            pdb_entry = {}
            for chain_name, chain_group in itertools.groupby(pdb_group, get_chain):
                chain_entry = [chainfile.strip() for chainfile in list(chain_group)]
                pdb_entry[chain_name] = chain_entry
            whole_pdb[pdb_name] = pdb_entry

    with open(index_file_pkl, 'wb') as f:
        pickle.dump(whole_pdb, f)


