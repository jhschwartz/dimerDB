import tarfile
import itertools
import os
import glob 


def tar_pdb_chains(indexfile, chains_lib, target_dir):
    '''
    Function to create a .tar.gz for groups of chain .pdb
    files. Specifically, uses the two middle characters of 
    pdb entry IDs to group, so we get files like aa.tar.gz 
    or e6.tar.gz.

    :param indexfile: str, path to a sorted textfile which
                      lists all pdb files contained within
                      chains_lib, one on each line.
                        This file includes divs, e.g. 
                        "aa/1aa1-a1-m1-cA.pdb"
    :param chains_lib: str, path to the directory which
                       contains the divs that contain single
                       chain pdb files in the naming scheme
                       of "1aa1-a1-m1-cA.pdb"
    :param target_dir: str, path to the directory in which
                       the .tar.gz files will be created,
                       i.e. <target_dir>/aa.tar.gz

    :returns: list, strings which are the paths to the created
              .tar.gz files (or, specifically, just the ones
              which exist in target_dir at end of func.)
    '''

    os.makedirs(target_dir, exist_ok=True)
    
    with open(indexfile, 'r') as f:
        chainfiles = []
        for line in f:
            cfile = os.path.basename(line.strip())
            chainfiles.append(cfile)

    middle_two = lambda filename: filename[1:3]
    for key, group in itertools.groupby(chainfiles, key=middle_two):
        tarname = os.path.join(target_dir, f'{key}.tar.gz')
        with tarfile.open(tarname, 'w:gz') as tar:
            for pdb in group:
                if not pdb.endswith('.pdb'):
                    continue
                filename = os.path.basename(pdb)
                path = os.path.join(chains_lib, key, filename)
                if not os.path.exists(path):
                    raise FileNotFoundError('creating tar for pdb entry'  \
                                          +f'group {key} failed upon add' \
                                          +f'of {pdb} because path {path}'\
                                          +f'does not exist!')
                tar.add(path)

    found_tars = glob.glob(f'{target_dir}/*.tar.gz')
    return found_tars
    



