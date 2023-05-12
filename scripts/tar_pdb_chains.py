import tarfile
import itertools
import os
import glob 
import name_pdb


def tar_pdb_chains(indexfile, div, chains_lib, target_dir):
    '''
    Function to create a .tar.gz for a group of chain .pdb
    files. Specifically, uses div parameter to look for
    pdb files in chains_lib via indexfile which have div 
    as the pdb entry's middle two char. 

    :param indexfile: str, path to a sorted textfile which
                      lists all pdb files contained within
                      chains_lib, one on each line.
                        This file includes divs, e.g. 
                        "aa/1aa1-a1-m1-cA.pdb"
    :param div: str, the div which defines the group of 
                pdb chains, specifically by the middle two
                chars of the pdb entry id of each chain
    :param chains_lib: str, path to the directory which
                       contains the divs that contain single
                       chain pdb files in the naming scheme
                       of "1aa1-a1-m1-cA.pdb"
    :param target_dir: str, path to the directory in which
                       the .tar.gz files will be created,
                       i.e. <target_dir>/aa.tar.gz

    :returns: lists, strings which are the names of the pdb
              files put into the created tar in target_dir.
    '''

    os.makedirs(target_dir, exist_ok=True)
    
    tar_path = os.path.join(target_dir, f'{div}.tar.gz')

    chains_tarred = []
    with open(indexfile, 'r') as f, tarfile.open(tar_path, 'w:gz') as tar:
        for line in f:
            chain_name = line.strip()
            if chain_name.endswith('.pdb') and name_pdb.get_div(chain_name) == div:
                chain_source =  name_pdb.name_pdb_file(
                                    *name_pdb.read_chain_names(chain_name),
                                    lib_path=chains_lib,
                                    allow_nonexist=False
                                    )
                simple_name = os.path.basename(chain_source)
                tar.add(chain_source, arcname=simple_name)
                chains_tarred.append(simple_name)

    return chains_tarred

            
    



