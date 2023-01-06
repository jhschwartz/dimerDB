import re
import argparse
import shutil
import glob
import os


def _bundle_name_from_mappingline(line):
    pattern = '\-(bundle[\d]*)\.pdb\:'
    bundle = re.findall(pattern, line)[0]
    return bundle


def _bundle_name_from_filename(filename):
    pattern = '\-(bundle[\d]*).{1}\.'
    bundle = re.findall(pattern, filename)[0]
    return bundle


def _chain_name_from_filename(filename):
    pattern = '\-bundle[\d]*(.{1})\.pdb'
    old_chain = re.findall(pattern, filename)[0]
    return old_chain  


def _read_mappings(mapping_file):
    mappings = {}
    with open(mapping_file, 'r') as mf:
        line = mf.readline()
        while line and not 'bundle' in line:
            line = mf.readline() 
        current_bundle = _bundle_name_from_mappingline(line)
        mappings[current_bundle] = {}
        while line := mf.readline():
            if line.strip() == '':
                continue
            elif 'bundle' in line:
                current_bundle = _bundle_name_from_mappingline(line)
                mappings[current_bundle] = {}
            else:
                remapped_chain, original_chain = line.split()
                mappings[current_bundle][remapped_chain] = original_chain
        return mappings




def _derive_new_name(chainfile, mappings):
    pdb_base_id = chainfile[:4]
    bundle_name = _bundle_name_from_filename(chainfile)
    old_chain = _chain_name_from_filename(chainfile)
    new_chain = mappings[bundle_name][old_chain]
    new_name = f'{pdb_base_id}{new_chain}.pdb'
    return new_name



def _file_is_single_chain(filename, mappings):
    if filename == '7w5z-pdb-bundle3Z.pdb':
        print('hello')
    bundles = mappings.keys()
    bundle_of_file = None
    for b in sorted(bundles):
        if b in filename:
            bundle_of_file = b
    if bundle_of_file == None:
        raise ValueError(f'file {filename} does not contain a valid bundle name for given mappings')
    if filename.split(bundle_of_file)[-1] == '.pdb':
        return False # not just one chain because there is no chain ID after the bundle name
    return True # chain ID exists after bundle name



def _rename_chainfiles_in_dir(dir_, pdb_id):
    mapping_file = f'{dir_}/{pdb_id}-chain-id-mapping.txt'
    mappings = _read_mappings(mapping_file)

    
    # file names will be like "7w5z-pdb-bundle1A.pdb" for bundle1 chain A, etc.
    pdb_files = glob.glob(f'{dir_}/*.pdb')
    for fn in pdb_files:
        if not _file_is_single_chain(fn, mappings):
            continue
        new_name = _derive_new_name(os.path.basename(fn), mappings)
        os.makedirs(f'{dir_}/split_renamed', exist_ok=True)
        shutil.move(fn, f'{dir_}/split_renamed/{new_name}')
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--lib', help='the locations of the oversized pdb lib, after splits done')
    parser.add_argument('--pdb-id', help='the pdb base code we are going to deal with')
    args = parser.parse_args()

    div = args.pdb_id[1:3]
    prot_dir = f'{args.lib}/{div}/{args.pdb_id}'

    # avoid empty pdb archives until RCSB fixes it  
    if not args.pdb_id in ['7nwh', '7nwi', '7nwg']:
        _rename_chainfiles_in_dir(prot_dir, args.pdb_id)


