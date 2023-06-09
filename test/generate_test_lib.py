import glob
import os
import shutil
import sys
import subprocess

sys.path.append('../scripts')
from generate_rcsb_index import generate_rcsb_index
from parallel_convert_split_cif import parallel_convert_split_rename_cifs



def generate_test_lib(source_lib, test_entries_file, dest_lib):
    with open(test_entries_file, 'r') as f:
        test_entries = [l.strip() for l in f]

    for te in test_entries:
        div = te[1:3]
        os.makedirs(f'{dest_lib}/rcsb/{div}', exist_ok=True)
        for cif in glob.iglob(f'{source_lib}/rcsb/{div}/{te}*.cif'):
            shutil.copy(cif, f'{dest_lib}/rcsb/{div}')

    # ungzip
    #subprocess.run(f'../scripts/parallel_ungzip_all.sh {dest_lib}/rcsb 8 100', shell=True, check=True)

    # split and name
    parallel_convert_split_rename_cifs(f'{dest_lib}/rcsb', '../bin/USalign/cif2pdb', 8)

    # generate index
    generate_rcsb_index(f'{dest_lib}/rcsb', f'{dest_lib}/pdb_index.txt')

