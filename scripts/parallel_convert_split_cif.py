'''
parallel_convert_split_cif.py - script to call cif2pdb (from USalign of the Pyle Lab)
    upon mmcif protein assemblies to both 1) convert each to .pdb format and 2) split
    each assembly into individual chains as defined by chainID and model number. 
    Functions below perform this recursively on all .cif files in a given parent
    directory and additionally rename the resulting .pdb files according to a standard
    scheme for this pipeline.

Written by Jacob Schwartz (jaschwa@umich.edu) in Janurary/February 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

These functions are unitested by test/test_parallel_convert_split_cif.py and are
    passing as of 2/1/2023.
This work requires python >= 3.8
'''
import os
import shutil
import glob
import re
import argparse
import subprocess
from multiprocessing import Pool
import tempfile
from name_pdb import name_pdb_file


def run_cif2pdb(cif: str, outdir: str, prefix: str, cif2pdb_exe: str) -> None:
    with open(os.devnull, 'w') as NULL: # NULL to silence stdout of cif2pdb
        subprocess.run(f'{cif2pdb_exe} -split 1 -mol 1 {cif} {outdir}/{prefix}', shell=True, check=True, stdout=NULL, stderr=subprocess.STDOUT)


def rename_resulting_pdbs(dir_: str, assembly_ID: str, prefix: str) -> None:
    pdbfiles = glob.glob(f'{dir_}/{prefix}*.pdb')
    new_paths = []
    if len(pdbfiles) < 1:
        raise RuntimeError(f'was not able to find the resulting files of cif2pdb in {dir_} with prefix {prefix}.')
    pattern_model1 = fr'{prefix}.+?\.pdb'
    pattern_model_other = fr'{prefix}.+?-\d+?\.pdb'
    for fn in pdbfiles:
        base_fn = os.path.basename(fn)
        if re.search(pattern_model_other, fn):
            model = re.findall(fr'{prefix}.+?-(\d+?)\.pdb', base_fn)[0]
            chain = re.findall(fr'{prefix}(.+?)-\d+?\.pdb', base_fn)[0]
        elif re.search(pattern_model1, fn):
            model = 1
            chain = re.findall(fr'{prefix}(.+?)\.pdb', base_fn)[0]   
        else:
            raise ValueError(f'the encountered pdb file, {fn}, which is supposed to be a cif2pdb result of prefix {prefix} in dir {dir_}, does not search match expected pattern.')
        
        new_fn = os.path.join(dir_, name_pdb_file(pdb_base=prefix, assembly=assembly_ID, model=model, chain=chain))
        os.rename(fn, new_fn)
        new_paths.append(new_fn)
    return new_paths


def process_helper(cif_fn: str, cif2pdb_exe: str) -> None:
    pdb_base = os.path.basename(cif_fn)[:4]
    assembly_ID = re.findall(r'assembly(\d+)\.cif', cif_fn)[0]
    dir_target = os.path.dirname(os.path.realpath(cif_fn))
    with tempfile.TemporaryDirectory() as td:
        run_cif2pdb(cif=cif_fn, outdir=td, prefix=pdb_base, cif2pdb_exe=cif2pdb_exe)
        named_pdbs = rename_resulting_pdbs(dir_=td, assembly_ID=assembly_ID, prefix=pdb_base)
        for f in named_pdbs:
            shutil.move(f, dir_target)



def parallel_convert_split_rename_cifs(parent_dir: str, cif2pdb_exe: str, threads: int) -> None:
    cif_iterator = glob.iglob(f'{parent_dir}/**/*-assembly*.cif', recursive=True)
    args = ((cif, cif2pdb_exe) for cif in cif_iterator)
    with Pool(processes=threads) as p:
        p.starmap(process_helper, args)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--parent-dir', required=True)
    parser.add_argument('-e', '--cif2pdb-exe-path', required=True)
    parser.add_argument('-t', '--threads', type=int, required=True)
    args = parser.parse_args()
    parallel_convert_split_rename_cifs(args.parent_dir, args.cif2pdb_exe_path, args.threads)

