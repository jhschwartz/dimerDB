import os
import glob
import re
import argparse
import subprocess
from multiprocessing import Pool


def run_cif2pdb(target: str, outprefix: str, cif2pdb_exe: str) -> None:
    outdir = os.path.dirname(target)
    with open(os.devnull, 'w') as NULL:
        subprocess.run(f'{cif2pdb_exe} -split 1 -mol 1 {target} {outdir}/{outprefix}', shell=True, check=True, stdout=NULL, stderr=subprocess.STDOUT)


def rename_resulting_pdbs(dir_: str, prefix: str) -> None:
    pdbfiles = glob.glob(f'{dir_}/{prefix}*.pdb')
    if len(pdbfiles) < 1:
        raise RuntimeError(f'was not able to find the resulting files of cif2pdb in {dir_} with prefix {prefix}.')
    pattern_model1 = fr'{prefix}.+?\.pdb'
    pattern_model_other = fr'{prefix}.+?-\d+?\.pdb'
    for fn in pdbfiles:
        base_fn = os.path.basename(fn)
        if re.search(pattern_model_other, fn):
            chain = re.findall(fr'{prefix}(.+?)-\d+?\.pdb', base_fn)[0]
            model = re.findall(fr'{prefix}.+?-(\d+?)\.pdb', base_fn)[0]
            new_fn = f'{dir_}/{prefix}-{chain}-{model}.pdb'
            os.rename(fn, new_fn)
        elif re.search(pattern_model1, fn):
            chain = re.findall(fr'{prefix}(.+?)\.pdb', base_fn)[0]   
            new_fn = f'{dir_}/{prefix}-{chain}-1.pdb'
            os.rename(fn, new_fn)
        else:
            raise ValueError(f'the encountered pdb file, {fn}, which is supposed to be a cif2pdb result of prefix {prefix} in dir {dir_}, does not search match expected pattern.')
        


def process_helper(cif_fn: str, cif2pdb_exe: str) -> None:
    prefix = os.path.basename(cif_fn)[:4]
    dir_ = os.path.dirname(os.path.realpath(cif_fn))
    run_cif2pdb(cif_fn, prefix, cif2pdb_exe)
    rename_resulting_pdbs(dir_, prefix)


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

