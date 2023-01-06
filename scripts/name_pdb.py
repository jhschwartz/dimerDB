import os





def name_pdb_file(pdb_base, chain, lib_path):
    div = pdb_base[1:3]
    path = f'{lib_path}/rcsb_pdb/{div}/{pdb_base}{chain}.pdb'
    if not os.path.exists(path):
        # try oversized
        path = f'{lib_path}/rcsb_oversized/{div}/{pdb_base}/split_renamed/{pdb_base}{chain}.pdb'
        if not os.path.exists(path):
            raise FileNotFoundError(f'was unable to find pdb file {path} for {pdb_base}/{chain}.') 
    return path



def dimer2pdbs(dimer_name, pdb_base_dir):
    pdb0, chain0 = dimer_name.split('-')[0].split('_')
    pdb1, chain1 = dimer_name.split('-')[1].split('_')
    return name_pdb_file(pdb0, chain0, pdb_base_dir), name_pdb_file(pdb1, chain1, pdb_base_dir)
    

