import glob 
import shutil

pdbs = [p.split('/')[-1] for p in glob.glob('*/*')]

source_parent = '/nfs/turbo/umms-petefred/jaschwa/dimerDB/scripts/test/data/fake_lib/rcsb_pdb/'


for pdb in pdbs:
    div = pdb[1:3]
    cifs = glob.glob(f'{source_parent}/{div}/{pdb}*.cif')
    chains = glob.glob(f'{source_parent}/{div}/{pdb}*.pdb')

    for cif in cifs:
        shutil.copy(cif, f'{div}/{pdb}/assemblies/') 
    for chain in chains:
        shutil.copy(chain, f'{div}/{pdb}/chains/') 
