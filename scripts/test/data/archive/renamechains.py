import glob 
import os


chainfiles = glob.glob(f'**/chains/*.pdb', recursive=True)

def convert(f):
    pdb, chain, model = f.replace('.pdb', '').split('-')
    return f'{pdb}-a1-m{model}-c{chain}.pdb'

newfiles = [convert(f) for f in chainfiles]

for a,b in zip(chainfiles, newfiles):
    os.rename(a,b)
