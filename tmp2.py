import numpy as np

file = 'logs/expand_clean_uniparc2others/expand_clean_uniparc2others-4141859.out'

with open(file) as f:
    lines = [f.strip() for f in f.readlines() if f.startswith('enc')]
    scores = [float(l.split()[5].split('=')[1]) for l in lines]
    parcs = [l.split()[3].split('=')[1][:-1] for l in lines]
    chains = [l.split()[4].split('=')[1][:-1] for l in lines]
    

parcs = [f'lib/uniparc/{p[-2:]}/{p}.fasta' for p in parcs]
chains = [f'lib/rcsb/{c[1:3]}/{c.replace("_","-")}.pdb' for c in chains]



print('len', len(scores))
print('mean', np.mean(scores))

#np.save('scores.npy', scores)



#for i, score in enumerate(scores):
#    if score > 0.85:
#        print(lines[i])


import tempfile
import subprocess
from scripts.read_fasta import read_prot_from_fasta

def get_pdb_seq(pdb):
    with tempfile.TemporaryDirectory() as td:
        cmd = f'bin/USalign/pdb2fasta {pdb} > {td}/fasta'
        subprocess.run(f'bin/USalign/pdb2fasta {pdb} > {td}/fasta', shell=True, text=True, capture_output=True)
        _, seq = next(read_prot_from_fasta(f'{td}/fasta'))
        return seq



good = []
short = 0
allX = 0
manyX = 0
i = 0
for p, c, s in zip(parcs, chains, scores):
    if i%100==0:
        print(i)
    i+=1
    pdb_seq = get_pdb_seq(c)
    if len(pdb_seq) < 25:
        short += 1
        continue
    elif all([e=='X' for e in pdb_seq]):
        allX += 1
        continue
    elif sum([e=='X' for e in pdb_seq]) > len(pdb_seq)*0.3:
        manyX += 1
        continue
    good.append( (p,c,s) )



print('tot', len(parcs))
print('too short', short)
print('all unknown', allX)
print('30% < unk < 100%', manyX)
print('leftover', len(parcs)-short-allX-manyX)


import pickle
with open('tmp.out.pkl','wb') as f:
    pickle.dump(good, f)



good_scores = [t[2] for t in good]
print('new mean',np.mean(good_scores))








