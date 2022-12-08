import itertools
from scripts.read_fasta import read_prot_from_fasta
import subprocess
import time
import numpy as np
import re

fasta = 'data/50seqs.fasta'
seqs = [s for h, s in read_prot_from_fasta(fasta)]
pairs = list(itertools.combinations(seqs, 2))

NW = '/nfs/amino-home/zhng/zhanglab_programs/NWalign/align'

Npair = len(pairs)
Nprot = len(seqs)

Lavg = np.mean([len(s) for s in seqs])
Lmax = np.max([len(s) for s in seqs])
Lmin = np.min([len(s) for s in seqs])
Lstd = np.std([len(s) for s in seqs])


t_0 = time.time()
for seq1, seq2 in pairs:
    result = subprocess.run(f'{NW} {seq1} {seq2} 3', capture_output=True, text=True, shell=True)
    if result.stderr != '':
        raise RuntimeError(f'NW did not work: {result.stderr}')
    seq_id = re.findall('Sequence identity\:\s+([\.\d]+)\s', result.stdout)[0]
    print(seq_id)

elapsed = time.time() - t_0
avg = elapsed / Npair

print(f'Ran {Npair} pairs in {elapsed} seconds total or {avg} seconds each.')
print(f'These {Nprot} sequences had an average length of {Lavg} (min = {Lmin}, max = {Lmax}, std = {Lstd}).')
