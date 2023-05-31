import subprocess
from operator import xor

def calc_nwalign_glocal(USnw, pdb1, pdb2, refshorter=True, reflonger=False):
    if not xor(refshorter, reflonger):
        raise ValueError('you must choose exactly one of refshorter and reflonger to be True')
    exe = f'{USnw} -glocal 2 -het 1 -outfmt 2 {pdb1} {pdb2}'
    result = subprocess.run(exe, text=True, capture_output=True, shell=True)
    if result.stderr != '':
        raise RuntimeError(f'nwalign for seqs {pdb1} and {pdb2} failed with stderr: {result.stderr}')
    spl = result.stdout.split('\n')[1].split()
    L1, L2 = spl[5], spl[6]
    index_ID_refshorter = 3
    index_ID_reflonger = 2
    if int(L1) < int(L2):
        index_ID_refshorter = 2
        index_ID_reflonger = 3
    if refshorter:
        return float(spl[index_ID_refshorter])
    return float(spl[index_ID_reflonger])

