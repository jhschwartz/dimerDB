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




# TODO: switch tests of calc_nwalign to match glocal USalign NWalign implementation

#def calc_nwalign(nw, fasta1: str, fasta2: str) -> float:
#    exe = f'{nw} {fasta1} {fasta2}'
#    result = subprocess.run(exe, text=True, capture_output=True, shell=True)
#    if result.stderr != '':
#        raise RuntimeError(f'nwalign for seqs {fasta1} and {fasta2} failed with stderr: {result.stderr}')
#    # The line we want looks like "Sequence identity:    0.625 (=   5/   8)"
#    identity = float(re.findall(r'Sequence identity:\s+([\d\.]+)', result.stdout)[0])
#    return identity
#
#
#@contextlib.contextmanager
#def fasta_of_pdb(pdbfile, pdb2fasta):
#    with tempfile.TemporaryDirectory() as td:
#        exe = f'{pdb2fasta} {pdbfile} > {td}/fasta'
#        result = subprocess.run(exe, text=True, capture_output=True, shell=True, check=True)
#        if result.stderr != '':
#            raise RuntimeError(f'pdb2fasta for pdb file {pdbfile} failed with stderr: {result.stderr}')
#        yield f'{td}/fasta'
#
#
#def nw_fasta_to_pdb(fastafile, pdbfile, nw, pdb2fasta):
#    with fasta_of_pdb(pdbfile, pdb2fasta) as fasta2:
#        return calc_nwalign(nw, fastafile, fasta2)
#
#
