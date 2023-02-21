import subprocess
import contextlib
import tempfile
import re

def calc_nwalign(nw, fasta1: str, fasta2: str) -> float:
    exe = f'{nw} {fasta1} {fasta2}'
    result = subprocess.run(exe, text=True, capture_output=True, shell=True)
    if result.stderr != '':
        raise RuntimeError(f'nwalign for seqs {fasta1} and {fasta2} failed with stderr: {result.stderr}')
    # The line we want looks like "Sequence identity:    0.625 (=   5/   8)"
    identity = float(re.findall(r'Sequence identity:\s+([\d\.]+)', result.stdout)[0])
    return identity


@contextlib.contextmanager
def fasta_of_pdb(pdbfile, pdb2fasta):
    with tempfile.TemporaryDirectory() as td:
        exe = f'{pdb2fasta} {pdbfile} > {td}/fasta'
        result = subprocess.run(exe, text=True, capture_output=True, shell=True, check=True)
        if result.stderr != '':
            raise RuntimeError(f'pdb2fasta for pdb file {pdbfile} failed with stderr: {result.stderr}')
        yield f'{td}/fasta'


def nw_fasta_to_pdb(fastafile, pdbfile, nw, pdb2fasta):
    with fasta_of_pdb(pdbfile, pdb2fasta) as fasta2:
        return calc_nwalign(nw, fastafile, fasta2)


