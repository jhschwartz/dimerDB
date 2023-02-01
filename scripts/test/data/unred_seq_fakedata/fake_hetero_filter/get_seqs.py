import yaml
import subprocess
import os 

with open('../heterodimers.yaml') as f:
    data = yaml.safe_load(f)


dimers = ['UPI0000052DE7-UPI0000052FA5',\
    'UPI0006C47D63-UPI0006BF9280',\
    'UPI0000052DE7-UPI0000052FA5',\
    'UPI0006BF9280-UPI0006C47D63',\
    'UPI000012DB37-UPI00003BB37D',\
    'UPI000012E856-UPI0003C96339',\
    'UPI000209B964-UPI0003CD15F4',\
    'UPI0003C96339-UPI000012E856']


for pair in dimers:
    n1, n2 = pair.split('-')
    div = n1[-1] + '-' + n2[-1]
    dir_ = f'{div}/{pair}'
    os.makedirs(dir_, exist_ok=True)
    out = f'{div}/{pair}/seq1.fasta'
    subprocess.Popen(f'wget -O {out} https://rest.uniprot.org/uniparc/{n1}.fasta'.split(), stdout=None, stderr=None)
    out = f'{div}/{pair}/seq2.fasta'
    subprocess.Popen(f'wget -O {out} https://rest.uniprot.org/uniparc/{n2}.fasta'.split(), stdout=None, stderr=None)

