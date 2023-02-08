import yaml

import sys
sys.path.append('../scripts')
from check_chains_contact import check_chains_contact
from name_pdb import dimer2pdbs

labels_loc = '../bin/labels'
sys.path.append(labels_loc)
import labels as label_funcs

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-u', '--uniparc', required=True) 
parser.add_argument('-o', '--outfile', required=True) 
args = parser.parse_args()

name = args.uniparc
div = name[-2:]
inyaml = f'/nfs/turbo/umms-petefred/jaschwa/dimerDB/intermediates/homodimer_filtering/{div}/{name}/initial.yaml'

lib = '../lib'
dist_thresh = 8
count_thresh = 10

with open(inyaml) as f:
    data = yaml.safe_load(f)

structure_pairs = data[name] 

print('total chain pairing possibilities:', len(structure_pairs))

good = []
for sp in structure_pairs:
    pdb1, pdb2 = dimer2pdbs(sp, lib)
    print(sp, end='    ')
    if check_chains_contact(pdb1, pdb2, label_funcs, dist_thresh, count_thresh):
        good.append(sp)
        print('O')
    else:
        print('X')

print('-----')

with open(args.outfile, 'w') as f:
    for g in good:
        f.write(f'{g}\n')

print('resulting pairs in contact:', len(good))

