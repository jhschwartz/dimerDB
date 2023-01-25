import yaml

import sys
sys.path.append('../scripts')
from check_chains_contact import check_chains_contact
from name_pdb import dimer2pdbs

labels_loc = '../bin/labels'
sys.path.append(labels_loc)
import labels as label_funcs

inyaml = '/nfs/turbo/umms-petefred/jaschwa/dimerDB/intermediates/homodimer_filtering/9B/UPI000013389B/initial.yaml'
name = 'UPI000013389B'

lib = '../lib'
dist_thresh = 8
count_thresh = 10

with open(inyaml) as f:
    data = yaml.safe_load(f)

structure_pairs = data[name] 

print('total chain pairing possibilities:', len(structure_pairs))

count = 0
for sp in structure_pairs:
    pdb1, pdb2 = dimer2pdbs(sp, lib)
    if check_chains_contact(pdb1, pdb2, label_funcs, dist_thresh, count_thresh):
        print(sp) 
        count += 1

print('resulting pairs in contact:', count)

