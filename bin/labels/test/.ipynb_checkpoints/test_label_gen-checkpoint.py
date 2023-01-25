#!/home/jaschwa/.miniconda3/envs/SLIPPI-label/bin/python

import sys
import os
import subprocess
import numpy as np
import pickle

sys.path.append('..')
from gen_labels import readpdb, extract_labels, extract_labels_dimer






def dist(a, b):
    return np.sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2 + (b[2]-a[2])**2)







###########################################################
### Test set 1: extract_labels as function
res_dic, nums, _ = readpdb('example.pdb', 'A')
assert len(nums) == 164 - 12 + 1
assert len([k for k in res_dic if 'CA' in k]) == 164 - 12 + 1

labels = extract_labels(res_dic, nums) 

cb_28 = (18.091, 115.601, 29.003)
cb_40 = (27.956, 131.375, 32.707)
assert dist(cb_28, cb_40) == labels['dis'][28-1, 40-1]
assert dist(cb_28, cb_40) == labels['dis'][40-1, 28-1]

ca_53 = (28.899, 122.081, 17.434)
ca_84 = (45.679, 90.680, 3.112)
assert dist(ca_53, ca_84) == labels['ca_dis'][53-1, 84-1]
assert dist(ca_53, ca_84) == labels['ca_dis'][84-1, 53-1]

# check no CB at a glycine
assert np.all(np.isnan(labels['dis'][56-1, :]))

print('Test set 1 PASS')
###########################################################












###########################################################
### Test set 2: extract_labels_dimer as function
res_dic_1, nums_1, _ = readpdb('example.pdb', 'A')
assert len(nums_1) == 164 - 12 + 1
assert len([k for k in res_dic_1 if 'CA' in k]) == 164 - 12 + 1

res_dic_2, nums_2, _ = readpdb('example.pdb', 'B')
assert len(nums_2) == 162 - 12 + 1
assert len([k for k in res_dic_2 if 'CA' in k]) == 162 - 12 + 1

intrachain_labels_1, intrachain_labels_2, interchain_labels, L1, L2 = extract_labels_dimer(res_dic_1, nums_1, res_dic_2, nums_2)    

# intrachain 1 tests 
cb_28 = (18.091, 115.601, 29.003)
cb_40 = (27.956, 131.375, 32.707)
assert dist(cb_28, cb_40) == intrachain_labels_1['dis'][28-1, 40-1]
assert dist(cb_28, cb_40) == intrachain_labels_1['dis'][40-1, 28-1]
ca_53 = (28.899, 122.081, 17.434)
ca_84 = (45.679, 90.680, 3.112)
assert dist(ca_53, ca_84) == intrachain_labels_1['ca_dis'][53-1, 84-1]
assert dist(ca_53, ca_84) == intrachain_labels_1['ca_dis'][84-1, 53-1]
assert np.all(np.isnan(intrachain_labels_1['dis'][56-1, :]))
assert intrachain_labels_1['omega'].shape[0] == L1


# intrachain 2 tests 
cb_113 = (42.559, 87.515, 28.350)
cb_121 = (45.823, 107.893, 26.774)
assert dist(cb_113, cb_121) == intrachain_labels_2['dis'][113-1, 121-1]
assert dist(cb_113, cb_121) == intrachain_labels_2['dis'][121-1, 113-1]
ca_80 = (44.065, 91.424, 31.130)
ca_14 = (5.248, 99.759, 21.153)
assert dist(ca_80, ca_14) == intrachain_labels_2['ca_dis'][80-1, 14-1]
assert dist(ca_80, ca_14) == intrachain_labels_2['ca_dis'][14-1, 80-1]
assert np.all(np.isnan(intrachain_labels_2['dis'][39-1, :]))
assert intrachain_labels_2['phi'].shape[0] == L2

# interchain tests
cb_chain1_res35 = (24.633, 128.795, 31.361)
cb_chain2_res25 = (18.169, 109.894, 16.517)
assert dist(cb_chain1_res35, cb_chain2_res25) == interchain_labels['dis'][35-1, 25-1]
ca_chain1_res118 = (31.732, 94.784, 19.962)
ca_chain2_res27 = (19.435, 114.034, 22.680)
assert dist(ca_chain1_res118, ca_chain2_res27) == interchain_labels['ca_dis'][118-1, 27-1]
assert np.all(np.isnan(interchain_labels['dis'][140-1, :])) # chain1 gly at 140
assert np.all(np.isnan(interchain_labels['dis'][:, 56-1])) # chain2 gly at 56
assert interchain_labels['theta'].shape[0] == L1
assert interchain_labels['theta'].shape[1] == L2

print('Test set 2 PASS')
###########################################################









###########################################################
### Test set 3: gen_labels.py as script (dimer calc)
subprocess.check_call('../gen_labels.py example.pdb A B testout.tmp', shell=True)

with open('testout.tmp' , 'rb') as f:
    data = pickle.load(f)

os.remove('testout.tmp')
    
intrachain_labels_1 = data['intrachain_labels_1']
intrachain_labels_2 = data['intrachain_labels_2']
interchain_labels = data['interchain_labels']
L1 = data['L1']
L2 = data['L2']

# intrachain 1 tests 
cb_28 = (18.091, 115.601, 29.003)
cb_40 = (27.956, 131.375, 32.707)
assert dist(cb_28, cb_40) == intrachain_labels_1['dis'][28-1, 40-1]
assert dist(cb_28, cb_40) == intrachain_labels_1['dis'][40-1, 28-1]
ca_53 = (28.899, 122.081, 17.434)
ca_84 = (45.679, 90.680, 3.112)
assert dist(ca_53, ca_84) == intrachain_labels_1['ca_dis'][53-1, 84-1]
assert dist(ca_53, ca_84) == intrachain_labels_1['ca_dis'][84-1, 53-1]
assert np.all(np.isnan(intrachain_labels_1['dis'][56-1, 12:]))
assert intrachain_labels_1['theta'].shape[0] == L1


# intrachain 2 tests 
cb_113 = (42.559, 87.515, 28.350)
cb_121 = (45.823, 107.893, 26.774)
assert dist(cb_113, cb_121) == intrachain_labels_2['dis'][113-1, 121-1]
assert dist(cb_113, cb_121) == intrachain_labels_2['dis'][121-1, 113-1]
ca_80 = (44.065, 91.424, 31.130)
ca_14 = (5.248, 99.759, 21.153)
assert dist(ca_80, ca_14) == intrachain_labels_2['ca_dis'][80-1, 14-1]
assert dist(ca_80, ca_14) == intrachain_labels_2['ca_dis'][14-1, 80-1]
assert np.all(np.isnan(intrachain_labels_2['dis'][39-1, :]))
assert intrachain_labels_2['omega'].shape[0] == L2

# interchain tests
cb_chain1_res35 = (24.633, 128.795, 31.361)
cb_chain2_res25 = (18.169, 109.894, 16.517)
assert dist(cb_chain1_res35, cb_chain2_res25) == interchain_labels['dis'][35-1, 25-1]
ca_chain1_res118 = (31.732, 94.784, 19.962)
ca_chain2_res27 = (19.435, 114.034, 22.680)
assert dist(ca_chain1_res118, ca_chain2_res27) == interchain_labels['ca_dis'][118-1, 27-1]
assert np.all(np.isnan(interchain_labels['dis'][140-1, :])) # chain1 gly at 140
assert np.all(np.isnan(interchain_labels['dis'][:, 56-1])) # chain2 gly at 56
assert interchain_labels['phi'].shape[0] == L1
assert interchain_labels['phi'].shape[1] == L2

print('Test set 3 PASS')
###########################################################




### END
print('ALL TESTS PASS')