#!/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python
import numpy as np

import sys
sys.path.append('../scripts')
from unredundant import RedundantDimerStructures

config = {
    'paths': {
        'lib': '../lib',
        'mmalign_exe': '../bin/MMalign'
    }
}

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile') 
parser.add_argument('-c', '--cores')
parser.add_argument('-l', '--left')
parser.add_argument('-r', '--right')
args = parser.parse_args()

with open(args.infile) as f:
    all_dimers = [l.strip() for l in f.readlines()]

thresh = 0.5
structures = RedundantDimerStructures(all_dimers, thresh, config)
print(f'{str(len(structures.things))} dimers loaded')

structures.initiate_distance_matrix(num_workers=int(args.cores))

for thresh in np.linspace(float(args.left), float(args.right), 50):
    structures.threshold = thresh
    num_clusters = structures.initiate_clusters()
    print('threshold:', thresh, ',', 'clusters:', num_clusters)

