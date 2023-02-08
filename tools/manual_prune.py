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
parser.add_argument('-t', '--thresh')
args = parser.parse_args()

with open(args.infile) as f:
    all_dimers = [l.strip() for l in f.readlines()]

structures = RedundantDimerStructures(all_dimers, float(args.thresh), config)
print(f'{str(len(structures.things))} dimers loaded')

result = structures.prune_redundancy(num_workers=int(args.cores))
print(result)

