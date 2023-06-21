import sys
import itertools

dimers_list_path = sys.argv[1]

dimers = []

with open(dimers_list_path, 'r') as f:
    for line in f:
        dimer = line.strip()
        div = dimer[1:3]
        dimer = f'{div}/{dimer}'
        dimers.append(dimer)

for d1, d2 in itertools.combinations(dimers, 2):
    print(d1, d2)
