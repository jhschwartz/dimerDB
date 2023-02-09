import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i','--infile')
parser.add_argument('-o','--outfile')
args = parser.parse_args()

with open(args.infile, 'r') as fi, open(args.outfile, 'w') as fo:
    while line := fi.readline():
        s = line
        pattern = '.{4}_.{1}_.{1}'
        pattern2 = '.{4}_.{1}_.{1}'
        if re.search(pattern, line):
            for fname in re.findall(r'(.{4}_.{1}_.{1})', line):
                pdb, chain, model = fname.split('_')
                nname = f'{pdb}_a1_m{model}_c{chain}'
                line = line.replace(fname, nname)
        elif re.search(pattern2, line):
            for fname in re.findall(r'(.{4}-.{1}-.{1})', line):
                pdb, chain, model = fname.split('-')
                nname = f'{pdb}-a1-m{model}-c{chain}'
                line = line.replace(fname, nname)
        #if (s!=line): print(s,line)
        fo.write(line) 
