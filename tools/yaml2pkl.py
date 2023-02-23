import yaml
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--yaml', required=True)
parser.add_argument('-o', '--pkl', required=True)
args = parser.parse_args()

with open(args.yaml, 'r') as f:
    data = yaml.safe_load(f)

with open(args.pkl, 'wb') as f:
    pickle.dump(data, f)
