import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--pkl', required=True)
args = parser.parse_args()

with open(args.pkl, 'rb') as f:
    print(pickle.load(f))
