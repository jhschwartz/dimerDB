'''
download_uniparc_fasta.py - downloads from uniprot the fasta file of a given uniparc id

Written by Jacob Schwartz (jaschwa@umich.edu) in December 2022.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu
'''
import requests
import os
from name_fasta import uniparc_fasta
from multiprocessing import Pool


def download_fasta(uniparc, out_path):
    url = f'https://rest.uniprot.org/uniparc/{uniparc}.fasta'
    r = requests.get(url)
    if r.status_code != 200:
        raise ConnectionError(f'received status code {r.status_code} upon requesting the fasta of {uniparc} from {url}')
    os.makedirs(os.path.dirname(os.path.realpath(out_path)), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(r.text)


def download_many_fasta(uniparcs, lib_path, num_workers=1):
    outpaths = [uniparc_fasta(uniparc_id=u, lib_path=lib_path, allow_nonexist=True) for u in uniparcs]
    with Pool(processes=num_workers) as p:
        p.starmap(download_fasta, zip(uniparcs, outpaths))
    return outpaths

