'''
download_uniparc_fasta.py - downloads from uniprot the fasta file of a given uniparc id

Written by Jacob Schwartz (jaschwa@umich.edu) in December 2022.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

This function is unittested by test/test_download_uniparc_fasta.py and passing as of 1/12/2023.
This work requires python >= 3.8
'''
import requests

def download_fasta(uniparc, out_path):
    url = f'https://rest.uniprot.org/uniparc/{uniparc}.fasta'
    r = requests.get(url)
    if r.status_code != 200:
        raise ConnectionError(f'received status code {r.status_code} upon requesting the fasta of {uniparc} from {url}')
    with open(out_path, 'w') as f:
        f.write(r.text)

