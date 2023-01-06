import requests

def download_fasta(uniparc, out_path):
    url = f'https://rest.uniprot.org/uniparc/{uniparc}.fasta'
    r = requests.get(url)
    if r.status_code != 200:
        raise ConnectionError(f'received status code {r.status_code} upon requesting the fasta of {uniparc} from {url}')
    with open(out_path, 'w') as f:
        f.write(r.text)

