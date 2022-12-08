import yaml
import time

uniparc_yaml_fn = '../intermediates/uniparc2others.yaml'

t = time.time()
with open(uniparc_yaml_fn, 'r') as f:
    data = yaml.safe_load(f)

print('yaml loaded in ', time.time()-t, 'seconds')

assert not 'UPI00000BE9B0' in data

entry = data['UPI000002DB1C']
assert 'P05067' in entry['uniprot']
assert 'A0A140VJC8' in entry['uniprot']
assert '1AMB_A' in entry['pdb'] 
assert '1AAP_B' in entry['pdb'] 

print('success')
