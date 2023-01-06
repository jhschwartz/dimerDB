import yaml
import time


def check_uniparc2others(yamlfile):
    with open(yamlfile, 'r') as f:
        data = yaml.safe_load(f)

    assert not 'UPI00000BE9B0' in data

    entry = data['UPI000002DB1C']
    assert 'P05067' in entry['uniprot']
    assert 'A0A140VJC8' in entry['uniprot']
    assert '1amb_A' in entry['pdb'] 
    assert '1aap_B' in entry['pdb'] 
    print('ok')

def check_homodimers(yamlfile):
    with open(yamlfile, 'r') as f:
        data = yaml.safe_load(f)
    
    parc = 'UPI000016A4A9' # human aldolase
    assert parc in data
    assert '1xdl_A-1xdl_Y' in data[parc]
    assert not '1xdl_A-1xdm_B' in data[parc]
    print('ok')


def check_heterodimers(yamlfile):
    with open(yamlfile, 'r') as f:
        data = yaml.safe_load(f)
    
    parc = 'UPI000003EB27-UPI000013471A'
    pdb = '4jk1_E-4jk1_D'
    assert pdb in data[parc]
    print('ok')
   

if __name__ == '__main__':
    check_uniparc2others('../intermediates/uniparc2others.yaml')
    check_homodimers('../intermediates/all_homodimers.yaml')
    check_heterodimers('../intermediates/all_heterodimers.yaml')
