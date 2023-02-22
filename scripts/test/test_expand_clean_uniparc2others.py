import unittest
import sys
import pickle
import tempfile

sys.path.append('..')
from expand_clean_uniparc2others import _expand_chains_across_assemblies_models, \
    _compare_uniparc_to_chain, expand_clean_uniparc2others

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
test_data = f'{test_dir}/data/expand_clean_uniparc2others'
test_lib = f'{test_data}/lib'


config = {
    'database_settings': {
        'chain_min_seq_len': 30
    },
    'paths': {
        'uniparc_seqs': f'{test_lib}/uniparc',
        'lib': test_lib 
    }
}


class TestFilterUnitparc2othersSeqmatch(unittest.TestCase):
    def test_expand_chains_across_assemblies_models(self):
        chains = ['1on3_A', '1on3_B', '6k3g_B']
        with open(f'{test_lib}/rcsb_index.pkl', 'rb') as f:
            index = pickle.load(f)
        result = _expand_chains_across_assemblies_models(chains, index)
        expected = ['1on3_a1_m1_cA', '1on3_a1_m1_cB', '6k3g_a1_m1_cB', '6k3g_a1_m2_cB', '6k3g_a2_m2_cB']
        self.assertEqual(set(result), set(expected))


    def test_compare_uniparc_to_chain_identical(self):
        uniparc_id = 'UPI00000C0390'
        chain_name = '1on3_a1_m1_cA'
        result = _compare_uniparc_to_chain(uniparc_id=uniparc_id, chain_name=chain_name, config=config)
        expected = 1.0
        self.assertEqual(result, expected)

    
    def test_compare_uniparc_to_chain_notidentical(self):
        uniparc_id = 'UPI000013594A'
        chain_name = '1on3_a1_m1_cA'
        result = _compare_uniparc_to_chain(uniparc_id=uniparc_id, chain_name=chain_name, config=config)
        expected = 0.161
        self.assertEqual(result, expected)


    def test_expand_clean_uniparc2others(self):
        inyaml = f'{test_data}/uniparc2others.yaml'
        with tempfile.TemporaryDirectory() as td:
            outyaml = f'{td}/out.yaml'
            expand_clean_uniparc2others(inyaml=inyaml, outyaml=outyaml, config=config)
            with open(outyaml, 'r') as f:
                result = yaml.safe_load(f)
        expected = {
            'UPI0000002FB4': ['1cv7_A'],
            'UPI00000000C1': ['3ek7_A', '1a29_A']
        }
        self.assertEqual(sorted(expected), sorted(result))
