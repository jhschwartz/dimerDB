import sys
sys.path.append('..')

import unittest
from tmscore_database import lookup_scores

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
db = f'{test_dir}/data/tmscore_database/db'

class TestLookupTMscores(unittest.TestCase):
    def test_aa_only(self):
        dimers = [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ'), 
                  ('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB'),
                  ('1aa9-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB'),
                  ('5aa9-a14-m8-cU7H_5aa9-a14-m19-c47', '8aau-a8-m1-cQ_8aau-a8-m1-cQ')
                ]

        expected =  [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ', 
                            '0.56', '0.77'), 
                     ('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB',
                            '0.84', '0.59'),
                     ('5aa9-a14-m8-cU7H_5aa9-a14-m19-c47', '8aau-a8-m1-cQ_8aau-a8-m1-cQ',
                            '0.00', '1.00')
                ]
       
        result = lookup_scores(dimer_pairs=dimers, db_path=db)

        self.assertEqual(sorted(list(result)), sorted(expected))


    
    def test_aa_bb(self):
        dimers = [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ'), 
                  ('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB'),
                  ('1aa9-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB'),
                  ('5aa9-a14-m8-cU7H_5aa9-a14-m19-c47', '8aau-a8-m1-cQ_8aau-a8-m1-cQ'),
                  ('1bb5-a1-m8-cA_1bb5-a1-m9-cA', '1bb5-a1-m8-cB_1bb5-a1-m9-cB'),
                  ('1bb5-a1-m8-cA_1bb5-a1-m9-cA', '2qq6-a1-m1-cA_2qq6-a1-m1-cB'),
                  ('2bb6-a2-m1-cF_2bb6-a2-m1-cG', '4aa9-a1-m1-cA_4aa9-a1-m1-cB'),
                  ('3bb7-a1-m1-cN_3bb7-a1-m1-cM', '8bb1-a1-m1-cA_8bb1-a1-m1-cB'),
                  ('3bb7-a1-m1-cN_3bb7-a1-m1-cM', '8bb1-a1-m1-cC_8bb1-a1-m1-cD')
                ]

        expected =  [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ', 
                            '0.56', '0.77'), 
                     ('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB',
                            '0.84', '0.59'),
                     ('5aa9-a14-m8-cU7H_5aa9-a14-m19-c47', '8aau-a8-m1-cQ_8aau-a8-m1-cQ',
                            '0.00', '1.00'),

                     ('1bb5-a1-m8-cA_1bb5-a1-m9-cA', '1bb5-a1-m8-cB_1bb5-a1-m9-cB', 
                            '0.99', '0.01'), 
                     ('2bb6-a2-m1-cF_2bb6-a2-m1-cG', '4aa9-a1-m1-cA_4aa9-a1-m1-cB',
                            '0.50', '0.50'),
                     ('3bb7-a1-m1-cN_3bb7-a1-m1-cM', '8bb1-a1-m1-cA_8bb1-a1-m1-cB',
                            '0.04', '0.23'),

                ]
       
        result = lookup_scores(dimer_pairs=dimers, db_path=db)

        self.assertEqual(sorted(list(result)), sorted(expected))

        
    def test_none_found(self):
        dimers = [('1zz2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ'), 
                  ('1zz2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB'),
                  ('1zz9-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB')
            ]

        expected = []

        result = lookup_scores(dimer_pairs=dimers, db_path=db)

        self.assertEqual(list(result), expected)

