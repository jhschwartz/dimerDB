import sys
sys.path.append('..')

import unittest
import os
import tempfile
import filecmp
import shutil
from tmscore_database import TMDB

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
db_expected = f'{test_dir}/data/tmscore_database/db'
db_update = f'{test_dir}/data/tmscore_database/db_for_update'


class TestUpdateTMscores(unittest.TestCase):
    def test_init_database_aa_only(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            data = [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ', 
                            '0.56', '0.77'), 
                     ('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB',
                            '0.84', '0.59'),
                     ('5aa9-a14-m8-cU7H_5aa9-a14-m19-c47', '8aau-a8-m1-cQ_8aau-a8-m1-cQ',
                            '0.00', '1.00')
                ]

            db = TMDB(tmp_db_path)
            db.update(data)

            result_file = os.path.join(tmp_db_path, 'aa.tsv')

            expected_file = os.path.join(db_expected, 'aa.tsv')

            self.assertTrue(filecmp.cmp(result_file, expected_file, shallow=False))


    def test_init_database_aa_bb(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            data = [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ', 
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
                            '0.04', '0.23')
                ]

            db = TMDB(tmp_db_path)
            db.update(data)

            result_file_aa = os.path.join(tmp_db_path, 'aa.tsv')
            result_file_bb = os.path.join(tmp_db_path, 'bb.tsv')

            expected_file_aa = os.path.join(db_expected, 'aa.tsv')
            expected_file_bb = os.path.join(db_expected, 'bb.tsv')

            self.assertTrue(filecmp.cmp(result_file_aa, expected_file_aa, shallow=False))
            self.assertTrue(filecmp.cmp(result_file_bb, expected_file_bb, shallow=False))


    def test_insert_scores_aa_only(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            shutil.copy(os.path.join(db_update, 'aa.tsv'), tmp_db_path)
            
            aafile = os.path.join(tmp_db_path, 'aa.tsv')
            expected_file = os.path.join(db_expected, 'aa.tsv')
            self.assertFalse(filecmp.cmp(aafile, expected_file))

            new_data = [('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB',
                            '0.84', '0.59')]

            db = TMDB(tmp_db_path)
            db.update(new_data)

            self.assertTrue(filecmp.cmp(aafile, expected_file))


    def test_insert_scores_aa_bb(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            shutil.copy(os.path.join(db_update, 'aa.tsv'), tmp_db_path)
            shutil.copy(os.path.join(db_update, 'bb.tsv'), tmp_db_path)
            
            aafile = os.path.join(tmp_db_path, 'aa.tsv')
            bbfile = os.path.join(tmp_db_path, 'bb.tsv')
            expected_aa = os.path.join(db_expected, 'aa.tsv')
            expected_bb = os.path.join(db_expected, 'bb.tsv')
            self.assertFalse(filecmp.cmp(aafile, expected_aa))
            self.assertFalse(filecmp.cmp(bbfile, expected_bb))

            new_data = [('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB',
                            '0.84', '0.59'),
                        ('1bb5-a1-m8-cA_1bb5-a1-m9-cA', '1bb5-a1-m8-cB_1bb5-a1-m9-cB',
                            '0.99', '0.01')
                    ]

            db = TMDB(tmp_db_path)
            db.update(new_data)

            self.assertTrue(filecmp.cmp(aafile, expected_aa))
            self.assertTrue(filecmp.cmp(bbfile, expected_bb))


    def test_change_existing_score(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            shutil.copy(os.path.join(db_update, 'aa.tsv'), tmp_db_path)
            
            change_data = [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ', 
                            '0.66', '0.11')]
            
            aafile = os.path.join(tmp_db_path, 'aa.tsv')
            expected_first_line_items = list(change_data[0])

            with open(aafile, 'r') as f:
                firstlineitems = f.readline().split()
                self.assertNotEqual(firstlineitems, expected_first_line_items)

            db = TMDB(tmp_db_path)
            db.update(change_data)

            with open(aafile, 'r') as f:
                firstlineitems = f.readline().split()
                self.assertEqual(firstlineitems, expected_first_line_items)



class TestGetTMscores(unittest.TestCase):
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
       
        db = TMDB(db_expected)
        result = db.get(dimers)

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
       
        db = TMDB(db_expected)
        result = db.get(dimers)

        self.assertEqual(sorted(list(result)), sorted(expected))

        
    def test_none_found(self):
        dimers = [('1zz2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ'), 
                  ('1zz2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB'),
                  ('1zz9-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB')
            ]

        expected = []

        db = TMDB(db_expected)
        result = db.get(dimers)

        self.assertEqual(list(result), expected)


