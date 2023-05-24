import sys
sys.path.append('..')

import unittest
import os
import tempfile
import filecmp
import shutil
from tmscore_database import update_db

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
db = f'{test_dir}/data/tmscore_database/db'
db_update = f'{test_dir}/data/tmscore_database/db_for_update'


class TestUpdateTMscores(unittest.TestCase):
    def test_init_database_aa_only(self):
        with tempfile.TemporaryDirectory() as tmp_db:
            data = [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ', 
                            '0.56', '0.77'), 
                     ('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB',
                            '0.84', '0.59'),
                     ('5aa9-a14-m8-cU7H_5aa9-a14-m19-c47', '8aau-a8-m1-cQ_8aau-a8-m1-cQ',
                            '0.00', '1.00')
                ]

            update_db(new_pairs_scores=data, db_path=tmp_db)
            result_file = os.path.join(tmp_db, 'aa.tsv')

            expected_file = os.path.join(db, 'aa.tsv')

            self.assertTrue(filecmp.cmp(result_file, expected_file, shallow=False))


    def test_init_database_aa_bb(self):
        with tempfile.TemporaryDirectory() as tmp_db:
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

            update_db(new_pairs_scores=data, db_path=tmp_db)
            result_file_aa = os.path.join(tmp_db, 'aa.tsv')
            result_file_bb = os.path.join(tmp_db, 'bb.tsv')

            expected_file_aa = os.path.join(db, 'aa.tsv')
            expected_file_bb = os.path.join(db, 'bb.tsv')

            self.assertTrue(filecmp.cmp(result_file_aa, expected_file_aa, shallow=False))
            self.assertTrue(filecmp.cmp(result_file_bb, expected_file_bb, shallow=False))


    def test_insert_scores_aa_only(self):
        with tempfile.TemporaryDirectory() as tmp_db:
            shutil.copy(os.path.join(db_update, 'aa.tsv'), tmp_db)
            
            aafile = os.path.join(tmp_db, 'aa.tsv')
            expected_file = os.path.join(db, 'aa.tsv')
            self.assertFalse(filecmp.cmp(aafile, expected_file))

            new_data = [('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB',
                            '0.84', '0.59')]

            update_db(new_pairs_scores=new_data, db_path=tmp_db)

            self.assertTrue(filecmp.cmp(aafile, expected_file))


    def test_insert_scores_aa_bb(self):
        with tempfile.TemporaryDirectory() as tmp_db:
            shutil.copy(os.path.join(db_update, 'aa.tsv'), tmp_db)
            shutil.copy(os.path.join(db_update, 'bb.tsv'), tmp_db)
            
            aafile = os.path.join(tmp_db, 'aa.tsv')
            bbfile = os.path.join(tmp_db, 'bb.tsv')
            expected_aa = os.path.join(db, 'aa.tsv')
            expected_bb = os.path.join(db, 'bb.tsv')
            self.assertFalse(filecmp.cmp(aafile, expected_aa))
            self.assertFalse(filecmp.cmp(bbfile, expected_bb))

            new_data = [('1aa2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB',
                            '0.84', '0.59'),
                        ('1bb5-a1-m8-cA_1bb5-a1-m9-cA', '1bb5-a1-m8-cB_1bb5-a1-m9-cB',
                            '0.99', '0.01')
                    ]

            update_db(new_pairs_scores=new_data, db_path=tmp_db)

            self.assertTrue(filecmp.cmp(aafile, expected_aa))
            self.assertTrue(filecmp.cmp(bbfile, expected_bb))


    def test_change_existing_score(self):
        with tempfile.TemporaryDirectory() as tmp_db:
            shutil.copy(os.path.join(db_update, 'aa.tsv'), tmp_db)
            
            change_data = [('1aa2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ', 
                            '0.66', '0.11')]
            
            aafile = os.path.join(tmp_db, 'aa.tsv')
            expected_first_line_items = list(change_data[0])

            with open(aafile, 'r') as f:
                firstlineitems = f.readline().split()
                self.assertNotEqual(firstlineitems, expected_first_line_items)

            update_db(new_pairs_scores=change_data, db_path=tmp_db)

            with open(aafile, 'r') as f:
                firstlineitems = f.readline().split()
                self.assertEqual(firstlineitems, expected_first_line_items)

