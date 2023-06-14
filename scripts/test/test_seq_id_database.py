import sys
sys.path.append('..')

import unittest
import os
import tempfile
import filecmp
import shutil
from seq_id_database import SeqidDB

import pathlib
test_dir = pathlib.Path(__file__).parent.resolve()
db_expected = f'{test_dir}/data/seq_id_database/db'
db_update = f'{test_dir}/data/seq_id_database/db_for_update'



class TestUpdateSeqIDs(unittest.TestCase):
    def test_init_database_a_only(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            data = [ ('1az2-a1-m1-cA', '1az2-a1-m2-cA', 0.01),
                     ('5az9-a14-m8-cU7H', '5az9-a14-m19-c47', 0.56),
                     ('2aze-a1-m1-cA', '2aze-a1-m1-cB', 1.0),
                ]

            db = SeqidDB(tmp_db_path)
            db.update(data)

            result_file = os.path.join(tmp_db_path, 'a.tsv')
            expected_file = os.path.join(db_expected, 'a.tsv')

            self.assertTrue(filecmp.cmp(result_file, expected_file, shallow=False))


    def test_init_database_a_b(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            data = [ ('1az2-a1-m1-cA', '1az2-a1-m2-cA', 0.01),
                     ('2bz6-a2-m1-cF', '2bz6-a2-m1-cG', 0.301),
                     ('2aze-a1-m1-cA', '2aze-a1-m1-cB', 1.0),
                     ('5az9-a14-m8-cU7H', '5az9-a14-m19-c47', 0.56),
                     ('1bz5-a1-m8-cA', '1bz5-a1-m9-cA', 0.999),
                     ('3bz7-a1-m1-cN', '3bz7-a1-m1-cM', 0.526)
                ]

            db = SeqidDB(tmp_db_path)
            db.update(data)

            result_file_a = os.path.join(tmp_db_path, 'a.tsv')
            expected_file_a = os.path.join(db_expected, 'a.tsv')

            result_file_b = os.path.join(tmp_db_path, 'b.tsv')
            expected_file_b = os.path.join(db_expected, 'b.tsv')

            self.assertTrue(filecmp.cmp(result_file_a, expected_file_a, shallow=False))
            self.assertTrue(filecmp.cmp(result_file_b, expected_file_b, shallow=False))


    def test_update_database_a_only(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            shutil.copy(os.path.join(db_update, 'a.tsv'), tmp_db_path)
            
            a_file = os.path.join(tmp_db_path, 'a.tsv')
            expected_file = os.path.join(db_expected, 'a.tsv')
            self.assertFalse(filecmp.cmp(a_file, expected_file))

            new_data = [ ('2aze-a1-m1-cA', '2aze-a1-m1-cB', 1.0) ]

            db = SeqidDB(tmp_db_path)
            db.update(new_data)
            
            self.assertTrue(filecmp.cmp(a_file, expected_file))


    def test_update_database_a_b(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            shutil.copy(os.path.join(db_update, 'a.tsv'), tmp_db_path)
            shutil.copy(os.path.join(db_update, 'b.tsv'), tmp_db_path)
            
            a_file = os.path.join(tmp_db_path, 'a.tsv')
            a_expected_file = os.path.join(db_expected, 'a.tsv')
            self.assertFalse(filecmp.cmp(a_file, a_expected_file))

            b_file = os.path.join(tmp_db_path, 'b.tsv')
            b_expected_file = os.path.join(db_expected, 'b.tsv')
            self.assertFalse(filecmp.cmp(b_file, b_expected_file))
            
            new_data = [ ('3bz7-a1-m1-cN', '3bz7-a1-m1-cM', 0.526),
                         ('2aze-a1-m1-cA', '2aze-a1-m1-cB', 1.0) 
                    ]

            db = SeqidDB(tmp_db_path)
            db.update(new_data)

            self.assertTrue(filecmp.cmp(a_file, a_expected_file))
            self.assertTrue(filecmp.cmp(b_file, b_expected_file))


    def test_change_contact(self):
        with tempfile.TemporaryDirectory() as tmp_db_path:
            shutil.copy(os.path.join(db_update, 'a.tsv'), tmp_db_path)
            
            change_data = [ ('1az2-a1-m1-cA', '1az2-a1-m2-cA', 0.62), # change val
                            ('2aze-a1-m1-cA', '2aze-a1-m1-cB', 1.0)  # new entry
                    ]
            
            a_file = os.path.join(tmp_db_path, 'a.tsv')

            with open(a_file, 'r') as f:
                filedata = f.read()
                self.assertFalse('2aze' in filedata)
                self.assertFalse('1az2-a1-m1-cA\t1az2-a1-m2-cA\t0.62' in filedata)
                self.assertTrue('1az2-a1-m1-cA\t1az2-a1-m2-cA\t0.01' in filedata)

            db = SeqidDB(tmp_db_path)
            db.update(change_data)

            with open(a_file, 'r') as f:
                filedata = f.read()
                self.assertTrue('2aze-a1-m1-cA\t2aze-a1-m1-cB\t1.0' in filedata)
                self.assertTrue('1az2-a1-m1-cA\t1az2-a1-m2-cA\t0.62' in filedata)
                self.assertFalse('1az2-a1-m1-cA\t1az2-a1-m2-cA\t0.01' in filedata)




class TestGetSeqIDs(unittest.TestCase):
    def test_get_a_only(self):
        dimers = [ ('1az2-a1-m1-cA', '1az2-a1-m2-cA'),
                   ('2aze-a1-m1-cA', '2aze-a1-m1-cB'),
                   ('5az9-a14-m8-cU7H', '5az9-a14-m19-c47')
            ]

        expected = [ ('1az2-a1-m1-cA', '1az2-a1-m2-cA', 0.01),
                     ('2aze-a1-m1-cA', '2aze-a1-m1-cB', 1.0),
                     ('5az9-a14-m8-cU7H', '5az9-a14-m19-c47', 0.56)
            ]

        db = SeqidDB(db_expected)
        result = db.get(dimers)

        self.assertEqual(sorted(list(result)), sorted(expected))


    def test_get_a_b(self):
        dimers = [ ('1az2-a1-m1-cA', '1az2-a1-m2-cA'),
                   ('2bz6-a2-m1-cF', '2bz6-a2-m1-cG'),
                   ('1bz5-a1-m8-cA', '1bz5-a1-m9-cA'),
                   ('2aze-a1-m1-cA', '2aze-a1-m1-cB'),
                   ('3bz7-a1-m1-cN', '3bz7-a1-m1-cM'),
                   ('5az9-a14-m8-cU7H', '5az9-a14-m19-c47'),
            ]
        
        expected = [ ('1az2-a1-m1-cA', '1az2-a1-m2-cA', 0.01),
                     ('2aze-a1-m1-cA', '2aze-a1-m1-cB', 1.0),
                     ('5az9-a14-m8-cU7H', '5az9-a14-m19-c47', 0.56),
                     ('1bz5-a1-m8-cA', '1bz5-a1-m9-cA', 0.999),
                     ('2bz6-a2-m1-cF', '2bz6-a2-m1-cG', 0.301),
                     ('3bz7-a1-m1-cN', '3bz7-a1-m1-cM', 0.526)
            ]

        db = SeqidDB(db_expected)
        result = db.get(dimers)

        self.assertEqual(sorted(list(result)), sorted(expected))



    def test_get_nonefound(self):
        dimers = [('1zz2-a1-m1-cA_1aa2-a1-m2-cA', '1aa7-a1-m1-cF_1aa7-a1-m1-cJ'),
                  ('1zz2-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB'),
                  ('1zz9-a1-m1-cA_1aa2-a1-m3-cA', '2aae-a1-m1-cA_2aae-a1-m1-cB')
            ]

        expected = []

        db = SeqidDB(db_expected)
        result = db.get(dimers)

        self.assertEqual(list(result), expected)


