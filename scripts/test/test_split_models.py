import unittest
import subprocess
import tempfile
import shutil

exe = "../split_models.sh"
ref_assembly = ''
ref_model1 = ''
ref_model2 = ''
ref_name = ''

class TestSplitModels(unittest.TestCase):
    def test_split_2gqq_into_models(self):
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(ref_assembly, td)
            subprocess.run(f'{exe} {ref_assembly} {ref_name} {exe}', shell=True, check=True)
            file1_result = f'{td}/{ref_name}-model1.pdb'
            file2_result = f'{td}/{ref_name}-model2.pdb'
            model1_correct = ??
            model2_correct = ??
            self.assertTrue(model1_correct)
            self.assertTrue(model2_correct)
               

if __name__ == '__main__':
    unittest.main()
