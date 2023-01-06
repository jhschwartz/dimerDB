import unittest
import subprocess
import tempfile
import gzip

exe = "../parallel_ungzip_all.sh"

class TestParallelUngzip(unittest.TestCase):
    def test_parallel_ungzip_divisible(self):
        with tempfile.TemporaryDirectory() as td:
            for i in range(10):
                path = f'{td}/{str(i)}.gz'
                with gzip.open(path, 'wt') as f:
                    f.write('Hello world!')
            subprocess.run(f'{exe} {td} 2 2', shell=True, check=True)
            for i in range(10):
                path = f'{td}/{str(i)}'
                with open(path, 'r') as f:
                    self.assertTrue('Hello world!' in f.read())

    def test_parallel_ungzip_indivisible(self):
        with tempfile.TemporaryDirectory() as td:
            for i in range(10):
                path = f'{td}/{str(i)}.gz'
                with gzip.open(path, 'wt') as f:
                    f.write('Hello world!')
            subprocess.run(f'{exe} {td} 2 3', shell=True, check=True)
            for i in range(10):
                path = f'{td}/{str(i)}'
                with open(path, 'r') as f:
                    self.assertTrue('Hello world!' in f.read())
                
if __name__ == '__main__':
    unittest.main()
