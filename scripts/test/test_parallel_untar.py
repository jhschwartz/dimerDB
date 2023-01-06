import unittest
import subprocess
import tempfile
import glob 
import os
import shutil

exe_tar = os.path.realpath('../parallel_untar_all.sh')
exe_gzip = os.path.realpath('../parallel_ungzip_all.sh')

        


class TestParallelUntar(unittest.TestCase):
    def test_parallel_untar_divisible(self):
        with tempfile.TemporaryDirectory() as td:
            os.chdir(td)

            names = [str(i) for i in range(10)]
            
            # 1: create directories with files in them and tar each dir
            for dir_ in names:
                os.mkdir(dir_)
                file_ = f'{dir_}/file.txt'
                with open(file_, 'w') as f:
                    f.write('Hello world!')
                subprocess.run(f'tar -czf {dir_}.tar.gz {dir_}', shell=True, check=True)
            tars = glob.glob(f'*.tar.gz')
            self.assertEqual(len(tars), len(names))

            # 2: delete the source directories, leaving only the tars 
            for dir_ in names:  
                shutil.rmtree(dir_)
            
            # 3: use script to extract everything and check
            subprocess.run(f'{exe_tar} . 2 2', shell=True, check=True)
            for dir_ in names:
                file_ = f'{dir_}/file.txt'
                with open(file_, 'r') as f:
                    self.assertTrue('Hello world!' in f.read())

            

            

    def test_parallel_untar_indivisible(self):
        with tempfile.TemporaryDirectory() as td:
            os.chdir(td)

            names = [str(i) for i in range(100)]
            
            # 1: create directories with files in them and tar each dir
            for dir_ in names:
                os.mkdir(dir_)
                file_ = f'{dir_}/file.txt'
                with open(file_, 'w') as f:
                    f.write('Hello world!')
                subprocess.run(f'tar -czf {dir_}.tar.gz {dir_}', shell=True, check=True)
            tars = glob.glob(f'*.tar.gz')
            self.assertEqual(len(tars), len(names))

            # 2: delete the source directories, leaving only the tars 
            for dir_ in names:  
                shutil.rmtree(dir_)
            
            # 3: use script to extract everything and check
            subprocess.run(f'{exe_tar} . 2 13', shell=True, check=True)
            for dir_ in names:
                file_ = f'{dir_}/file.txt'
                with open(file_, 'r') as f:
                    self.assertTrue('Hello world!' in f.read())


if __name__ == '__main__':
    unittest.main()
