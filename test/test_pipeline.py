import unittest
import subprocess
from generate_test_lib import generate_test_lib

class TestPipeline(unittest.TestCase):

    
    def setUp(self):
        # set a tmpdir here...
        
   #     generate_test_lib(  source_lib='/nfs/turbo/umms-zcx/jaschwa/dimerDB-dev/lib',
  #                          test_entries_file='test_entries.txt', 
 #                           dest_lib='lib'  )

        #subprocess.run('./prep_test.sh', shell=True, check=True)

        subprocess.run('../run_test_smk.sh', shell=True, check=True, cwd='tmp')


    def tearDown(self):
        pass
        # cleanup would go here! 

        # ...and delete the tmpdir here


    def test_end_to_end(self):
        # some real tests would go here

        self.assertTrue(True)

