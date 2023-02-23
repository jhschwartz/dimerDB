import sys
sys.path.append('../scripts')

from derive_all_possible_homodimers import homodimers

file_in = 'tmpsmall_uniparc2others.pkl' 
file_out =  'out.yaml'
lib = '../lib'

@profile
def main():
    homodimers(file_in, file_out, lib)
   
main()
