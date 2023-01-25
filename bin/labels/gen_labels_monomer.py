#!/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python

import sys
import pickle
import labels as label_funcs

if __name__ == '__main__':
    if not len(sys.argv) == 4:
        print('Expected use: ./gen_labels.py input1.pdb seq.txt output_filename')
        print('The output is a pickled dictionary with keys "intrachain_labels" (dict), "L" (int), "seq" (str)')
        print('The three dicts in the output dict each contain several label maps "omega", "theta", "phi", "dis" (cb), "ca_dis".')
        print('These maps are zero-indexed, so if you want residue 12 you will have to look at index 11, for example.')
        exit(-1)
    
    chain1_fn = sys.argv[1]
    seq_fn = sys.argv[2]
    outfile = sys.argv[3]
        
    # checks on args
    if not chain1_fn.endswith('.pdb'):
        raise ValueError('expected input file to end in ".pdb" but it does not.')
        
    # read seq file
    with open(seq_fn, 'r') as seq_fs:
        seq_lines = seq_fs.readlines()

        if len(seq_lines) == 2 and seq_lines[0].startswith('>'):
            seq = seq_lines[1].strip()
        elif len(seq_lines) == 1:
            seq = seq_lines[0].strip()
        else:
            raise ValueError('seq is not valid.')

    
    # read chains 1 and 2
    res_dic_1, nums_1, pdb_seq = label_funcs.read_chain(chain1_fn)
  
    if seq != pdb_seq:
        raise ValueError('the seq read from pdb 1 does not match the seq.txt for chain1')

    L = len(seq)
    intrachain_labels = label_funcs.extract_labels(res_dic_1, nums_1)

    data = {
        'intrachain_labels': intrachain_labels,
        'L': L,
        'seq': seq,
    }
    
    with open(outfile, 'wb') as f:
        pickle.dump(data, f)

