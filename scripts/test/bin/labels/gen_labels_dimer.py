#!/nfs/turbo/umms-petefred/jaschwa/HDPRED/bin/python

import sys
import pickle
import labels as label_funcs

if __name__ == '__main__':
    if not len(sys.argv) == 6:
        print('Expected use: ./gen_labels.py input1.pdb seq1.txt input2.pdb seq2.txt output_filename')
        print('The output is a pickled dictionary with keys "intrachain_labels_1" (dict), "interchain_labels_2" (dict), "interchain_labels" (dict), "L1" (int), "L2" (int).')
        print('The three dicts in the output dict each contain several label maps "omega", "theta", "phi", "dis" (cb), "ca_dis".')
        print('These maps are zero-indexed, so if you want residue 12 you will have to look at index 11, for example.')
        exit(-1)
    
    chain1_fn = sys.argv[1]
    seq1_fn = sys.argv[2]
    chain2_fn = sys.argv[3]
    seq2_fn = sys.argv[4]
    outfile = sys.argv[5]
        
    # checks on args
    if not chain1_fn.endswith('.pdb') or not chain2_fn.endswith('.pdb'):
        raise ValueError('expected input files to end in ".pdb", but at least one does not.')
        
    # read seq files 1 and 2
    with open(seq1_fn, 'r') as seq1_fs, open(seq2_fn, 'r') as seq2_fs:
        seq1_lines = seq1_fs.readlines()
        seq2_lines = seq2_fs.readlines()

        def extract_seq(lines):
            if len(lines) == 2 and lines[0].startswith('>'):
                return lines[1].strip()
            elif len(lines) == 1:
                return lines[0].strip()
            else:
                raise ValueError('seq is not valid.')

        seq1 = extract_seq(seq1_lines)
        seq2 = extract_seq(seq2_lines)
    
    # read chains 1 and 2
    res_dic_1, nums_1, pdb_seq1 = label_funcs.read_chain(chain1_fn)
    res_dic_2, nums_2, pdb_seq2 = label_funcs.read_chain(chain2_fn)
  
    if seq1 != pdb_seq1:
        raise ValueError('the seq read from pdb 1 does not match the seq.txt for chain1')
    if seq2 != pdb_seq2:
        raise ValueError('the seq read from pdb 2 does not match the seq.txt for chain2')

    intrachain_labels_1, intrachain_labels_2, interchain_labels, L1, L2 = label_funcs.extract_labels_dimer(res_dic_1, nums_1, res_dic_2, nums_2)    
    
    data = {
        'intrachain_labels_1': intrachain_labels_1,
        'intrachain_labels_2': intrachain_labels_2,
        'interchain_labels': interchain_labels,
        'L1': L1,
        'L2': L2,
        'seq1': seq1,
        'seq2': seq2
    }
    
    with open(outfile, 'wb') as f:
        pickle.dump(data, f)

