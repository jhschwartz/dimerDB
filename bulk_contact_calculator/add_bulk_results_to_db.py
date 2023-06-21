#!/expanse/lustre/projects/mia174/jaschwa/dimerDB/env/bin/python

import os
import glob
import sys
sys.path.append('../scripts')

from contact_database import ContactDB
from sortedcontainers import SortedSet
import name_pdb

contacts_db_path = '../lib/contactsdb'
bulk_results_dir = '../bulk_contacts_out_20230616'

db = ContactDB(contacts_db_path, max_buffer=100000)

c_tot = 0
p_tot = 0

infiles = sorted(glob.glob(f'{bulk_results_dir}/split_*'))
outfiles = sorted(glob.glob(f'{bulk_results_dir}/out_split_*'))

if len(infiles) != len(outfiles):
    raise ValueError

for i in range(len(infiles)):
    infile = infiles[i]
    outfile = outfiles[i]

    split_name_in = infile.split('_')[-1]
    split_name_out = outfile.split('_')[-1]
    if split_name_in != split_name_out:
        raise ValueError(f'split names not matching, in: {infile}, out: {outfile}')

    with open(infile, 'r') as fi, \
                    open(outfile, 'r') as fo:
        for inline, outline in zip(fi, fo):

            c1, c2 = inline.split()
            contacting = outline.split()[0] == '1'

            db.auto_buffer( (c1, c2, contacting) )

            p_tot += 1

            if contacting:
                c_tot += 1

        print(f'running total of {c_tot} / {p_tot} contacts.')


db.update_clear_buffer()


print(f'final total of {c_tot} / {p_tot} contacts.')

#for div_path in glob.iglob('intermediates/check_pairs/*'):
#    div = div_path.split('/')[-1]
#
#    pairsfile = os.path.join(div_path, 'intra_assembly_chain_pairs.txt')
#    contactsfile = os.path.join(div_path, 'dimers.txt')
#
#    if not os.path.exists(pairsfile) or not os.path.exists(contactsfile):
#        print(f'no contacts for {div}')
#        continue
#
#    num_checks = 0
#
#    is_contact = SortedSet()
#    with open(contactsfile, 'r') as f:
#        for line in f:
#            is_contact.add(line.strip())
#    with open(pairsfile, 'r') as f:
#        for line in f:
#            dimer = line.strip()
#            c1, c2 = name_pdb.dimer2chains(dimer)
#            contacting = dimer in is_contact
#            db.auto_buffer( (c1, c2, contacting) )
#            num_checks += 1
#    
#    print(f'{len(is_contact)} / {num_checks} pairs are contacts for {div}.')
#    
#    c_tot += len(is_contact)
#    p_tot += num_checks
#
#
#db.update_clear_buffer()
#
#print(f'total of {c_tot} / {p_tot} contacts.')
    


