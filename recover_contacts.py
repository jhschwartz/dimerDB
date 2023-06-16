import os
import glob
import sys
sys.path.append('scripts')

from contact_database import ContactDB
from sortedcontainers import SortedSet
import name_pdb

contacts_db_path = 'lib/contactsdb'

db = ContactDB(contacts_db_path, max_buffer=100000)

c_tot = 0
p_tot = 0

for div_path in glob.iglob('intermediates/check_pairs/*'):
    div = div_path.split('/')[-1]

    pairsfile = os.path.join(div_path, 'intra_assembly_chain_pairs.txt')
    contactsfile = os.path.join(div_path, 'dimers.txt')

    if not os.path.exists(pairsfile) or not os.path.exists(contactsfile):
        print(f'no contacts for {div}')
        continue

    num_checks = 0

    is_contact = SortedSet()
    with open(contactsfile, 'r') as f:
        for line in f:
            is_contact.add(line.strip())
    with open(pairsfile, 'r') as f:
        for line in f:
            dimer = line.strip()
            c1, c2 = name_pdb.dimer2chains(dimer)
            contacting = dimer in is_contact
            db.auto_buffer( (c1, c2, contacting) )
            num_checks += 1
    
    print(f'{len(is_contact)} / {num_checks} pairs are contacts for {div}.')
    
    c_tot += len(is_contact)
    p_tot += num_checks


db.update_clear_buffer()

print(f'total of {c_tot} / {p_tot} contacts.')
    


