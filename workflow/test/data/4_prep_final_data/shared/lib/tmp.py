import subprocess
import sys
#sys.path.append('/nfs/turbo/umms-petefred/jaschwa/dimerDB-optimize-dev/bin/check_contact')
#from check_contact import count_contact_many


dimers  = [   \
('rcsb/ql/6qle-a1-m1-cH.pdb', 'rcsb/ql/6qle-a1-m1-cI.pdb'), \
('rcsb/ql/6qle-a1-m1-cH.pdb', 'rcsb/ql/6qle-a1-m1-cO.pdb'), \
('rcsb/ql/6qle-a1-m1-cK.pdb', 'rcsb/ql/6qle-a1-m1-cP.pdb'), \
('rcsb/ql/6qle-a1-m1-cK.pdb', 'rcsb/ql/6qle-a1-m1-cQ.pdb'), \
('rcsb/ql/6qle-a1-m1-cK.pdb', 'rcsb/ql/6qle-a1-m1-cY.pdb'), \
('rcsb/ql/6qle-a1-m1-cK.pdb', 'rcsb/ql/6qle-a1-m1-cZ.pdb'), \
('rcsb/ql/6qle-a1-m1-cL.pdb', 'rcsb/ql/6qle-a1-m1-cI.pdb'), \
('rcsb/ql/6qle-a1-m1-cU.pdb', 'rcsb/ql/6qle-a1-m1-cY.pdb'), \
('rcsb/ql/6qle-a1-m1-cU.pdb', 'rcsb/ql/6qle-a1-m2-cN.pdb'), \
('rcsb/ql/6qle-a1-m1-cU.pdb', 'rcsb/ql/6qle-a4-m2-cN.pdb'), \
('rcsb/a6/7a6w-a1-m1-cAAA.pdb', 'rcsb/a6/7a6w-a2-m42-cAAA.pdb'), \
('rcsb/a6/7a6w-a1-m42-cAAA.pdb','rcsb/a6/7a6w-a2-m42-cAAA.pdb')]


#print(count_contact_many(pairs_list=dimers, thresh_max_dist=8, thresh_min_pairs=10, exe_path='/nfs/turbo/umms-petefred/jaschwa/dimerDB-optimize-dev/bin/check_contact/check_contact.exe'))

for dimer in dimers:
    subprocess.run(f'/nfs/turbo/umms-petefred/jaschwa/dimerDB-optimize-dev/bin/USalign/NWalign {dimer[0]} {dimer[1]} -het 1 -glocal 2 -outfmt 2 | grep rcsb', check=True, shell=True)
