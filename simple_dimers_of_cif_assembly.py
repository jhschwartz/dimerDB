import os
import sys
import tempfile
import shutil
import glob
import itertools
from scipy.stats import gmean
import random

sys.path.append('scripts')
from parallel_convert_split_cif import parallel_convert_split_rename_cifs
from align_tools import calc_nwalign_glocal
from wrap_tmscore import calculate_dimers_TM_score

sys.path.append('bin/check_contact')
from check_contact import check_contact_many

def routine(ciffile, td):


    shutil.copy(ciffile, td)

    parallel_convert_split_rename_cifs(td, 'bin/USalign/cif2pdb', 1)

    pdbfiles = []
    for f in glob.glob(f'{td}/*.pdb'):
        pdbfiles.append(f)

    np = []
    for p in pdbfiles:
        n = p+","+str(random.randint(0,9999))
        os.rename(p,n)
        np.append(n)
    pdbfiles = np

    pairs = list(itertools.combinations(pdbfiles, 2))

    contacting = check_contact_many(pairs, 8, 10)

    ids = []
    for pair in pairs:
        id_ = calc_nwalign_glocal('bin/USalign/NWalign', pair[0], pair[1])
        ids.append(id_)

    homodimers = []
    for pair, contact, id_ in zip(pairs, contacting, ids):
        if contact and id_ >= 0.98:
            homodimers.append(pair)
    return homodimers

def compare(homodimers):

    returnable = []
    dimer_pairs = itertools.combinations(homodimers, 2)
    for dimer_pair in dimer_pairs:
        d1 = dimer_pair[0][0].split('/')[-1] + '_' + dimer_pair[0][1].split('/')[-1]
        d2 = dimer_pair[1][0].split('/')[-1] + '_' + dimer_pair[1][1].split('/')[-1]
        if d1.split('_')[0].split('-')[1] != d1.split('_')[1].split('-')[1]:
            continue
        if d2.split('_')[0].split('-')[1] != d2.split('_')[1].split('-')[1]:
            continue
        
        score1, score2 = calculate_dimers_TM_score(*dimer_pair, 'bin/USalign/USalign')
        
        print((d1, d2, score1, score2))

        if d1 <= d2:
            returnable.append( (d1, d2, score1, score2) )
        else:
            returnable.append( (d2, d1, score2, score1) )

    return sorted(returnable)

if __name__ == '__main__':
    code = sys.argv[1]
    path = f'intermediates/sub0/no_sample/assemblies/{code[1:3]}/{code}-ass*'
    assemblies = glob.glob(path)
    with tempfile.TemporaryDirectory() as td:
        homodimers = []
        for ass in assemblies:
            homodimers += routine(ass, td)
    
        results = compare(homodimers)

    print('-------')

    for k, group in itertools.groupby(results, key=lambda r: r[0]):
        scores = []
        for r in group:
            scores.append(float(r[-1]))
        print(k, gmean(scores))



