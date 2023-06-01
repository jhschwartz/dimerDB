import name_pdb
from scipy.stats import gmean
import math


def dimers_of_lowest_distance_to_others(group, dist_mat, all_dimers):
    '''
    From a group of dimers, returns the dimer that is most like all
    the other dimers, meaning that it has the lowest geometric mean
    distance to all the others in the group compared to all the others
    in the group.

    :param group: list[str], each element is a dimer name as defined
                    by name_pdb; this list is the group from which
                    we wish to select a dimer most like the group.
    :param dist_mat: np.ndarray, two-dimensional, gives the distances
                        between all dimers in all_dimers according
                        to the order of all_dimers, where distance
                        is on a scale of 0-1 and 1 is far/dissimilar
    :param all_dimers: list[str], each element is a dimer name as defined
                    by name_pdb; this list is the ordered set of dimers
                    to which group belongs and to which dist_mat corresponds.

    :returns dimers: list[str], each element is a dimer name as defined
                    by name_pdb; this list gives the dimers of group that
                    are most like all the others in group. Usually, this 
                    list has only one string element, but in the case of a
                    tie of geometric mean distances of one dimer to others
                    in group, this list can have multiple elements.
    '''
    gmean_dists_to_others = []

    for this_dimer in group:
        
        this_dimer_dists = []
        
        for other_dimer in group:
            
            if other_dimer == this_dimer:
                continue

            this_index = all_dimers.index(this_dimer)
            other_index = all_dimers.index(other_dimer)

            this2other_dist = dist_mat[this_index, other_index]
            if this2other_dist <= 0:
                this2other_dist = 0.00000001
            this_dimer_dists.append(this2other_dist)

        mean = gmean(this_dimer_dists)
        gmean_dists_to_others.append(mean)

    min_dist = min(gmean_dists_to_others)

    dimers = []
    for mean, dimer in zip(gmean_dists_to_others, group):
        if math.isclose(mean, min_dist):
            dimers.append(dimer)
    
    return dimers



def _num_residues(pdbfile):
    residues_found = set()
    with open(pdbfile, 'r') as f:
        while line := f.readline():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_kind = line[12:16].strip()
                res_num = line[22:26].strip()
                chain = line[21]
                residue_tuple = (res_num, chain)
                if atom_kind in ['CA', 'CB']:
                    residues_found.add(residue_tuple)
    return len(residues_found)


def dimer_coverage(dimer, lib_path):
    mono0, mono1 = name_pdb.dimer2pdbs(dimer, lib_path)
    num0 = _num_residues(mono0)
    num1 = _num_residues(mono1)
    return gmean([num0, num1])
  



def _get_chain_resolu(chain_name, resolu_path): 
    pdb_base, _, _, _ = name_pdb.read_chain_names(chain_name)
    with open(resolu_path, 'r') as f:
        # skip header
        for line in f:
            if line.startswith('-----'):
                break
        # return first match
        for line in f:
            if line.startswith(pdb_base.upper()):
                res = line.split()[2]
                return float(res)
        raise KeyError(f'{pdb_base} not found in resolu file at {resolu_path}')





def get_dimer_avg_resolu(dimer, resolu_path):
    c1, c2 = name_pdb.dimer2chains(dimer)
    r1 = _get_chain_resolu(c1, resolu_path)
    r2 = _get_chain_resolu(c2, resolu_path)
    if r1 == -1.0 or r2 == -1.0:
        return -1
    return gmean([r1, r2])




def dimer_is_xray(dimer, entry_methods_file):
    c1, _ = name_pdb.dimer2chains(dimer)
    pdb, _, _ , _ = name_pdb.read_chain_names(c1)
    
    with open(entry_methods_file, 'r') as f:
        for line in f:
            if line.startswith(pdb):
                _, _, method = line.split()
                if method == 'diffraction':
                    return True
                return False
    
    raise ValueError(f'reached end of {entry_methods_file} and could not find entry {pdb} of {dimer}')


