import numpy as np

# set this to choose how to handle label values which cannot be calculated (for example, if an atom does not exist to calculate a label value.)
INVALID_LABEL = np.nan
ALLOW_MISSING_CB = True



aacode = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



# helper to read only atom lines from the pdb file "fn"
def read_atom_lines(fn):
    with open(fn, 'r') as f:
        lines = f.readlines()
    lines = [line.strip() for line in lines if line.startswith('ATOM')]
    return lines


def residues_equal(r1, r2):
    if r1.keys() != r2.keys():
        return False
    for k, v in r1.items():
        if isinstance(v, np.ndarray):
            if not np.array_equal(v, r2[k]):
                return False
        elif v != r2[k]:
            return False
    return True        

# create a dictionary of residue dictionaries from atom lines
def collect_residues_from_lines(lines):
    residues = []
    curr_num = None
    this_residue = {}

    for atom_line in lines:
        residue_num = atom_line[22:26].strip()

        # reset residue upon encountering new residue number
        if curr_num and residue_num != curr_num:
           residues.append(this_residue)
           this_residue = {}

        res_name = atom_line[17:20].strip()
        atom_name = atom_line[12:16].strip()
        atom_coor = np.zeros([1,3])
        atom_coor[0,0] = float(atom_line[30:38])
        atom_coor[0,1] = float(atom_line[38:46])
        atom_coor[0,2] = float(atom_line[46:54])
        this_residue['res_name'] = res_name
        this_residue[atom_name] = atom_coor

        # save number so we can see if another residue is encountered later
        curr_num = residue_num


    if not residues_equal(residues[-1], this_residue) and this_residue != {}:
        residues.append(this_residue)


    return residues






# filter out residues which lack CAs and raise an exception 
# if any non-GLY residue lacks a CB
def filter_valid_residues(residues):
    residues = [residue for residue in residues if 'CA' in residue]

    if not ALLOW_MISSING_CB and any('CB' not in residue and residue['res_name'] != 'GLY' for residue in residues):
        raise ValueError('Encountered non-Glycine missing a beta carbon. This is not allowed!')

    return residues




# convert the dictionary of residue dictionaries to a singular dictionary.
# this resdic format matches deeppotential and stores keys of the form "XX_YY" where XX is the residue number and YY is the atom type, e.g. "12_CA" would be the residue 12 c-alpha and "55_N" would be the residue 55 nitrogen. Values are coordinates.
def convert_to_resdic_format(residues):
    residues_dic = {}
    nums = []
    seq = []

    for res_num in range(len(residues)):
        nums.append(res_num)
        residue = residues[res_num]
        seq.append(residue['res_name'])
        for atom_name, coordinate in residue.items():
            if atom_name != 'res_name':
                residues_dic[str(res_num)+'_'+atom_name] = coordinate

    return residues_dic, nums, seq    








def read_chain(chain_fn):
    '''
    This function extracts all atoms from a pdb file. It skips atoms which belong to residues which lack a C-alpha atom.
    
    :param chain_fn: the name of the chain's pdb file we are extracting from. Expected that only one chain exists in this file.

    :return: residue_dic, a dictionary which has keys of the form '12_CA' or '56_N' to specify the 12th C-alpha atom and the 56th nitrogen atom, respectively, and values which are (x,y,z) corrdinate lists of size (1, 3). All atom types encountered are included, not just N and CA. Note that residue numbers are 1-indexed.
    :return: nums, a list of arbitrary order of strings of integers representing all resiudes which are encountered for at least one atom type. Note that residue numbers are 1-indexed.
    :return: seq_dic, a dictionary mapping residue numbers to their residue name. Keys are integers for residues numbers and values are 3-letter amino acid codes. Note that residue numbers are 1-indexed.
    '''
   
    # STEP 1: read in ATOM lines only
    lines = read_atom_lines(chain_fn)

    # STEP 2: Read residues.
    # The format is a list of dictionaries, where each member dictionary of the list
    # corresponds to one valid residue. The member dictionaries have keys which correspond
    # to atom names and values which are coordinates. 
    residues = collect_residues_from_lines(lines)

    # STEP 3: Filter out residues which do not have CA atoms. Raise exception if any non-GLY which passes the filter lacks a CB.
    residues = filter_valid_residues(residues)

    # STEP 4: Convert to legacy residue dictionary format.
    residue_dic, nums, seq = convert_to_resdic_format(residues) 

    # STEP 5: Make seq a string of one-letter amino acid codes
    seq = ''.join([aacode.get(aa_three, 'X') for aa_three in seq])

    return residue_dic, nums, seq





def compute_dis(residue_i, residue_j, res_dic, atomtypei, atomtypej):
    '''
    Compute euclidean distance between residues i and j for the given atom types.
    
    :param residue_i: integer of first residue number, 1-indexed.
    :param residue_j: integer of second residue number, 1-indexed.
    :param res_dic: the residue dictionary from the function readpdb.
    :param atomtypei: string of first atom type, like 'CA' or 'CB' or 'N', etc.
    :param atomtypej: string of second atom type, like 'CA' or 'CB' or 'N', etc.
    
    :return: dist, a numpy float representing the calculated distance.
    '''
    
    name_i = str(residue_i)+'_'+atomtypei
    name_j = str(residue_j)+'_'+atomtypej
    
    if not name_i in res_dic or not name_j in res_dic:
        return INVALID_LABEL
    
    resi=res_dic[name_i]
    resj=res_dic[name_j]
    dist = np.sqrt( ((resi-resj)[0]**2).sum() )
    
    return dist





def compute_omg(i, j, res_dic):
    '''
    Compute omega angle between residues i and j.
    
    :param residue_i: integer of first residue number, 1-indexed.
    :param residue_j: integer of second residue number, 1-indexed.
    :param res_dic: the residue dictionary from the function readpdb.
    
    :return: the omega angle as a numpy float.
    '''
    
    try:
        p0 = res_dic[str(i)+'_CA'][0]
        p1 = res_dic[str(i)+'_CB'][0]
        p2 = res_dic[str(j)+'_CB'][0]
        p3 = res_dic[str(j)+'_CA'][0]
    except KeyError:
        return INVALID_LABEL
        
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)+1e-9
    
    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)





def compute_theta(i, j, res_dic):
    '''
    Compute theta angle between residues i and j.
    
    :param residue_i: integer of first residue number, 1-indexed.
    :param residue_j: integer of second residue number, 1-indexed.
    :param res_dic: the residue dictionary from the function readpdb.
    
    :return: the theta angle as a numpy float.
    '''
    
    try:
        p0 = res_dic[str(i)+'_N'][0]
        p1 = res_dic[str(i)+'_CA'][0]
        p2 = res_dic[str(i)+'_CB'][0]
        p3 = res_dic[str(j)+'_CB'][0]
    except KeyError:
        return INVALID_LABEL

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)+1e-9
    
    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    
    return np.arctan2(y, x)   





def unit_vector(vector):
    '''
    Returns the unit vector of vector, where vector is a 1-D numpy array.
    '''
    
    return vector / (np.linalg.norm(vector)+1e-9)





def compute_phi(i, j, res_dic):
    '''
    Compute phi angle between residues i and j.
    
    :param residue_i: integer of first residue number, 1-indexed.
    :param residue_j: integer of second residue number, 1-indexed.
    :param res_dic: the residue dictionary from the function readpdb.
    
    :return: the phi angle as a numpy float.
    '''
    
    try:
        a1 = res_dic[str(i)+'_CA'][0]
        a2 = res_dic[str(i)+'_CB'][0]
        a3 = res_dic[str(j)+'_CB'][0]
    except KeyError:
        return INVALID_LABEL
    
    v1 = a2 - a1
    v2 = a2 - a3
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))     








def extract_labels(res_dic, nums, ca_distance_only=False):
    '''
    Generate ground-truth labels for the structure represented by chain in pdbfile.
    
    :param res_dic: the residue dictionary from the function readpdb.
    :param nums: the list of strings of residue numbers from the function readpdb.
    
    :return: labels, a dictionary containing the label maps for omega, theta, and phi angles, as well as CA and CB distances, and masks for N and CB. 
        NOTE: the maps in labels are zero-indexed, unlike pdb files. So if we want the relationship between residues 12 and 57, for example, we'd look at omega[12-1, 57-1]
        keys of labels:
            'omega': 2D numpy array (LxL) representing omega angles between i and j
            'theta': 2D numpy array (LxL) representing theta angles between i and j
            'phi': 2D numpy array (LxL) representing phi angles between i and j
            'dis': 2D numpy array (LxL) representing CB distances between i and j
            'ca_dis': 2D numpy array (LxL) representing CA distances between i and j
            'mask_cb': 2D numpy array (LxL) representing masked-out CB (i.e. where there is no CB, usually because of glycine). Mask is 1 where CB exists and 0 where it does not.
            'mask_n': 2D numpy array (LxL) representing masked-out N (i.e. where there is no N). Mask is 1 where N exists and 0 where it does not.
        
    '''
    
    nums_ints = [int(n) for n in nums]
    nums_ints.sort()
    max_res = nums_ints[-1] + 1
    
    labels = {}
    
    labels['omega'] = np.full((max_res, max_res), np.nan)
    labels['theta'] = np.full((max_res, max_res), np.nan)
    labels['phi'] = np.full((max_res, max_res), np.nan)
    labels['dis'] = np.full((max_res, max_res), np.nan)
    labels['ca_dis'] = np.full((max_res, max_res), np.nan)

    labels['mask_cb'] = np.ones((max_res, max_res))
    labels['mask_n'] = np.ones((max_res, max_res))
    for i in range(max_res):
        if not str(i) in nums or not str(i)+'_CB' in res_dic:
            labels['mask_cb'][i,:] = 0
            labels['mask_cb'][:,i] = 0
        if not str(i) in nums or not str(i)+'_N' in res_dic:
            labels['mask_n'][i,:] = 0            
            labels['mask_n'][:,i] = 0
   
    for residue_i in nums_ints:
        for residue_j in nums_ints:
            if not ca_distance_only:
                labels['omega'][residue_i,residue_j] = compute_omg(residue_i, residue_j, res_dic)
                labels['theta'][residue_i,residue_j] = compute_theta(residue_i, residue_j, res_dic)
                labels['phi'][residue_i,residue_j] = compute_phi(residue_i, residue_j, res_dic)
                labels['dis'][residue_i,residue_j] = compute_dis(residue_i, residue_j, res_dic, 'CB', 'CB')
            labels['ca_dis'][residue_i,residue_j] = compute_dis(residue_i, residue_j, res_dic, 'CA', 'CA')
            
    return labels







def extract_labels_dimer(res_dic_1, nums_1, res_dic_2, nums_2, ca_distance_only=False):
    '''
    Function to extract dimeric labels from chain1 and chain2 in the same pdbfile.
    Note that these chains HAVE to exist in the same structure as otherwise distances and such would be arbitrary.
    This function works by considering chain1 and chain2 as one concatenated chain, computing labels using extract_labels upon 
    that concatenated chain, then pulling out the two intrachain regions and the interchain region for the returned labels.
    
    :param res_dic_1: the residue dictionary from the function readpdb which describes chain1.
    :param nums_1: the list of strings of residue numbers from the function readpdb which describe chain1.
    :param res_dic_2: the residue dictionary from the function readpdb which describes chain2.
    :param nums_2: the list of strings of residue numbers from the function readpdb which describe chain2.
    
    :return: labels_intra_1, the labels dictionary from function extract_labels which describes the intrachain maps of chain1
    :return: labels_intra_2, the labels dictionary from function extract_labels which describes the intrachain maps of chain2
    :return: labels_inter_12, the labels dictionary from function extract_labels which describes the interchain maps between chains 1 and 2
    :return: L1, int, the length of chain1
    :return: L2, int, the length of chain2
    '''
    
    # figure out L1
    nums_ints_1 = [int(n) for n in nums_1]
    nums_ints_1.sort()
    L1 = nums_ints_1[-1] + 1
    
    # figure out L2
    nums_ints_2 = [int(n) for n in nums_2]
    nums_ints_2.sort()
    L2 = nums_ints_2[-1] + 1
    
    
    # keep chain1 numberings in the merged representation
    res_dic_merged = res_dic_1.copy()
    nums_merged = set(nums_1.copy()) 
    
    # Now, we have to re-number the residues of chain2 to represent the two chains as one concatenated chain.
    for k, v in res_dic_2.items():
        k_split = k.split('_') # e.g. '12_CA' --> ['12', 'CA']
        old_i = int(k_split[0]) # e.g. '12'
        atom_kind = k_split[1] # e.g. 'CA'
        new_i = L1 + old_i # e.g 12 + 148 --> 160
        new_k = str(new_i) + '_' + atom_kind # e.g. '160' + '_' + 'CA' --> '160_CA'
        res_dic_merged[new_k] = v # e.g. res_dic_merged['160_CA'] = [some coords]
        nums_merged.add(str(new_i)) # e.g. nums_merged.add('160')
    
    # check that nums_merged is right
    nums_merged = list(nums_merged)
    if max(map(int, nums_merged)) + 1 != L1 + L2:
        raise ValueError('Read of two chains failed to maintain proper lengths.')
    
    # compute labels by considering the dimer one chain
    double_labels = extract_labels(res_dic_merged, nums_merged, ca_distance_only) 
    
    # pull out the 2 intrachain label dictionaries and the 1 interchain label dictionary using L1 and L2
    labels_intra_1 = {}
    labels_intra_2 = {}
    labels_inter_12 = {}
    for k, v in double_labels.items():
        labels_intra_1[k] = v[:L1, :L1]
        labels_intra_2[k] = v[L1:, L1:]
        labels_inter_12[k] = v[:L1, L1:]
        
    
    return labels_intra_1, labels_intra_2, labels_inter_12, L1, L2


    
