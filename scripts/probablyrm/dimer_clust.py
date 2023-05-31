import numpy as np
from sklearn.cluster import AgglomerativeClustering
import sys
from scipy.stats import gmean
import subprocess
from multiprocessing import Pool



def cluster_from_mat(mat):
    if len(self.mat) == 1:
        return [0]
    cluster = AgglomerativeClustering(  distance_threshold=self.threshold,
                                        affinity='precomputed',
                                        linkage='complete',
                                        n_clusters=None
                                                                            )
    return cluster.fit(mat).labels_




def get_cluster(dimers, cluster_labels, cluster_index):
    for dimer, index in zip(dimers, cluster_labels):
        if cluster_label == cluster_index:
            yield dimer





def dimer_most_like_cluster(cluster_index, dimers, cluster_labels, mat):
    mean_dist_to_others = []

    for index, dimer in enumerate(dimers):
        if cluster_labels[index] != cluster_index:
            continue

        # dists to other dimers on left of dimer
        for other_index in range(index):
            dist = mat[index, other_index]
            dist = max(dist, 0.00000001) # avoid gmean divide by zero
            dists.append(dist)

        # dists to other dimers on right of dimer
        for other_index in range(index+1, len(dimers)):
            dist = mat[index, other_index]
            dist = max(dist, 0.00000001) # avoid gmean divide by zero
            dists.append(dist)

        mean = gmean(dists)
        mean_dist_to_others.append((mean, dimer))
    
    min_ = None
    dimers_to_return = []
    for mean, dimer in mean_dist_to_others:
        if min_ is None or mean < min_:
            min_ = mean
            dimers_to_return = [ dimer ]
        elif mean == min_:
            dimers_to_return.append(dimer)
    return dimers_to_return








    #@staticmethod
    #@contextlib.contextmanager
    #def _tmp_assembly_file(pdb_path_0: str, pdb_path_1: str):
    #    with open(pdb_path_0, 'r') as f0, open(pdb_path_1, 'r') as f1:
    #        dimer_data = f0.read() + '\n' + f1.read()
    #    
    #    with tempfile.NamedTemporaryFile(mode='w') as f:
    #        f.write(dimer_data)
    #        f.seek(0)
    #        yield f.name


    #def distance(self, dimer0_name: str, dimer1_name: str) -> float:
    #    dimer0_pdb0, dimer0_pdb1 = dimer2pdbs(dimer0_name, self.config['paths']['lib'])
    #    dimer1_pdb0, dimer1_pdb1 = dimer2pdbs(dimer1_name, self.config['paths']['lib'])
    #    with RedundantDimerStructures._tmp_assembly_file(dimer0_pdb0, dimer0_pdb1) as d0file, \
    #            RedundantDimerStructures._tmp_assembly_file(dimer1_pdb0, dimer1_pdb1) as d1file:
    #        cmd = f'{self.config["paths"]["usalign"]} -mm 1 -ter 1 {d0file} {d1file}'
    #        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    #        if result.stderr != '':
    #            raise RuntimeError(result.stderr)

    #        TMscores = re.findall(r'TM-score=\s([\d\.]+)', result.stdout)
    #        if len(TMscores) != 2:
    #            raise ValueError(f'was unable to retrieve both TMscores, instead found {len(TMscores)} scores.')
    #        score = max([float(tm) for tm in TMscores])

    #        return 1 - score

    
    @staticmethod
    def _num_residues(pdbfile: str) -> int:
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


    def _dimer_coverage(self, dimer_name: str) -> float:
        mono0, mono1 = dimer2pdbs(dimer_name, self.config['paths']['lib'])
        num0 = RedundantDimerStructures._num_residues(mono0)
        num1 = RedundantDimerStructures._num_residues(mono1)
        return gmean([num0, num1])


    def _get_chain_resolu(self, chain_name):
        pdb_base, _, _, _ = read_chain_names(chain_name)
        with open(self.config['paths']['resolu_file'], 'r') as f:
            # skip header
            for line in f:
                if line.startswith('-----'):
                    break
            # return first match
            for line in f:
                if line.startswith(pdb_base.upper()):
                    res = line.split()[2]
                    return float(res)
            raise KeyError(f'{pdb_base} not found in resolu file at {self.config["paths"]["resolu_file"]}')


    def _get_dimer_resolu(self, dimer_name):
        d1 = name_pdb.dimer2chains(dimer_name)[0]
        d1, d2 = dimer_name.split('-')
        r1 = self._get_chain_resolu(d1)
        r2 = self._get_chain_resolu(d2)
        if r1 == -1.0 or r2 == -1.0:
            return -1
        return gmean([r1, r2])




# TODO: move to smk file! 
    def representative(self, cluster: list, prefer_xray=True) -> str:
        # Criterion 1: keep only structures with greater than half the maximum coverage,
        #               where coverage is defined as the gemoetric mean number of residues
        #               in two chains
        coverages = [self._dimer_coverage(dimer_name) for dimer_name in cluster]
        cutoff = max(coverages)/2 
        choices = []
        for index, dimer in enumerate(cluster):
            if coverages[index] >= cutoff:
                choices.append(dimer)
        
        # end if no tie
        if len(choices) == 1:
            return choices[0]

        resolus = {}

        # Criterion 1.5: preference xray structures if prefer_xray is True
        #       the actual rule: if there exists an xray choice, then preference xray choices
        if prefer_xray:
            
            xray_choice_exists = False

            # extract resolus 
            for dimer in choices:
                resolus[dimer] = self._get_dimer_avg_resolu(dimer)
                if resolus[dimer] != -1.0:
                    xray_choice_exists = True
            
            # if any xray exist, they become choices
            if xray_choice_exists:
                choices = [c for c in choices if resolus[c] != -1.0]
                
                # end if no tie
                if len(choices) == 1:
                    return choices[0]


        # Criterion 2: pick dimers with lowest average distance to remaining
        choices = super()._dimers_of_lowest_distance_to_others(choices)

        # end if no tie
        if len(choices) == 1:
            return choices[0]

        # Criterion 2.5: preference best xray structure if prefer_xray
        #           actual rule: if an xray choice exists, choose the best resolution xray choice
        if prefer_xray:
            # Note, only xray results would make it this far.
            leftover_resolus = [resolus[c] for c in choices]
            lowest = sorted(leftover_resolus)[0]
            best_choices = []
            for choice, resolu in zip(choices, leftover_resolus):
                if resolu == lowest:
                    best_choices.append(choice)
            choices = best_choices


        # Criterion 3: pick the newest, assuming higher in alphabetical order is newer
        return sorted(choices, key=lambda dimer: dimer.split('-')[1], reverse=True)[0]





