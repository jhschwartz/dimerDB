'''
unredundant.py - classes to remove redundancy in data using clustering. Specifically, below an abstract class
                 RedundantThings is defined, which is meant to remove redundancy in a list of some datatype.
                 RedundantThings uses hierarchial agglomerative clustering to define clusters. Concrete
                 subclasses require only two abstract methods be implemented:
                    - distance: this method defines the difference between two "things"
                    - representative: this method defines how we select a representative "thing" from a 
                                      cluster. Note that only this representative will remain in the
                                      resulting data and all others in a cluster are ignored as redundant.

Written by Jacob Schwartz (jaschwa@umich.edu) November 2022 - January 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu

This function is unittested by:
    - test/test_unredundant_generic.py - tests abstract class RedundantThings, passing as of 1/12/2023.
    - test/test_unredundant_dimer_structures.py - tests RedundantDimerStructures, passing as of 1/12/2023.
    - test/test_unredundant_dimer_seqs.py - tests RedundantSeqsHomodimer, passing as of 1/25/2023.
This work requires python >= 3.8
'''

from abc import ABC, abstractmethod
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import sys
import requests
import re
import os
import shutil
import tempfile
from scipy.stats import gmean
import contextlib
import subprocess
import pickle
import threading
from multiprocessing import Pool
from read_fasta import read_prot_from_fasta
from name_pdb import dimer2pdbs, read_chain_names
from name_fasta import uniparc_fasta
from align_tools import calc_nwalign


class RedundantThings:
    def __init__(self, things: list, threshold: float):
        self.things = things
        self.threshold = threshold
        self.distance_matrix = None


    def prune_redundancy(self, num_workers=1, calc_dist_matrix=True, rep_extra_kwargs={}) -> list:
        if calc_dist_matrix:
            self.initiate_distance_matrix(num_workers)
        elif self.distance_matrix is None:
            raise ValueError('need a distance matrix to prune redundancy')
        num_clusters = self.initiate_clusters()
        self.non_redundant_things = []
        for cluster_index in range(num_clusters):
            cluster = self.retrieve_cluster(cluster_index)
            rep = self.representative(cluster=cluster, **rep_extra_kwargs)
            self.non_redundant_things.append(rep)
        return self.non_redundant_things


    def save_dist_matrix(self, savefile: str) -> None:
        if self.distance_matrix is None:
            raise ValueError('cannot save nonexistent distance matrix!')
        np.save(savefile, self.distance_matrix)


    def load_dist_matrix(self, loadfile: str) -> None:
        self.distance_matrix = np.load(loadfile)


    def _distance_thread_helper(self, i: int, j: int) -> float:
        return self.distance(self.things[i], self.things[j])   

    def initiate_distance_matrix(self, num_workers=1) -> None:
        N = len(self.things)
        self.distance_matrix = np.zeros((N,N))
       
        args = []
        for i in range(N):
            for j in range(i+1, N):
                args.append((i,j))

        with Pool(processes=num_workers) as p:
            distances = p.starmap(self._distance_thread_helper, args)
        
        for distance, (i, j) in zip(distances, args):
            self.distance_matrix[i,j] = distance
            self.distance_matrix[j,i] = distance




    def initiate_clusters(self) -> None:
        if len(self.distance_matrix) == 1:
            self.things_cluster_labels = [0]
            return 1 
        cluster = AgglomerativeClustering(distance_threshold=self.threshold,
                                          affinity='precomputed',
                                          linkage='complete',
                                          n_clusters=None
                                        )
        self.things_cluster_labels = cluster.fit(self.distance_matrix).labels_
        return len(set(self.things_cluster_labels))


    def cluster_index_of_thing(self, thing) -> int:
        for i, t in enumerate(self.things):
            if t == thing:
                return self.things_cluster_labels[i]


    def retrieve_cluster(self, cluster_index: int) -> list:
        return [t for i, t in enumerate(self.things) if self.things_cluster_labels[i] == cluster_index]


    def _things_of_lowest_distance_to_others(self, cluster: list) -> list:
        mean_dist_to_others = []
        for index, thing in enumerate(cluster):
            dists = [self.distance(thing, other_thing) for other_thing in cluster]
            dists.pop(index) # take out self distance

            for i, d in enumerate(dists):
                # avoid gmean divide by zero
                if d == 0:
                    dists[i] = 0.00000001
            
            geometric_mean = gmean(dists)
            mean_dist_to_others.append(geometric_mean)
        
        things_to_return = []
        for i, thing in enumerate(cluster):
            if mean_dist_to_others[i] == min(mean_dist_to_others):
                things_to_return.append(thing)
        return things_to_return


    @abstractmethod
    def distance(self, thing1, thing2) -> float:
        pass


    @abstractmethod
    def representative(self, cluster: list) -> str:
        pass



class RedundantDimerStructures(RedundantThings):
    def __init__(self, dimer_names: list, threshold: float, config: str):
        super().__init__(things=dimer_names, threshold=threshold)
        self.config = config


    @staticmethod
    @contextlib.contextmanager
    def _tmp_assembly_file(pdb_path_0: str, pdb_path_1: str):
        with open(pdb_path_0, 'r') as f0, open(pdb_path_1, 'r') as f1:
            dimer_data = f0.read() + '\n' + f1.read()
        
        with tempfile.NamedTemporaryFile(mode='w') as f:
            f.write(dimer_data)
            f.seek(0)
            yield f.name


    def distance(self, dimer0_name: str, dimer1_name: str) -> float:
        dimer0_pdb0, dimer0_pdb1 = dimer2pdbs(dimer0_name, self.config['paths']['lib'])
        dimer1_pdb0, dimer1_pdb1 = dimer2pdbs(dimer1_name, self.config['paths']['lib'])
        with RedundantDimerStructures._tmp_assembly_file(dimer0_pdb0, dimer0_pdb1) as d0file, \
                RedundantDimerStructures._tmp_assembly_file(dimer1_pdb0, dimer1_pdb1) as d1file:
            cmd = f'{self.config["paths"]["usalign"]} -mm 1 -ter 1 {d0file} {d1file}'
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if result.stderr != '':
                raise RuntimeError(result.stderr)

            TMscores = re.findall(r'TM-score=\s([\d\.]+)', result.stdout)
            if len(TMscores) != 2:
                raise ValueError(f'was unable to retrieve both TMscores, instead found {len(TMscores)} scores.')
            score = max([float(tm) for tm in TMscores])

            return 1 - score

    
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


    def _get_dimer_avg_resolu(self, dimer_name):
        d1, d2 = dimer_name.split('-')
        r1 = self._get_chain_resolu(d1)
        r2 = self._get_chain_resolu(d2)
        if r1 == -1.0 or r2 == -1.0:
            return -1
        return gmean([r1, r2])


    def representative(self, cluster: list, prefer_xray=False) -> str:
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
        choices = super()._things_of_lowest_distance_to_others(choices)

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





class RedundantSeqs(RedundantThings): 
    def __init__(self, dimer_names, datadict, threshold, config):
        super().__init__(things=dimer_names, threshold=threshold)
        self.dimers = datadict
        self.config = config
        self.nw = self.config['paths']['nwalign']


        
    
    @staticmethod
    def max_both_ways_nw(nw, fasta1: str, fasta2: str) -> float:
        left = calc_nwalign(nw, fasta1, fasta2)
        right = calc_nwalign(nw, fasta2, fasta1)
        if left > right:
            return left
        return right


    @abstractmethod
    def _count_dimer_length(self, dimer_name: str) -> int:
        pass


    def representative(self, cluster: list) -> str:
        # Criterion 1: pick seq with most nonredundant structures
        counter = {}
        for dimer_name in cluster:
            counter[dimer_name] = len(self.dimers[dimer_name])
        max_count = max(counter.values())
        best_dimer_names = [dimer_name for dimer_name, count in counter.items() if count == max_count]


        # end if no tie
        if len(best_dimer_names) == 1:
            return best_dimer_names[0]

        # Criterion 2: pick seq most like others
        best_dimer_names = super()._things_of_lowest_distance_to_others(best_dimer_names)

        # end if no tie
        if len(best_dimer_names) == 1:
            return best_dimer_names[0]

        # Criterion 3: pick the longest seq
        best_dimer_names = []
        max_len = -1
        for dimer_name in cluster:
            L = self._count_dimer_length(dimer_name)
            if L > max_len:
                best_dimer_names = [ dimer_name ]
                max_len = L
            elif L == max_len:
                best_dimer_names.append(dimer_name)

        # end if no tie
        if len(best_dimer_names) == 1:
            return best_dimer_names[0]

        # Criterion 4: pick the first seq name alphabetically
        return sorted(best_dimer_names)[0]


class RedundantSeqsHomodimer(RedundantSeqs):
    def __init__(self, dimer_names, datadict, threshold, config):
        super().__init__(dimer_names, datadict, threshold, config)
        self.lib = self.config['paths']['lib']

    
    def _count_dimer_length(self, dimer_name: str) -> int:
        fasta = uniparc_fasta(uniparc_id=dimer_name, lib_path=self.lib)
        _, seq = next(read_prot_from_fasta(fasta))
        return 2*len(seq)
    
    
    def distance(self, dimer1name: str, dimer2name: str) -> float:
        fasta1 = uniparc_fasta(uniparc_id=dimer1name, lib_path=self.lib)
        fasta2 = uniparc_fasta(uniparc_id=dimer2name, lib_path=self.lib)
        nw_val = super().max_both_ways_nw(self.nw, fasta1, fasta2)
        return 1 - nw_val
