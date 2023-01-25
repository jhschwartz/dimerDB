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
    - test/test_unredundant_dimer_seqs.py - tests RedundantSeqsHomodimer and RedundantSeqsHeterodimer, untested as of 1/12/23.
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
import yaml
import threading
from multiprocessing import Pool
from read_fasta import read_prot_from_fasta
from name_pdb import dimer2pdbs
from name_fasta import get_homodimer_fasta, get_heterodimer_fasta


class RedundantThings:
    def __init__(self, things: list, threshold: float):
        self.things = things
        self.threshold = threshold


    def prune_redundancy(self, num_workers=1) -> list:
        self.initiate_distance_matrix(num_workers)
        num_clusters = self.initiate_clusters()
        self.non_redundant_things = []
        for cluster_index in range(num_clusters):
            cluster = self.retrieve_cluster(cluster_index)
            rep = self.representative(cluster)
            self.non_redundant_things.append(rep)
        return self.non_redundant_things



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
            cmd = f'{self.config["paths"]["mmalign_exe"]} {d0file} {d1file}'
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


    def representative(self, cluster: list) -> str:
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

        # Criterion 2: pick dimers with lowest average distance to remaining
        choices = super()._things_of_lowest_distance_to_others(choices)

        # end if no tie
        if len(choices) == 1:
            return choices[0]

        # Criterion 3: pick the newest, assuming higher in alphabetical order is newer
        return sorted(choices, key=lambda dimer_tuple: dimer_tuple[0])[-1]





class RedundantSeqs(RedundantThings): 
    def __init__(self, dimer_names, yamlfile, threshold, config):
        super().__init__(things=dimer_names, threshold=threshold)
        self.config = config
        with open(yamlfile, 'r') as f:
            self.dimers = yaml.safe_load(f)
        self.nw = self.config['paths']['nwalign']


    @staticmethod
    def _calc_nw(nw, fasta1: str, fasta2: str) -> float:
        exe = f'{nw} {fasta1} {fasta2}'
        result = subprocess.run(exe, text=True, capture_output=True, shell=True)
        if result.stderr != '':
            raise RuntimeError(f'nwalign for seqs {fasta1} and {fasta2} failed with stderr: {result.stderr}')
        # The line we want looks like "Sequence identity:    0.625 (=   5/   8)"
        identity = float(re.findall(r'Sequence identity:\s+([\d\.]+)', result.stdout)[0])
        return identity
        
    
    @staticmethod
    def max_both_ways_nw(nw, fasta1: str, fasta2: str) -> float:
        left = RedundantSeqs._calc_nw(nw, fasta1, fasta2)
        right = RedundantSeqs._calc_nw(nw, fasta2, fasta1)
        if left > right:
            return left
        return right


    @abstractmethod
    def _count_dimer_length(self, dimer_name: str) -> int:
        pass


    def representative(self, seq_cluster: list) -> str:
        # Criterion 1: pick seq with most nonredundant structures
        counter = {}
        for dimer_name in seq_cluster:
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
        for dimer_name in seq_cluster:
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
    def __init__(self, dimer_names, yamlfile, threshold, config):
        super().__init__(dimer_names, yamlfile, threshold, config)
        self.fd = self.config['paths']['intermediates_homodimer_filtering']

    
    def _count_dimer_length(self, dimer_name: str) -> int:
        fasta = get_homodimer_fasta(uniparc_id=dimer_name, filtering_dir=self.fd)
        _, seq = next(read_prot_from_fasta(fasta))
        return 2*len(seq)
    
    
    def distance(self, dimer1name: str, dimer2name: str) -> float:
        fasta1 = get_homodimer_fasta(uniparc_id=dimer1name, filtering_dir=self.fd)
        fasta2 = get_homodimer_fasta(uniparc_id=dimer2name, filtering_dir=self.fd)
        nw_val = super().max_both_ways_nw(self.nw, fasta1, fasta2)
        return 1 - nw_val


class RedundantSeqsHeterodimer(RedundantSeqs):
    def __init__(self, dimer_names, yamlfile, threshold, config):
        super().__init__(dimer_names, yamlfile, threshold, config)
        self.fd = self.config['paths']['intermediates_heterodimer_filtering']
  

    def _count_dimer_length(self, dimer_name: str) -> int:
        uniparc1, uniparc2 = dimer_name.split('-')
        fasta1, fasta2 = get_heterodimer_fasta(dimer_name=dimer_name, filtering_dir=self.fd)
        _, seq1 = next(read_prot_from_fasta(fasta1))
        _, seq2 = next(read_prot_from_fasta(fasta2))
        return len(seq1) + len(seq2)


    def distance(self, dimer1name: str, dimer2name: str) -> float:
        fasta1A, fasta1B = get_heterodimer_fasta(dimer_name=dimer1name, filtering_dir=self.fd)
        fasta2A, fasta2B = get_heterodimer_fasta(dimer_name=dimer2name, filtering_dir=self.fd)

        # must first figure out if the correspondance is 1A->2A or 1A->2B (and 1B->2B or 1B->2A)
        #   we call the 1A->2A / 1B->2B case "forwards"
        #   we call the 1A->2B / 1B->2A case "backwards"
        nw_1A_2A = super().max_both_ways_nw(self.nw, fasta1A, fasta2A)
        nw_1A_2B = super().max_both_ways_nw(self.nw, fasta1A, fasta2B)
        nw_1B_2A = super().max_both_ways_nw(self.nw, fasta1B, fasta2A)
        nw_1B_2B = super().max_both_ways_nw(self.nw, fasta1B, fasta2B)

        forward_tot = nw_1A_2A + nw_1B_2B
        backward_tot = nw_1A_2B + nw_1B_2A

        if forward_tot > backward_tot:
            nw_I = nw_1A_2A
            nw_II = nw_1B_2B
        else:
            nw_I = nw_1A_2B
            nw_II = nw_1B_2A

        avg_nw = (nw_I + nw_II) / 2
        return 1 - avg_nw


