from abc import ABC, abstractmethod
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import sys
import requests
import re
import os
import tempfile
from scipy.stats import gmean
import contextlib
import subprocess
import yaml
from name_pdb import dimer2pdbs
import threading
from multiprocessing import Pool


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




    def _distance_thread_helper(self, i, j):
        distance = self.distance(self.things[i], self.things[j])
        self.distance_matrix[i,j] = distance
        self.distance_matrix[j,i] = distance
        return self.distance_matrix 
    
    def _distance_thread_helper_2(self, i, j):
        return self.distance(self.things[i], self.things[j])

    
    def initiate_distance_matrix(self, num_workers=1) -> None:
        N = len(self.things)
        self.distance_matrix = np.zeros((N,N))
       
        ## single thread method
       # for i in range(N):
       #     for j in range(i, N):
       #         if i == j:
       #             continue
       #         distance = self.distance(self.things[i], self.things[j])
       #         self.distance_matrix[i,j] = distance
       #         self.distance_matrix[j,i] = distance



        ## multithreaded method 3
        args = []
        for i in range(N):
            for j in range(i, N):
                if i == j:
                    continue
                args.append((i,j))

        with Pool(processes=num_workers) as p:
            distances = p.starmap(self._distance_thread_helper_2, args)

        for distance, (i, j) in zip(distances, args):
            self.distance_matrix[i,j] = distance
            self.distance_matrix[j,i] = distance

        ## multithreaded method 2
        #args = []
        #for i in range(N):
        #    for j in range(i, N):
        #        if i == j:
        #            continue
        #        args.append((i,j))

        #with Pool(processes=num_workers) as p:
        #    print(p.starmap(self._distance_thread_helper, args))

        
        # multithreaded method 1
        #pool = []
        #for i in range(N):
        #    for j in range(i, N):
        #        if i == j:
        #            continue
        #        thread = threading.Thread(target=self._distance_thread_helper, args=(i, j))
        #        pool.append(thread)
        #
        #for thread in pool:
        #    thread.start()
        #for thread in pool:
        #    thread.join()


    def initiate_clusters(self) -> None:
        cluster = AgglomerativeClustering(distance_threshold=self.threshold,
                                          affinity='precomputed',
                                          linkage='complete',
                                          n_clusters=None
                                        )
        self.things_cluster_labels = cluster.fit(self.distance_matrix).labels_
        return len(set(self.things_cluster_labels))


    def cluster_index_of_thing(self, thing):
        for i, t in enumerate(self.things):
            if t == thing:
                return self.things_cluster_labels[i]


    def retrieve_cluster(self, cluster_index: int) -> list:
        return [t for i, t in enumerate(self.things) if self.things_cluster_labels[i] == cluster_index]


    def _things_of_lowest_distance_to_others(self, cluster):
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


    @staticmethod
    @abstractmethod
    def representative(cluster: list):
        pass



class RedundantDimers(RedundantThings):
    def __init__(self, dimer_names, threshold, config):
        super().__init__(things=dimer_names, threshold=threshold)
        self.config = config


    @staticmethod
    @contextlib.contextmanager
    def _tmp_assembly_file(pdb_path_0, pdb_path_1):
        with open(pdb_path_0, 'r') as f0, open(pdb_path_1, 'r') as f1:
            dimer_data = f0.read() + '\n' + f1.read()
        
        with tempfile.NamedTemporaryFile(mode='w') as f:
            f.write(dimer_data)
            f.seek(0)
            yield f.name


    def distance(self, dimer0_name, dimer1_name):
        dimer0_pdb0, dimer0_pdb1 = dimer2pdbs(dimer0_name, self.config['paths']['lib'])
        dimer1_pdb0, dimer1_pdb1 = dimer2pdbs(dimer1_name, self.config['paths']['lib'])
        with RedundantDimers._tmp_assembly_file(dimer0_pdb0, dimer0_pdb1) as d0file, \
                RedundantDimers._tmp_assembly_file(dimer1_pdb0, dimer1_pdb1) as d1file:
            cmd = f'{self.config["paths"]["mmalign_exe"]} {d0file} {d1file}'
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if result.stderr != '':
                raise RuntimeError(result.stderr)

            TMscores = re.findall('TM-score=\s([\d\.]+)', result.stdout)
            if len(TMscores) != 2:
                raise ValueError(f'was unable to retrieve both TMscores, instead found {len(TMscores)} scores.')
            score = max([float(tm) for tm in TMscores])

            return 1 - score

    
    @staticmethod
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


    def _dimer_coverage(self, dimer_name):
        mono0, mono1 = dimer2pdbs(dimer_name, self.config['paths']['lib'])
        num0 = RedundantDimers._num_residues(mono0)
        num1 = RedundantDimers._num_residues(mono1)
        return gmean([num0, num1])


    def representative(self, cluster):
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
        return sorted(choices, key=lambda dimer_tuple: dimer_tuple[0].split('_')[0])[-1]





class RedundantSeqs(RedundantThings):
    def __init__(self, seq_names, threshold, homodimers_dict):
        super().__init__(things=seq_names, threshold=threshold)
        # self.homodimers_yaml_file = homodimers_yaml_file
        self.homodimers_dict = homodimers_dict
        self.tmpfiles = {}


    def _retrieve_seqfile(self, seq_name):
        if seq_name in self.tmpfiles:
            return self.tmpfiles[seq_name]
        url = 'https://rest.uniprot.org/uniprotkb/{seq_name}.fasta'
        result = requests.get(url)
        if result.status_code != 200:
            raise ConnectionError(f'Failed to retrieve fasta for {seqname}.\
                                    Attempted url was {url} and response was {result.text}')
        fasta_data = result.text
        fasta_file = f'{self.tmp_dir}/{seq_name}.fasta'
        with open(fasta_file, 'w') as f:
            f.write(fasta_data)
        self.tmpfiles[seq_name] = fasta_file
        return fasta_file



    @classmethod
    def distance(cls, seq1name, seq2name):
        seq1file = cls._retrieve_seqfile(seq1name)
        seq2file = cls._retrieve_seqfile(seq2name)

        exe = f'{NWALIGN} {seq1file} {seq2file}'
        result = subprocess.run(exe, text=True, capture_output=True, shell=True)
        if result.stderr != '':
            raise RuntimeError(f'nwalign for seqs {seq1name} and {seq2name} failed with stderr: {result.stderr}')
        # The line we want looks like "Sequence identity:    0.625 (=   5/   8)"
        identity = float(re.findall('Sequence identity:=\s([\d\.]+)', result.stdout)[0])
        return identity


    # _homodimer_yaml_count_next_uniprot(self):
    #     with open(self.homodimers_yaml_file, 'r') as f:
    #         uniprot = f.readline().strip().split(':')[0]
    #         count = 0
    #         while line := f.readline().strip():
    #             if ':' in line:
    #                 yield uniprot, count
    #                 count = 0
    #                 uniprot = line.split(':')[0]
    #             else:
    #                 count += 1



    @staticmethod
    def representative(seq_cluster):
        # Criterion 1: pick seq with most nonredundant structures

        # count_dir = {}
        # for uniprot, count in self._homodimer_yaml_count_next_uniprot():
        #     for seq_name in seq_cluster:
        #         if seq_name in uniprot:
        #             count_dict[seq_name] = count
        # if len(count_dict.items()) != len(seq_cluster):
        #     raise ValueError('unable to count structures for each sequence.')

        # max_count = np.max(count_dict.values())
        # best_seqs = [name for name, count in count_dict.items if count == max_count]


        best_seqs = []
        max_structs = -1
        for seq_name in seq_cluster:
            num_structs = len(self.homodimers_dict[seq_name])
            if num_structs > max_structs:
                best_seqs = [ seq_name ]
                max_structs = num_structs
            elif num_structs == max_structs:
                best_seqs.append(seq_name)

        # end if no tie
        if len(best_seqs) == 1:
            return best_seqs[0]

        # Criterion 2: pick seq most like others
        best_seqs = super()._things_of_lowest_distance_to_others(best_seqs)

        # end if no tie
        if len(best_seqs) == 1:
            return best_seqs[0]

        # Criterion 3: pick the longest seq
        lengths = []
        for seq_name in best_seqs:
            seq_fasta = self._retrieve_seqfile(seq_name)
            with open(seq_fasta, 'r') as f:
                lines = [line.strip() for line in f.readlines() if not '>' in line]
                seq = ''.join(lines)
                lengths.append(len(seq))
        max_len = np.max(lengths)
        best_seqs = best_seqs[np.where(lengths == max_len)]

        # end if no tie
        if len(best_seqs) == 1:
            return best_seqs[0]

        # Criterion 4: pick the first seq name alphabetically
        return sorted(seq_names)[0]

