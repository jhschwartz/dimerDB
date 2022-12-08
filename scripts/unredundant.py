from abc import ABC, abstractmethod
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import sys
import random
import requests
import glob
import re
import os
import tempfile

MMALIGN='../bin/MMalign'
SPLITPDB = '../bin/PDBParser/split_chain'
NWALIGN = '../../bin/NWalign/align'



class RedundantThings:
    def __init__(self, things: list, threshold: float):
        self.things = things
        self.threshold = threshold
        self.tmp_dir = f'/tmp/{str(random.randint(1,999999))}'

    def __enter__(self):
        return self


    def __exit__(self):
        /


    def prune_redundancy(self, num_workers: int = 1) -> list:
        self.initiate_distance_matrix(num_workers)
        num_clusters = self.initiate_clusters()
        self.non_redundant_things = []
        for cluster_index in range(num_clusters):
            cluster = self.retrieve_cluster(cluster_index)
            rep = self.representative(cluster)
            self.non_redundant_things.append(rep)
        self.rm_all_temp_files()
        return self.non_redundant_things


    def initiate_distance_matrix(self, num_workers: int = 1) -> None:
        # TODO: parallelize with num_workers
        N = len(self.things)
        self.distance_matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(i, N):
                if i == j:
                    continue
                distance = self.distance(self.things[i], self.things[j])
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


    def cluster_index_of_thing(self, thing):
        for i, t in enumerate(self.things):
            if t == thing:
                return self.things_cluster_labels[i]


    def retrieve_cluster(self, cluster_index: int) -> list:
        return [t for i, t in enumerate(self.things) if self.things_cluster_labels[i] == cluster_index]


    @staticmethod
    def _things_of_lowest_distance_to_others(cluster):
        lowest_mean_distance = sys.maxsize
        reps = []
        assert False # fixme
        thing_indecies = [self.things.index(thing) for thing in cluster] # shouldn't look at all things but rather cluster only
        for ti in thing_indecies:                                               # needs a lot of work...what exactly should be happening here?? 
            distances = [self.distance_matrix[other_index] \
                            for other_index in thing_indecies \
                            if ti != other_index]
            avg = np.mean(distances)
            if avg < lowest_mean_distance:
                lowest_mean_distance = avg
                reps = [ self.things[di] ]
            elif avg == lowest_mean_distance:
                reps.append(self.things[di])
        assert False # this has to be a geometric mean
        return reps


    @abstractmethod
    def rm_all_temp_files(self) -> int:
        pass


    @classmethod
    @abstractmethod
    def distance(cls, thing1, thing2) -> float:
        pass


    @staticmethod
    @abstractmethod
    def representative(cluster: list):
        pass



class RedundantDimers(RedundantThings):
    def __init__(self, dimer_tuples, threshold):
        super().__init__(things=dimer_tuples, threshold=threshold)
        self.tmpfiles = self._download_files


    def _filename(self, pdb, chain):
        return f'{self.tmp_dir}/{pdb}_{chain}.pdb'


    @staticmethod
    def _download_pdb(pdb):
        baseurl = 'https://files.wwpdb.org/pub/pdb/data/biounit/PDB/divided/'
        divider = pdb[1:3]
        url = f'{baseurl}/{divider}/{pdb}.pdb1.gz'
        result = requests.get(url)
        if result.status_code != 200:
            raise ConnectionError(f'Failed to retrieve pdb assembly for {pdb}.\
                                    Attempted url was {url} and response was {result.text}')
        pdb_data = result.text
        filename = self._filename(pdb)
        with open(_filename, 'w') as f:
            f.write(pdb_data)
        return filename



    @staticmethod
    def _split_pdb_file(file):
        exe = f'{SPLITPDB} {file}'
        result = subprocess.run(exe, text=True, capture_output=True, shell=True)
        if result.stderr != '':
            raise RuntimeError(f'split_chain for {file} failed with stderr: {result.stderr}')
        file_base_path = file.split('.')[0]
        chain_files = []
        for f in glob.glob(f'{file_base_path}*'):
            if f == file:
                continue
            chain = f.replace(file)
            chain_tuple = (chain, f)
            chain_files.append(chain_tuple)
        return chain_files


    def _download_files(self):
        files = {} # 'pdb': {'A': '/tmp/whatever'}
        for name1, name2 in self.things:
            pdb = name1.split('_')[0]
            if name2.split('_')[0] != pdb:
                raise ValueError(f'chain1 and chain2 of dimer {name1}-{name2} do not \
                                    match in pdb code. This should never happen.')
            file = self._download_pdb(pdb)
            files[pdb]['assembly'] = file
            chains_files = _split_pdb_file(file)
            for chain, chainfile in chains_files:
                files[pdb][chain] = chainfile
        return files


    @staticmethod
    def _dimer_assembly_file(dimer_tuple):
        name1, name2 = dimer_tuple
        pbb = name1.split('_')[0]
        chain1 = name1.split('_')[1]
        chain2 = name2.split('_')[1]
        if chain1+chain2 in self.tmpfiles[pdb]:
            return self.tmpfiles[pdb][chain1+chain2]
        file1 = self.tmpfiles[pdb][chain1]
        file2 = self.tmpfiles[pdb][chain2]
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            dimer_data = f1.read() + '\n' + f2.read()
        outfile = self._filename(pdb, chain1+chain2)
        with open(outfile, 'w') as f:
            f.write(dimer_data)
        self.tmpfiles[pdb][chain1+chain2] = outfile
        return outfile


    def rm_all_temp_files(self) -> int:
        count = 0
        for pdb, files in self.tmpfiles.items():
            for file in files:
                os.remove(file)
                count += 1
        self.tmpfiles = {}
        return count


    @classmethod
    def distance(cls, dimer1_tuple, dimer2_tuple):
        dimer1file = cls._dimer_assembly_file(dimer1_tuple)
        dimer2file = cls._dimer_assembly_file(dimer2_tuple)

        cmd = f'{MMALIGN} {dimer1file} {dimer2file}'
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        if result.stderr != '':
            raise RuntimeError(result.stderr)

        TMscores = re.findall('TM-score=\s([\d\.]+)', result.stdout)
        if len(TMscores) != 2:
            raise ValueError(f'was unable to retrieve both TMscores, instead found {len(TMscores)} scores.')
        score = max([float(tm) for tm in TMscores])

        return 1 - score


    @staticmethod
    def representative(cluster):
        # Criterion 1: keep 25th percentile and up of coverage
        coverages = [self._count_coverage(dimer) for dimer in cluster]
        cutoff = np.percentile(coverages, 25)
        choices = cluster[np.where(coverages >= cutoff)]

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


    def rm_all_temp_files(self) -> int:
        count = 0
        for seq_name, file in self.tmpfiles:
            os.remove(file)
            count += 1
        return count


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









# def example():
#     seq_names = ['ABC123', 'XYZ987', 'ERT456', 'QWE399']
#     threshold = 0.7
#     red_seqs = RedundantSeqs(seq_names, threshold)
#     nonredundant_seq_names = red_seqs.prune_redundancy(num_workers=4)






