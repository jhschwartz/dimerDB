from simple_database import AbstractSimpleDatabase
import name_pdb
import os


class ContactDB(AbstractSimpleDatabase):
    '''
    pair is like: (chain1, chain2, True/False)
    line: "chain1, chain2, 1/0"
    '''

    def group_key(self, pair):
        # group by 2 characters, the middle two like "ab" in 1abc
        chain1, chain2, = pair[:2]
        if chain1 < chain2:
            return chain1[1:3]
        return chain2[1:3]


    def sort_key(self, pair):
        char_chain1 = self.group_key(pair)
        return (char_chain1, *pair)


    def pair_to_line(self, pair):
        chain1, chain2, in_contact = pair
        int_contact = 0
        if in_contact:
            int_contact = 1
        return f'{chain1}\t{chain2}\t{int_contact}'


    def line_to_pair(self, line):
        chain1, chain2, in_contact = line.split()
        in_contact = bool(in_contact)
        return (chain1, chain2, in_contact)


    def compare_pair_names(self, pair1, pair2):
        chains1 = pair1[:2]
        chains2 = pair2[:2]
        if chains1 == chains2:
            return 0
        ordered = sorted([chains1, chains2])
        index_chains1 = ordered.index(chains1)
        if index_chains1 == 0:
            return -1
        return 1


    def name_group_file(self, group_name):
        return os.path.join(self.db_path, f'{group_name}.tsv')


    def line_to_output(self, line):
        out = line.split()
        out[-1] = out[-1] == '1'
        return tuple(out)



