from simple_database import AbstractSimpleDatabase
import name_pdb
import os


class SeqidDB(AbstractSimpleDatabase):
    '''
    pair is like: (chain1, chain2, float)
    line: "chain1, chain2, float"
    '''

    def group_key(self, pair):
        # group by 1 character, the second, like "a" in 1abc
        chain1, chain2, = pair[:2]
        if chain1 < chain2:
            return chain1[1]
        return chain2[1]


    def sort_key(self, pair):
        char_chain1 = self.group_key(pair)
        return (char_chain1, *pair)


    def pair_to_line(self, pair):
        chain1, chain2, seq_id = pair
        return f'{chain1}\t{chain2}\t{seq_id}'


    def line_to_pair(self, line):
        chain1, chain2, seq_id = line.split()
        return (chain1, chain2, float(seq_id))


    def compare_pair_names(self, pair1, pair2):
        name1 = ''.join([*pair1[:2]])
        name2 = ''.join([*pair2[:2]])
        if name1 < name2:
            return -1
        elif name1 > name2:
            return 1
        return 0


    def name_group_file(self, group_name):
        return os.path.join(self.db_path, f'{group_name}.tsv')


    def line_to_output(self, line):
        return self.line_to_pair(line)


