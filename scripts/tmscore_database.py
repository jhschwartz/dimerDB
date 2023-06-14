from simple_database import AbstractSimpleDatabase
import name_pdb
import os


class TMDB(AbstractSimpleDatabase):
    '''
    pair is like: (dimer1, dimer2, tmscore_1to2, tmscore_2to1)
    line: "dimer1, dimer2, value1, value2"
    '''

    def group_key(self, pair):
        dimer1, dimer2 = pair[:2]
        if dimer1 < dimer2:
            return name_pdb.get_div(name_pdb.dimer2chains(dimer1)[0])
        return name_pdb.get_div(name_pdb.dimer2chains(dimer2)[0])


    def sort_key(self, pair):
        div_dimer1 = self.group_key(pair)
        return (div_dimer1, *pair)


    def pair_to_line(self, pair):
        dimer1, dimer2, value1, value2 = pair
        return f'{dimer1}\t{dimer2}\t{value1}\t{value2}'


    def line_to_pair(self, line):
        dimer1, dimer2, value1, value2 = line.split()
        return (dimer1, dimer2, value1, value2)


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
        return tuple(line.split())


