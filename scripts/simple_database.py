import os
import itertools
import tempfile
from name_pdb import get_div
import abc


class AbstractSimpleDatabase(abc.ABC):

    def __init__(self, db_path, max_buffer=1000):
        self.db_path = db_path
        self.buffer = []
        self.max_buffer = max_buffer

    

    @abc.abstractmethod
    def group_key(self, pair):
        pass


    @abc.abstractmethod
    def sort_key(self, pair):
        pass

    
    @abc.abstractmethod
    def pair_to_line(self, pair):
        pass


    @abc.abstractmethod
    def line_to_pair(self, line):
        pass


    @abc.abstractmethod
    def compare_pair_names(self, pair1, pair2):
        pass


    @abc.abstractmethod
    def name_group_file(self, group_name):
        pass

    
    # add a pair to the buffer for later update
    def buffer_update(self, name_value_pair):
        self.buffer.append(name_value_pair)



    # add a pair to the buffer, and if the buffer reaches size max_buffer, perform update
    def auto_buffer(self, name_value_pair):
        self.buffer_update(name_value_pair)
        if len(self.buffer) >= self.max_buffer:
            self.update_clear_buffer()



    # use all values in the buffer for update and clear it
    def update_clear_buffer(self):
        self.update(self.buffer)
        self.buffer = []



    def update(self, new_name_value_pairs): #, already_sorted=False):
        #if not already_sorted:
        sorted_pairs = sorted(new_name_value_pairs, key=self.sort_key)
        #else:
        #    sorted_pairs = new_name_value_pairs
        os.makedirs(self.db_path, exist_ok=True)

        for gkey, group in itertools.groupby(sorted_pairs, key=self.group_key):

            db_file = self.name_group_file(group_name=gkey)

            if not os.path.exists(db_file):
                with open(db_file, 'w') as f:
                    for pair in group:
                        line = self.pair_to_line(pair)
                        f.write(f'{line}\n')
        
            
            else:
                with tempfile.NamedTemporaryFile('w+t') as ftmp:
                    
                    with open(db_file, 'r') as fxx:
                        
                        line = fxx.readline()
                        new_pair = next(group, None)
                        
                        while new_pair and line:
                            
                            # new_pair is alphabetically first, so write it and advance new
                            if self.compare_pair_names(new_pair, self.line_to_pair(line)) < 0:
                                new_line = self.pair_to_line(new_pair)
                                ftmp.write(f'{new_line}\n')
                                new_pair = next(group, None)

                            # new and old values are same, so pick new and skip old
                            elif self.compare_pair_names(new_pair, self.line_to_pair(line)) == 0:
                                new_line = self.pair_to_line(new_pair)
                                ftmp.write(f'{new_line}\n')
                                new_pair = next(group, None)
                                line = fxx.readline()

                            # old value is alphabetically first, so write it and advance old
                            else:
                                ftmp.write(line)
                                line = fxx.readline()

                        # write any leftovers from old pair
                        while line:
                            ftmp.write(line)
                            line = fxx.readline()

                        # write any leftovers from new_pair
                        while new_pair:
                            new_line = self.pair_to_line(new_pair)
                            ftmp.write(f'{new_line}\n')
                            new_pair = next(group, None)

                    # now that the updated "xx.tsv" is in a tempfile, overwrite the original
                    with open(db_file, 'w') as fxx:
                        ftmp.seek(0)
                        for line in ftmp:
                            fxx.write(line)



    def get(self, names, strict=False):

        #if not already_sorted:
        unvalued_pairs = ((*n, '') for n in sorted(names, key=self.sort_key))
        #else:
        #    unvalued_pairs = ((*n, '') for n in names)
        
        for gkey, group in itertools.groupby(unvalued_pairs, key=self.group_key):
            
            db_file = self.name_group_file(group_name=gkey)

            if not os.path.exists(db_file):
                if strict:
                    raise ValueError(f'unable to open database file for group {gkey} at path {self.db_path}')
                else:
                    continue
       
            with open(db_file, 'r') as infile:

                current_pair = next(group, None)
                line = infile.readline()
                
                while line and current_pair:
                    # the file pointer is before the next item, so advance pointer
                    if self.compare_pair_names(self.line_to_pair(line), current_pair) < 0:
                        line = infile.readline()

                    # the file pointer is at the next item, so yield it and advance both
                    if self.compare_pair_names(self.line_to_pair(line), current_pair) == 0:
                        yield self.line_to_output(line)
                        current_pair = next(group, None)
                        line = infile.readline()

                    # the file pointer is beyond where the next item would be, 
                    # so the next item is not found and we move onto the next one
                    else:
                        current_pair = next(group, None)





