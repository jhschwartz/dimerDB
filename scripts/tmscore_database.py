'''
tmscore_database.py - functions to update calculated dimer-dimer tmscores to a 
                      tsv database, and additionally retreive these scores.
                      This operates efficiently to avoid a bottleneck in the
                      HomodimerDB pipeline. 

    How is data stored?
        TMscores are stored in files like "xx.tsv" in a folder we call "db_path".
        Each file like "xx.tsv" has 4 columns, without headers:
            [Dimer 1 Name]  [Dimer 2 Name]  [TMscore, D1 to D2] [TMscore, D2 to D1]
        The "xx" in "xx.tsv" is hash defined by the two middle characters of the
        PDB entry to which Dimer 1 belongs. For example, given Dimer 1 with name 
        "1aa1-a1-m1-cA_1aa1-a1-m1-cB", the hash would be "aa" and the data for any
        pair which has this dimer as Dimer 1 would be in "aa.tsv". Note that Dimer 1
        is always alphabetically first compared to Dimer 2.
        
        Some examples:
          - "1aa1-a1-m1-cA_1aa1-a1-m1-cB" and "1aa1-a1-m1-cC_1aa1-a1-m1-cD" ? 
              -> "1aa1-a1-m1-cA_1aa1-a1-m1-cB" is Dimer 1 and data is in "aa.tsv".

          - "2aa1-a1-m1-cA_2aa1-a1-m1-cB" and "1zz2-a1-m1-cA_1zz2-a1-m1-cB" ? 
              -> "1zz2-a1-m1-cA_1zz2-a1-m1-cB" is Dimer 1 and data is in "zz.tsv".

          - "3bb7-a1-m1-cA_3bb7-a1-m1-cB" and "5rr1-a7-m6-c7_5rr1-a7-m14-c7" ?
              -> "3bb7-a1-m1-cA_3bb7-a1-m1-cB" is Dimer 1 and data is in "bb.tsv"

'''
import os
import itertools
import tempfile
#from name_pdb import get_div
get_div = lambda name: name[1:3]



def update_db(new_pairs_scores, db_path):
    '''
    This function updates the database to insert new scores into .tsv files in
    their appropriate place to maintain sorting. If a .tsv file of interest does
    not exist, this will create it. So, this function can be used to create and
    initialize a database, or update an already existing database with new scores.

    Note also that this can change scores of pairs already existing in the database.

    :param new_pairs_scores: list of tuples to describe pairs and their tmscores
                             to insert/update.
                              - each tuple is of size 4: 
                                    (dimer1name, dimer2name, score1to2, score2to1)
                              - each tuple element is a string
                              - dimer1name is alphabetically before dimer2name
                              - dimer names are each named with chains sorted
    :param db_path: str, the path to the tmscore database already existing or to
                    be created.
    '''
    div_dimer1 = lambda item: get_div(item[0])
    new_data = sorted(new_pairs_scores, key=lambda t: (div_dimer1(t), *t))

    os.makedirs(db_path, exist_ok=True)

    for div, group in itertools.groupby(new_data, key=div_dimer1):
        
        db_file = os.path.join(db_path, div+'.tsv')
        
        if not os.path.exists(db_file):
            with open(db_file, 'w') as f:
                for pair_and_score in group:
                    f.write('{}\t{}\t{}\t{}\n'.format(*pair_and_score))
        
        else:
            with tempfile.NamedTemporaryFile('w+t') as ftmp:
                
                with open(db_file, 'r') as fxx:
                    
                    line = fxx.readline()
                    new_pair = next(group, None)
                    
                    while new_pair and line:
                        
                        new_d1, new_d2 = new_pair[:2]
                        old_d1, old_d2 = line.split()[:2]
                        
                        # new pair is alphabetically first, so write it and advance new
                        if (new_d1, new_d2) < (old_d1, old_d2):
                            ftmp.write('{}\t{}\t{}\t{}\n'.format(*new_pair))
                            new_pair = next(group, None)

                        # new and old pairs are same, so pick new and skip old
                        elif (new_d1, new_d2) == (old_d1, old_d2):
                            ftmp.write('{}\t{}\t{}\t{}\n'.format(*new_pair))
                            new_pair = next(group, None)
                            line = fxx.readline()

                        # old pair is alphabetically first, so write it and advance old
                        else:
                            ftmp.write('{}\t{}\t{}\t{}\n'.format(*line.split()))
                            line = fxx.readline()

                    # write any leftovers from old pair
                    while line:
                        ftmp.write('{}\t{}\t{}\t{}\n'.format(*line.split()))
                        line = fxx.readline()

                    # write any leftovers from new_pair
                    while new_pair:
                        ftmp.write('{}\t{}\t{}\t{}\n'.format(*new_pair))
                        new_pair = next(group, None)

                # now that the updated "xx.tsv" is in a tempfile, overwrite the original
                with open(db_file, 'w') as fxx:
                    ftmp.seek(0)
                    for line in ftmp:
                        fxx.write(line)





def lookup_scores(dimer_pairs, db_path):
    '''
    This function retrieves the tm scores of the input dimer_pairs
    that exist in the database. Any that are nonexistent are simply 
    skipped over. The return values are not guranteed to match input order.

    :param dimer_pairs: list of which each element a tuple to describe pairs to retrive
                        - tuple size two, each element a dimer name as defined in name_pdb.py
                        - these describe dimer pairs which need tm scores between them
                        - it is expected that the tuples themselves are sorted, and that
                          the dimers themselves, as elements of the tuples, have their
                          two chains sorted to make up each dimer's name
    :param db_path: str, the path to the tmscore database which contains files like "xx.txt"

    :yields: indimer1, indimer2, in1to2score, in2to1score
        indimer1: str, name of dimer1 of a pair found in the database
        indimer2: str, name of dimer2 of a pair found in the database
        in1to2score: str, tmscore of indimer1 to indimer2 as found in the database
        in2to1score: str, tmscore of indimer2 to indimer1 as found in the database
    '''


    div_dimer1 = lambda item: get_div(item[0])
    data = sorted(dimer_pairs, key=lambda t: (div_dimer1(t), *t))
    
    for div, group in itertools.groupby(data, key=div_dimer1):
   
        db_file = os.path.join(db_path, div+'.tsv')
        if not os.path.exists(db_file):
            continue
   
        with open(db_file, 'r') as infile:

            current_pair = next(group, None)
            line = infile.readline()
            
            while line and current_pair:

                indimer1, indimer2, in1to2score, in2to1score = line.split()

                # the file pointer is before the next item, so advance pointer
                if indimer1 < current_pair[0]:
                    line = infile.readline()

                # the file pointer is at the next item, so yield it and advance both
                elif (indimer1, indimer2) == current_pair:
                    yield indimer1, indimer2, in1to2score, in2to1score 
                    current_pair = next(group, None)
                    line = infile.readline()

                # the file pointer is beyond where the next item would be, 
                # so the next item is not found and we move onto the next one
                else:
                    current_pair = next(group, None)
            

