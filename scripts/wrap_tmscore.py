import os
import itertools
import tempfile
import subprocess
import name_pdb
from multiprocessing import Pool


def calculate_dimers_TM_score(dimer1, dimer2, usalign_exe):
    '''
    This function calculates the TMscore between two dimers 
    using USalign. The dimers are input as chain pairs with
    each chain in separate files.

    Note: since this function uses USalign "fast" mode, it
    only gives TMscores to the hundreths place (e.g. 0.91
    instead of 0.9123) because any more that 0.01 precision
    likely does not mean much.

    :param dimer1: tuple, two string paths to pdb chains that
                    form a protein dimer.
    :param dimer2: tuple, two string paths to pdb chains that
                    form a second protein dimer.
    :usalign_exe: str, path to the compiled USalign executable

    :returns: TM_2to1, TM_1to2
                TM_2to1: float, the TM score with reference to
                            dimer1, rounded to 0.01 precision
                TM_1to2: float, the TM score with reference to
                            dimer2, rounded to 0.01 precision
    '''
    
    with tempfile.TemporaryDirectory() as td:
        
        with open(os.path.join(td, 'd1.pdb'), 'w') as dimer1file:
            with open(dimer1[0], 'r') as d1c1:
                dimer1file.write(d1c1.read())
            dimer1file.write('\n')
            with open(dimer1[1], 'r') as d1c2:
                dimer1file.write(d1c2.read())

        with open(os.path.join(td, 'd2.pdb'), 'w') as dimer2file:
            with open(dimer2[0], 'r') as d2c1:
                dimer2file.write(d2c1.read())
            dimer2file.write('\n')
            with open(dimer2[1], 'r') as d2c2:
                dimer2file.write(d2c2.read())

        cmd = f'{usalign_exe} -fast -mm 1 -ter 1 -outfmt 2 {dimer1file.name} {dimer2file.name}'
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)

        if result.returncode != 0:
            raise RuntimeError(f'usalign execution failed for {dimer1file.name}, {dimer2file.name} with message {result.stderr}{result.stdout}')

        result_line = result.stdout.splitlines()[-1]
        TM_2to1, TM_1to2 = result_line.split()[2:4] 

        TM_2to1 = round(float(TM_2to1), 2)
        TM_1to2 = round(float(TM_1to2), 2)

        return TM_2to1, TM_1to2



def calculate_many_dimers_TM_score(dimer_pairs, usalign_exe, lib_path, cores=1):
    args = []
    for dp in dimer_pairs:
        dimer1_chain1, dimer1_chain2 = name_pdb.dimer2pdbs(dp[0], lib_path)
        dimer2_chain1, dimer2_chain2 = name_pdb.dimer2pdbs(dp[1], lib_path)
        dimer1_t = (dimer1_chain1, dimer1_chain2)
        dimer2_t = (dimer2_chain1, dimer2_chain2)
        args.append( (dimer1_t, dimer2_t, usalign_exe) )

    with Pool(processes=cores) as p:
        return p.starmap(calculate_dimers_TM_score, args)



