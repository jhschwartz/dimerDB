'''
expand_uniparc2others.py - functions to expand chains of the raw uniparc2others.yaml file, changing names from simple '1abc_A' to '1abc_a1_m1_cA' as listed in the rcsb index file. Additionally, drops the uniprot matches as they are no longer needed.

Input looks like:
    UP000000123:
        pdb:
            - 1abc_A
            - 2xyz_B
        uniprot:
            - P12345
            - O98765

Output looks like:
    UP000000123: [1abc_a1_m1_cA, 1abc_a1_m2_cA, 2xyz_a2_m4_cB]

Written by Jacob Schwartz (jaschwa@umich.edu) in February 2023.
Copyright Jacob Schwartz, developed for the Peter Freddolino Lab while employed at the University of Michigan.
https://freddolino-lab.med.umich.edu
'''

import yaml

